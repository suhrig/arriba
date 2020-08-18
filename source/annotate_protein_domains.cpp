#include <string>
#include <unordered_map>
#include "common.hpp"
#include "annotation.hpp"
#include "read_compressed_file.hpp"
#include "annotate_protein_domains.hpp"

using namespace std;

bool get_gff3_attribute(const string& attributes, const string& attribute_name, string& attribute_value) {

	// find start of attribute
	size_t start = attributes.find(attribute_name + "=");
	if (start == string::npos) {
		cerr << "WARNING: failed to extract " << attribute_name << " from line in GFF3 file: " << attributes << endl;
		return false;
	}
	start += attribute_name.size() + 1; // move to position after "="

	// find end of attribute
	size_t end = attributes.find(';', start);
	if (end == string::npos)
		end = attributes.size();

	attribute_value = attributes.substr(start, end - start);
	return true;
}

void load_protein_domains(const string& filename, const contigs_t& contigs, const gene_annotation_t& gene_annotation, const unordered_map<string,gene_t>& gene_names, protein_domain_annotation_t& protein_domain_annotation, protein_domain_annotation_index_t& protein_domain_annotation_index) {

	// make a map of gene id -> gene
	unordered_map<string,gene_t> gene_ids;
	for (auto gene = gene_annotation.begin(); gene != gene_annotation.end(); ++gene)
		gene_ids[gene->gene_id] = (gene_t) &(*gene);

	// read lines from GFF3 file
	autodecompress_file_t gff3_file(filename);
	string line;
	while (gff3_file.getline(line)) {
		if (!line.empty() && line[0] != '#') { // skip comment lines

			tsv_stream_t tsv(line);
			protein_domain_annotation_record_t protein_domain;
			string contig, strand, attributes, gene_name, gene_id, trash;

			// parse line
			tsv >> contig >> trash >> trash >> protein_domain.start >> protein_domain.end >> trash >> strand >> trash >> attributes;
			if (tsv.fail() || contig.empty() || strand.empty() || attributes.empty()) {
				cerr << "WARNING: failed to parse line in GFF3 file: " << line << endl;
				continue;
			}

			// extract gene name and protein domain name from attributes
			if (!get_gff3_attribute(attributes, "gene_name", gene_name) ||
			    !get_gff3_attribute(attributes, "gene_id", gene_id) ||
			    !get_gff3_attribute(attributes, "Name", protein_domain.name))
				continue;

			// convert string representation of contig to numeric ID
			contigs_t::const_iterator find_contig_by_name = contigs.find(removeChr(contig));
			if (find_contig_by_name == contigs.end()) {
				cerr << "WARNING: unknown contig: " << contig << endl;
				continue;
			}

			// decode hex representations of special characters (e.g., "%2C" for ",")
			for (string::size_type pos = protein_domain.name.find("%"); pos != string::npos; pos = protein_domain.name.find("%", pos+1)) {
				if (pos + 2 < protein_domain.name.size()) {
					char c1 = protein_domain.name[pos+1];
					char c2 = protein_domain.name[pos+2];
					if ((c1 >= '0' && c1 <= '9' || c1 >= 'a' && c1 <= 'f' || c1 >= 'A' && c1 <= 'F') &&
					    (c2 >= '0' && c2 <= '9' || c2 >= 'a' && c2 <= 'f' || c2 >= 'A' && c2 <= 'F')) {
						istringstream iss(protein_domain.name.substr(pos+1, 2));
						unsigned int c;
						iss >> hex >> c;
						protein_domain.name = protein_domain.name.substr(0, pos) + ((char) c) + protein_domain.name.substr(pos + 3);
					}
				}
			}
			// replace white space, non-printable characters, commas and pipes with underscores
			// (the latter two because they have a special meaning in the output format)
			for (string::size_type pos = 0; pos < protein_domain.name.size(); ++pos)
				if (protein_domain.name[pos] < '!' || protein_domain.name[pos] > '~' || protein_domain.name[pos] == ',' || protein_domain.name[pos] == '|')
					protein_domain.name[pos] = '_';

			// map protein domain to gene by gene ID or name
			string short_gene_id; // if Gencode, remove version number from gene ID
                        string::size_type trim_position = string::npos;
			if (gene_id.substr(0, 3) == "ENS" && (trim_position = gene_id.find_last_of('.', string::npos)) != string::npos)
				short_gene_id = gene_id.substr(0, trim_position);
			else
				short_gene_id = gene_id;
			auto find_gene_by_id = gene_ids.find(gene_id);
			if (find_gene_by_id == gene_ids.end()) {
				auto find_gene_by_name = gene_names.find(gene_name);
				if (find_gene_by_name == gene_names.end()) {
					cerr << "WARNING: unknown gene: " << gene_name << " " << gene_id << endl;
					continue;
				} else {
					protein_domain.gene = find_gene_by_name->second;
				}
			} else {
				protein_domain.gene = find_gene_by_id->second;
			}

			// make annotation record
			protein_domain.contig = find_contig_by_name->second;
			protein_domain.start--; // GFF3 files are one-based
			protein_domain.end--; // GFF3 files are one-based
			protein_domain.strand = (strand[0] == '+') ? FORWARD : REVERSE;
			protein_domain_annotation.push_back(protein_domain);
		}
	}

	if (protein_domain_annotation.empty()) {
		cerr << "ERROR: failed to parse GFF3 file" << endl;
		exit(1);
	}

	// index domains by coordinate
	make_annotation_index(protein_domain_annotation, protein_domain_annotation_index);
}


string annotate_retained_protein_domains(const contig_t contig, const position_t breakpoint, const strand_t predicted_strand, const bool predicted_strand_ambiguous, const gene_t gene, const direction_t direction, const protein_domain_annotation_index_t& protein_domain_annotation_index) {

	if (!gene->is_protein_coding)
		return "";

	if (predicted_strand_ambiguous || predicted_strand != gene->strand)
		return "";

	if (contig >= protein_domain_annotation_index.size())
		return "";

	// determine part of gene that is retained in the fusion
	position_t start = (direction == UPSTREAM) ? breakpoint : gene->start;
	position_t end = (direction == UPSTREAM) ? gene->end : breakpoint;

	// find all protein domains located in the retained part of gene
	set<string> retained_protein_domains;
	for (auto protein_domains = protein_domain_annotation_index[contig].lower_bound(start); protein_domains != protein_domain_annotation_index[contig].end() && protein_domains->first <= end; ++protein_domains)
		for (auto protein_domain = protein_domains->second.begin(); protein_domain != protein_domains->second.end(); ++protein_domain)
			if ((**protein_domain).gene == gene)
				retained_protein_domains.insert((**protein_domain).name);

	// concatenate tags to comma-separated string
	string result;
	for (auto protein_domain = retained_protein_domains.begin(); protein_domain != retained_protein_domains.end(); ++protein_domain) {
		if (!result.empty())
			result += ",";
		result += *protein_domain;
	}
	return result;
}

