#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <unordered_map>
#include "annotation.hpp"
#include "read_compressed_file.hpp"
#include "htslib/faidx.h"

using namespace std;

string removeChr(string contig) {
	if (contig.substr(0, 3) == "chr")
		contig = contig.substr(3);
	if (contig == "M")
		contig = "MT";
	return contig;
}

string addChr(string contig) {
	if (contig.substr(0, 3) != "chr")
		contig = "chr" + contig;
	if (contig == "chrMT")
		contig = "chrM";
	return contig;
}

void read_annotation_gtf(const string& filename, annotation_t& gene_annotation, annotation_t& exon_annotation, contigs_t& contigs) {
	stringstream gtf_file;
	autodecompress_file(filename, gtf_file);
	string line;
	while (getline(gtf_file, line)) {
		if (!line.empty() && line[0] != '#') { // skip comment lines

			istringstream iss(line);
			annotation_record_t annotation_record;
			string contig, strand, feature, trash;

			// parse line
			iss >> contig >> trash >> feature >> annotation_record.start >> annotation_record.end >> trash >> strand >> trash;
			if (contig.empty() || feature.empty() || strand.empty()) {
				cerr << "WARNING: failed to parse line in GTF file '" << filename << "': " << line << endl;
				continue;
			}
			getline(iss, annotation_record.name);
			size_t gene_name_start = annotation_record.name.find("gene_name \"");
			if (gene_name_start != string::npos)
				gene_name_start = annotation_record.name.find('"', gene_name_start);
			if (gene_name_start == string::npos) {
				cerr << "WARNING: failed to extract gene name from line in GTF file '" << filename << "': " << line << endl;
				continue;
			}
			gene_name_start++;
			size_t gene_name_end = annotation_record.name.find('"', gene_name_start);
			if (gene_name_end == string::npos) {
				cerr << "WARNING: failed to extract gene name from line in GTF file '" << filename << "': " << line << endl;
				continue;
			}
			annotation_record.name = annotation_record.name.substr(gene_name_start, gene_name_end - gene_name_start);

			contig = removeChr(contig);
			if (contigs.find(contig) == contigs.end()) {
				cerr << "WARNING: unknown contig in GTF file '" << filename << "': " << contig << endl;

			} else {
				annotation_record.contig = contigs[contig];
				annotation_record.start--; // GTF files are one-based
				annotation_record.end--; // GTF files are one-based
				annotation_record.strand = (strand[0] == '+') ? FORWARD : REVERSE;
				annotation_record.exonic_length = 0; // is calculated later
				annotation_record.is_dummy = false;
				if (feature == "gene") {
					annotation_record.id = gene_annotation.size();
					gene_annotation.push_back(annotation_record);
				} else if (feature == "exon" || feature == "UTR") {
					annotation_record.id = exon_annotation.size();
					exon_annotation.push_back(annotation_record);
				}
			}
		}
	}
}

// split overlapping genes into disjunct regions
// for each region, collect the genes overlapping it and store them in a gene set
// for example, if the annotation contains two regions:
// - gene1 chr1:10,000-20,000
// - gene2 chr1:12,000-13,000
// then the resulting index will contain the following items:
// - chr1:10,000-11,999 gene1
// - chr1:12,000-13,000 gene1+gene2
// - chr1:13,001-20,000 gene1
void make_annotation_index(const annotation_t& annotation, annotation_index_t& annotation_index, const contigs_t& contigs) {
	annotation_index.resize(contigs.size()); // create a contig_annotation_index_t for each contig
	for (annotation_t::const_iterator feature = annotation.begin(); feature != annotation.end(); ++feature) {

		contig_annotation_index_t::const_iterator overlapping_genes = annotation_index[feature->contig].lower_bound(feature->end);
		if (overlapping_genes == annotation_index[feature->contig].end())
			annotation_index[feature->contig][feature->end]; // this creates an empty gene set, if it does not exist yet
		else
			annotation_index[feature->contig][feature->end] = overlapping_genes->second;

		overlapping_genes = annotation_index[feature->contig].lower_bound(feature->start-1);
		if (overlapping_genes == annotation_index[feature->contig].end())
			annotation_index[feature->contig][feature->start-1]; // this creates an empty gene set, if it does not exist yet
		else
			annotation_index[feature->contig][feature->start-1] = overlapping_genes->second;

		// add the gene to all gene sets between start and end of the gene
		for (contig_annotation_index_t::iterator i = annotation_index[feature->contig].lower_bound(feature->end); i->first >= feature->start; --i)
			i->second.insert(feature->id);
	}
}

void gene_multiset_to_set(gene_multiset_t gene_multiset, gene_set_t& gene_set) {
	for (gene_multiset_t::iterator i = gene_multiset.begin(); i != gene_multiset.end(); i = gene_multiset.upper_bound(*i))
		gene_set.insert(*i);
}

gene_set_t gene_multiset_to_set(gene_multiset_t gene_multiset) {
	gene_set_t gene_set;
	gene_multiset_to_set(gene_multiset, gene_set);
	return gene_set;
}

void combine_annotations(const gene_set_t& genes1, const gene_set_t& genes2, gene_set_t& combined, bool make_union) {
	// when the two ends of a read map to different genes, the mapping is ambiguous
	// in this case, we try to resolve the ambiguity by taking the gene that both - start and end - overlap with
	set_intersection(genes1.begin(), genes1.end(), genes2.begin(), genes2.end(), inserter(combined, combined.begin()));
	if (combined.empty() && make_union)
		set_union(genes1.begin(), genes1.end(), genes2.begin(), genes2.end(), inserter(combined, combined.begin()));
}

void get_annotation_by_coordinate(const contig_t contig, const position_t start, const position_t end, gene_set_t& gene_set, annotation_index_t& annotation_index) {
//TODO support strand-specific libraries
//TODO make use of splicing
	if (contig < annotation_index.size()) {
		contig_annotation_index_t::iterator result_start = annotation_index[contig].lower_bound(start);
		gene_set_t empty_set;
		if (start == end) {
			if (result_start != annotation_index[contig].end())
				gene_multiset_to_set(result_start->second, gene_set);
			else
				gene_set = empty_set;
		} else {
			contig_annotation_index_t::iterator result_end = annotation_index[contig].lower_bound(end);
			combine_annotations(
				(result_start != annotation_index[contig].end()) ? gene_multiset_to_set(result_start->second) : empty_set,
				(result_end != annotation_index[contig].end()) ? gene_multiset_to_set(result_end->second) : empty_set,
				gene_set
			);
		}
	}
}

// when a read overlaps with multiple genes, this function returns the boundaries of the biggest one
void get_boundaries_of_biggest_gene(gene_set_t& genes, const annotation_t& gene_annotation, position_t& start, position_t& end) {
	start = -1;
	end = -1;
	for (gene_set_t::iterator i = genes.begin(); i != genes.end(); ++i) {
		if (start == -1 || start > gene_annotation[*i].start)
			start = gene_annotation[*i].start;
		if (end == -1 || end < gene_annotation[*i].end)
			end = gene_annotation[*i].end;
	}
}

void fetch_gene_sequences_from_fasta(const string& assembly_file_path, fusions_t& fusions, annotation_t& gene_annotation, const vector<string>& contigs_by_id) {

	// open FastA file and index
	faidx_t* genome_fa_index = fai_load(assembly_file_path.c_str());

	// fetch sequences of non-discarded fusions from FastA file
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (!i->second.filters.empty())
			continue; // skip discarded fusions

		vector<gene_t> genes = { i->second.gene1, i->second.gene2 };
		for (vector<gene_t>::iterator j = genes.begin(); j != genes.end(); ++j) {

			if (gene_annotation[*j].sequence.empty()) { // if we have not already fetched the sequence
				char* sequence;
				int length;
				string region_string = contigs_by_id[gene_annotation[*j].contig] + ":" + to_string(gene_annotation[*j].start+1) + "-" + to_string(gene_annotation[*j].end+1);
				sequence = fai_fetch(genome_fa_index, region_string.c_str(), &length);
				if (length == 0) { // fetching failed
					// maybe the FastA file has "chr" in contig names => try again
					sequence = fai_fetch(genome_fa_index, addChr(region_string).c_str(), &length);
					if (length == 0) { // fetching failed again => give up
						cerr << "WARNING: Failed to fetch sequence of feature '" << gene_annotation[*j].name << "' from '" << assembly_file_path << "'.";
						return;
					}
				}
				gene_annotation[*j].sequence = string(sequence, sequence + length);
				// convert to uppercase
				std::transform(gene_annotation[*j].sequence.begin(), gene_annotation[*j].sequence.end(), gene_annotation[*j].sequence.begin(), (int (*)(int))std::toupper);
				free(sequence);
			}
		}
	}

	// close FastA file and index
	fai_destroy(genome_fa_index);
}

void dna_to_reverse_complement(string& dna, string& reverse_complement) {
	if (!reverse_complement.empty())
		reverse_complement.clear();
	reverse_complement.reserve(dna.length());
	for (string::reverse_iterator i = dna.rbegin(); i != dna.rend(); ++i)
		if (*i == 'a') {
			reverse_complement += 't';
		} else if (*i == 't') {
			reverse_complement += 'a';
		} else if (*i == 'c') {
			reverse_complement += 'g';
		} else if (*i == 'g') {
			reverse_complement += 'c';
		} else if (*i == 'A') {
			reverse_complement += 'T';
		} else if (*i == 'T') {
			reverse_complement += 'A';
		} else if (*i == 'C') {
			reverse_complement += 'G';
		} else if (*i == 'G') {
			reverse_complement += 'C';
		} else if (*i == '[') {
			reverse_complement += ']';
		} else if (*i == ']') {
			reverse_complement += '[';
		} else
			reverse_complement += *i;
}

string dna_to_reverse_complement(string& dna) {
	string reverse_complement;
	dna_to_reverse_complement(dna, reverse_complement);
	return reverse_complement;
}

// check if a breakpoint is near an annotated splice site
bool is_breakpoint_spliced(const gene_t gene, const direction_t direction, const contig_t contig, const position_t breakpoint, const annotation_index_t& exon_annotation_index) {

	const unsigned int max_splice_site_distance = 2;

	// find overlapping exons
	contig_annotation_index_t::const_iterator genes_at_breakpoint = exon_annotation_index[contig].lower_bound(breakpoint);
	contig_annotation_index_t::const_iterator genes_before_breakpoint = genes_at_breakpoint;
	if (genes_before_breakpoint != exon_annotation_index[contig].begin() && exon_annotation_index[contig].size() > 0)
		--genes_before_breakpoint;

	if (genes_before_breakpoint != exon_annotation_index[contig].end() && breakpoint - genes_before_breakpoint->first <= max_splice_site_distance &&
	    (direction == UPSTREAM && genes_at_breakpoint != exon_annotation_index[contig].end() && genes_at_breakpoint->second.count(gene) > genes_before_breakpoint->second.count(gene) ||
	     direction == DOWNSTREAM && (genes_at_breakpoint == exon_annotation_index[contig].end() && genes_before_breakpoint->second.count(gene) > 0 || genes_at_breakpoint != exon_annotation_index[contig].end() && genes_before_breakpoint->second.count(gene) > genes_at_breakpoint->second.count(gene))))
		return true;

	contig_annotation_index_t::const_iterator genes_after_breakpoint = genes_at_breakpoint;
	if (genes_after_breakpoint != exon_annotation_index[contig].end())
		++genes_after_breakpoint;

	if (genes_at_breakpoint != exon_annotation_index[contig].end() && genes_at_breakpoint->first - breakpoint <= max_splice_site_distance &&
	    (direction == UPSTREAM && genes_after_breakpoint != exon_annotation_index[contig].end() && genes_after_breakpoint->second.count(gene) > genes_at_breakpoint->second.count(gene) ||
	     direction == DOWNSTREAM && (genes_after_breakpoint == exon_annotation_index[contig].end() && genes_at_breakpoint->second.count(gene) > 0 || genes_after_breakpoint != exon_annotation_index[contig].end() && genes_at_breakpoint->second.count(gene) > genes_after_breakpoint->second.count(gene))))
		return true;

	return false;
}

