#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <unordered_map>
#include "sam.h"
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

bool get_gtf_attribute(const string& attributes, const string& attribute_name, string& attribute_value) {

	// find start of attribute
	size_t start = attributes.find(attribute_name + " \"");
	if (start != string::npos)
		start = attributes.find('"', start);
	if (start == string::npos)
		return false;
	start++;

	// find end of attribute
	size_t end = attributes.find('"', start);
	if (end == string::npos)
		return false;
	attribute_value = attributes.substr(start, end - start);

	return true;
}

void read_annotation_gtf(const string& filename, const contigs_t& contigs, gene_annotation_t& gene_annotation, exon_annotation_t& exon_annotation, unordered_map<string,gene_t>& gene_names) {

	unordered_map<string,transcript_t> transcripts; // translates transcript IDs to numeric IDs

	unordered_map<exon_t,string> gene_name_by_exon; // remembers which exon belongs to which gene

	stringstream gtf_file;
	autodecompress_file(filename, gtf_file);
	string line;
	while (getline(gtf_file, line)) {
		if (!line.empty() && line[0] != '#') { // skip comment lines

			istringstream iss(line);
			annotation_record_t annotation_record;
			string contig, strand, feature, attributes, trash, gene_name;

			// parse line
			iss >> contig >> trash >> feature >> annotation_record.start >> annotation_record.end >> trash >> strand >> trash;
			if (contig.empty() || feature.empty() || strand.empty()) {
				cerr << "WARNING: failed to parse line in GTF file '" << filename << "': " << line << endl;
				continue;
			}
			getline(iss, attributes);

			// extract gene name from attributes
			if (!get_gtf_attribute(attributes, "gene_name", gene_name)) {
				cerr << "WARNING: failed to extract gene name from line in GTF file '" << filename << "': " << line << endl;
				continue;
			}

			contig = removeChr(contig);
			if (contigs.find(contig) == contigs.end()) {
				cerr << "WARNING: unknown contig in GTF file '" << filename << "': " << contig << endl;

			} else {

				// make annotation record
				annotation_record.contig = contigs.at(contig);
				annotation_record.start--; // GTF files are one-based
				annotation_record.end--; // GTF files are one-based
				annotation_record.strand = (strand[0] == '+') ? FORWARD : REVERSE;

				if (feature == "gene") {

					// make gene annotation record
					gene_annotation_record_t gene_annotation_record;
					gene_annotation_record.copy(annotation_record);

					gene_annotation_record.name = gene_name;
					gene_annotation_record.exonic_length = 0; // is calculated later
					gene_annotation_record.is_dummy = false;
					gene_annotation.push_back(gene_annotation_record);

				} else if (feature == "exon" || feature == "UTR") {

					// make exon annotation record
					exon_annotation_record_t exon_annotation_record;
					exon_annotation_record.copy(annotation_record);

					// extract transcript ID from attributes
					string transcript_id;
					if (!get_gtf_attribute(attributes, "transcript_id", transcript_id)) {
						cerr << "WARNING: failed to extract transcript ID from line in GTF file '" << filename << "': " << line << endl;
						continue;
					}
					exon_annotation_record.transcript = transcripts[transcript_id];
					if (exon_annotation_record.transcript == 0) // this is the first time we encounter this transcript ID => give it a numeric ID
						exon_annotation_record.transcript = transcripts[transcript_id] = transcripts.size();

					// give the record an ID (= a consecutive number)
					exon_annotation.push_back(exon_annotation_record);

					// remember name of gene, so we can make pointers later
					gene_name_by_exon[&(exon_annotation[exon_annotation.size()-1])] = gene_name;
				}
			}
		}
	}

	// assign exons to genes
	for (gene_annotation_t::iterator gene = gene_annotation.begin(); gene != gene_annotation.end(); ++gene)
		gene_names[gene->name] = &(*gene);
	for (exon_annotation_t::iterator exon = exon_annotation.begin(); exon != exon_annotation.end(); ++exon) {
		auto gene_name = gene_names.find(gene_name_by_exon[&(*exon)]);
		if (gene_name == gene_names.end()) {
			cerr << "ERROR: exon belongs to unknown gene: " << gene_name_by_exon[&(*exon)] << endl;
			exit(1);
		} else {
			exon->gene = gene_name->second;
		}
	}

}

gene_multiset_t get_genes_from_exons(const exon_multiset_t& exons) {
	gene_multiset_t result;
	for (exon_multiset_t::const_iterator exon = exons.begin(); exon != exons.end(); ++exon)
		result.insert((**exon).gene);
	return result;
}

// check if a breakpoint is near an annotated splice site
exon_multiset_t get_exons_from_splice_site(const gene_t gene, const direction_t direction, const contig_t contig, const position_t breakpoint, const exon_annotation_index_t& exon_annotation_index) {

	const unsigned int max_splice_site_distance = 2;

	// find overlapping exons
	exon_contig_annotation_index_t::const_iterator exons_at_breakpoint = exon_annotation_index[contig].lower_bound(breakpoint);
	exon_contig_annotation_index_t::const_iterator exons_before_breakpoint = exons_at_breakpoint;
	if (exons_before_breakpoint != exon_annotation_index[contig].begin() && exon_annotation_index[contig].size() > 0)
		--exons_before_breakpoint;

	// count how many of the overlapping exons belong to <gene>
	unsigned int genes_at_breakpoint = (exons_at_breakpoint != exon_annotation_index[contig].end()) ? get_genes_from_exons(exons_at_breakpoint->second).count(gene) : 0;
	unsigned int genes_before_breakpoint = (exons_before_breakpoint != exon_annotation_index[contig].end()) ? get_genes_from_exons(exons_before_breakpoint->second).count(gene) : 0;

	if (exons_before_breakpoint != exon_annotation_index[contig].end() && breakpoint - exons_before_breakpoint->first <= max_splice_site_distance) {
		if (direction == UPSTREAM && exons_at_breakpoint != exon_annotation_index[contig].end() && genes_at_breakpoint > genes_before_breakpoint)
			return exons_at_breakpoint->second;

		if (direction == DOWNSTREAM && (exons_at_breakpoint == exon_annotation_index[contig].end() && genes_before_breakpoint > 0 || exons_at_breakpoint != exon_annotation_index[contig].end() && genes_before_breakpoint > genes_at_breakpoint))
			return exons_before_breakpoint->second;
	}

	exon_contig_annotation_index_t::const_iterator exons_after_breakpoint = exons_at_breakpoint;
	if (exons_after_breakpoint != exon_annotation_index[contig].end())
		++exons_after_breakpoint;

	unsigned int genes_after_breakpoint = (exons_after_breakpoint != exon_annotation_index[contig].end()) ? get_genes_from_exons(exons_after_breakpoint->second).count(gene) : 0;

	if (exons_at_breakpoint != exon_annotation_index[contig].end() && exons_at_breakpoint->first - breakpoint <= max_splice_site_distance) {
		if (direction == UPSTREAM && exons_after_breakpoint != exon_annotation_index[contig].end() && genes_after_breakpoint > genes_at_breakpoint)
			return exons_after_breakpoint->second;

		if (direction == DOWNSTREAM && (exons_after_breakpoint == exon_annotation_index[contig].end() && genes_at_breakpoint > 0 || exons_after_breakpoint != exon_annotation_index[contig].end() && genes_at_breakpoint > genes_after_breakpoint))
			return exons_at_breakpoint->second;
	}

	exon_multiset_t empty_set;
	return empty_set;
}

// check if a breakpoint is near an annotated splice site
bool is_breakpoint_spliced(const gene_t gene, const direction_t direction, const contig_t contig, const position_t breakpoint, const exon_annotation_index_t& exon_annotation_index) {
	return get_exons_from_splice_site(gene, direction, contig, breakpoint, exon_annotation_index).size() > 0;
}

void get_annotation_by_alignment(const alignment_t& alignment, gene_set_t& gene_set, const exon_annotation_index_t& exon_annotation_index) {

	// first, try to annotate based on the boundaries (start+end) of the alignment
	exon_set_t exon_set;
	get_annotation_by_coordinate(alignment.contig, alignment.start, alignment.end, exon_set, exon_annotation_index);

	// translate exons to genes
	for (auto exon = exon_set.begin(); exon != exon_set.end(); ++exon)
		gene_set.insert((**exon).gene);

	// when the alignment overlaps with multiple genes, try to resolve the ambiguity by looking for splice sites that are specific for one gene
	if (gene_set.size() > 1) {
		// cycle through CIGAR string and look for introns
		gene_set_t gene_set_supported_by_splicing = gene_set;
		position_t reference_position = alignment.start;
		for (unsigned int i = 0; i < alignment.cigar.size() && gene_set_supported_by_splicing.size() > 1; ++i) {
			switch (alignment.cigar.operation(i)) {
				case BAM_CREF_SKIP:

					// whenever we hit an intron in the CIGAR string, check if it aligns with splice sites for each gene
					for (gene_set_t::iterator gene = gene_set_supported_by_splicing.begin(); gene != gene_set_supported_by_splicing.end();) {
						if (!is_breakpoint_spliced(*gene, DOWNSTREAM, alignment.contig, reference_position, exon_annotation_index) &&
						    !is_breakpoint_spliced(*gene, UPSTREAM, alignment.contig, reference_position + alignment.cigar.op_length(i), exon_annotation_index))
							gene = gene_set_supported_by_splicing.erase(gene);
						else
							++gene;
					}

					// none of the genes match the splice pattern of the alignment
					// => it's a tie, so we hope to find another intron in the CIGAR string
					if (gene_set_supported_by_splicing.empty())
						gene_set_supported_by_splicing = gene_set;

					// fall-through to next case
				case BAM_CMATCH:
				case BAM_CDEL:
					reference_position += alignment.cigar.op_length(i);
					break;
			}
		}

		// if none of the genes match the splice pattern of the alignment, return all genes,
		// otherwise, only return the matching genes
		if (!gene_set_supported_by_splicing.empty() && gene_set_supported_by_splicing.size() < gene_set.size())
			gene_set = gene_set_supported_by_splicing;
	}
}

// when a read overlaps with multiple genes, this function returns the boundaries of the biggest one
void get_boundaries_of_biggest_gene(gene_set_t& genes, position_t& start, position_t& end) {
	start = -1;
	end = -1;
	for (gene_set_t::iterator i = genes.begin(); i != genes.end(); ++i) {
		if (start == -1 || start > (**i).start)
			start = (**i).start;
		if (end == -1 || end < (**i).end)
			end = (**i).end;
	}
}

void fetch_gene_sequences_from_fasta(const string& assembly_file_path, fusions_t& fusions, const vector<string>& contigs_by_id) {

	// open FastA file and index
	faidx_t* genome_fa_index = fai_load(assembly_file_path.c_str());

	// fetch sequences of non-discarded fusions from FastA file
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (!i->second.filters.empty())
			continue; // skip discarded fusions

		vector<gene_t> genes = { i->second.gene1, i->second.gene2 };
		for (vector<gene_t>::iterator gene = genes.begin(); gene != genes.end(); ++gene) {

			if ((**gene).sequence.empty()) { // if we have not already fetched the sequence
				char* sequence;
				int length;
				string region_string = contigs_by_id[(**gene).contig] + ":" + to_string((**gene).start+1) + "-" + to_string((**gene).end+1);
				sequence = fai_fetch(genome_fa_index, region_string.c_str(), &length);
				if (length == 0) { // fetching failed
					// maybe the FastA file has "chr" in contig names => try again
					sequence = fai_fetch(genome_fa_index, addChr(region_string).c_str(), &length);
					if (length == 0) { // fetching failed again => give up
						cerr << "WARNING: Failed to fetch sequence of feature '" << (**gene).name << "' from '" << assembly_file_path << "'.";
						return;
					}
				}
				(**gene).sequence = string(sequence, sequence + length);
				// convert to uppercase
				std::transform((**gene).sequence.begin(), (**gene).sequence.end(), (**gene).sequence.begin(), (int (*)(int))std::toupper);
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


void get_unique_transcripts_from_exons(const exon_multiset_t& exons, set<transcript_t>& transcripts) {
	transcripts.clear();
	for (exon_multiset_t::const_iterator exon = exons.begin(); exon != exons.end(); ++exon)
		transcripts.insert((**exon).transcript);
}


// get the distance between two positions after splicing (i.e., ignoring introns)
int get_spliced_distance(const contig_t contig, position_t position1, position_t position2, direction_t direction1, direction_t direction2, const gene_t gene, const exon_annotation_index_t& exon_annotation_index) {

	if (position1 > position2) {
		swap(position1, position2);
		swap(direction1, direction2);
	}

	// find exon/intron of position1
	exon_contig_annotation_index_t::const_iterator p1 = exon_annotation_index[contig].lower_bound(position1);
	// find exon/intron of position2
	exon_contig_annotation_index_t::const_iterator p2 = exon_annotation_index[contig].lower_bound(position2);

	if (p1 == p2) // position1 and position2 are in the same intron/exon => distance = difference
		return position2 - position1;

	// the exon/intron where position1/position2 is located must only be counted partially
	// (i.e., from position1/position2 to next exon boundary)
	position_t distance = 0;
	exon_multiset_t exons_at_position1, exons_at_position2;
	exons_at_position1 = get_exons_from_splice_site(gene, direction1, contig, position1, exon_annotation_index);
	if (exons_at_position1.size() == 0) { // position is not at a splice site
		if (p1 != exon_annotation_index[contig].end())
			exons_at_position1 = p1->second;
		// add distance from position1 to next higher exon boundary
		distance += p1->first - position1;
	}
	exons_at_position2 = get_exons_from_splice_site(gene, direction2, contig, position2, exon_annotation_index);
	if (exons_at_position2.size() > 0) {
		--p2;
	} else {
		if (p2 != exon_annotation_index[contig].end())
			exons_at_position2 = p2->second;
		// add distance from position2 to next lower exon boundary
		distance += position2;
		--p2;
		distance -= p2->first;
	}

	// check if the positions are in exons which belong to the same transcript
	// if so, calculate the distance considering only the exons of this transcript
	// if not, calculate the distance considering all the exons between position1 and 2
	set<transcript_t> transcripts1, transcripts2;
	get_unique_transcripts_from_exons(exons_at_position1, transcripts1);
	get_unique_transcripts_from_exons(exons_at_position2, transcripts2);
	set<transcript_t> common_transcripts;
	set_intersection(transcripts1.begin(), transcripts1.end(), transcripts2.begin(), transcripts2.end(), inserter(common_transcripts, common_transcripts.begin()));
	unordered_map<transcript_t,int> distance_by_transcript;
	for (auto common_transcript = common_transcripts.begin(); common_transcript != common_transcripts.end(); ++common_transcript)
		distance_by_transcript[*common_transcript] = distance;

	// sum up exon sizes between position1 and position2
	while (p1 != p2) {
		position_t boundary = p1->first;
		++p1;
		if (get_genes_from_exons(p1->second).count(gene) > 0) {
			// calculate distance considering all exons
			distance += p1->first - boundary;
			// calculate distance considering only exons of transcripts common to position1 and 2
			set<transcript_t> transcripts_at_current_position;
			get_unique_transcripts_from_exons(p1->second, transcripts_at_current_position);
			for (auto transcript = transcripts_at_current_position.begin(); transcript != transcripts_at_current_position.end(); ++transcript) {
				auto d = distance_by_transcript.find(*transcript);
				if (d != distance_by_transcript.end())
					d->second += p1->first - boundary;
			}
		}
	}

	// select transcript with shortest distance
	for (auto d = distance_by_transcript.begin(); d != distance_by_transcript.end(); ++d)
		if (d->second < distance)
			distance = d->second;

	return distance;
}
