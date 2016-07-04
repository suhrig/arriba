#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <iostream>
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

void read_annotation_bed(const string& filename, annotation_t& annotation, contigs_t& contigs) {
	stringstream bed_file;
	autodecompress_file(filename, bed_file);
	string line;
	while (getline(bed_file, line)) {
		if (!line.empty() && line[0] != '#') { // skip comment lines

			istringstream iss(line);
			annotation_record_t annotation_record;
			string contig, strand, trash;

			// parse line
			iss >> contig >> annotation_record.start >> annotation_record.end >> annotation_record.name >> trash >> strand;
			if (contig.empty() || annotation_record.name.empty() || strand.empty()) {
				cerr << "WARNING: failed to parse line in BED file '" << filename << "': " << line << endl;
				continue;
			}

			contig = removeChr(contig);

			if (contigs.find(contig) == contigs.end()) {
				cerr << "WARNING: unknown contig in BED file '" << filename << "': " << contig << endl;
			} else {
				annotation_record.contig = contigs[contig];
				annotation_record.end--; // BED files are half open
				annotation_record.id = annotation.size();
				annotation_record.strand = (strand[0] == '+') ? FORWARD : REVERSE;
				annotation_record.exonic_length = 0; // is calculated later
				annotation_record.is_dummy = false;
				annotation.push_back(annotation_record);
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
void make_annotation_index(annotation_t annotation, annotation_index_t& annotation_index, const contigs_t& contigs) {

	annotation_index.resize(contigs.size()); // create a contig_annotation_index_t for each contig

	sort(annotation.rbegin(), annotation.rend()); // sort annotation by coordinate for traversal from end to beginning

	contig_t current_contig = annotation[0].contig; // keeps track of the contig we are currently processing
	position_t current_end = annotation[0].end; // keeps track of the end of the region we are currently processing

	struct overlapping_gene_t {
		gene_t gene;
		position_t start;
		// function to sort genes by start coordinate
		inline bool operator < (const overlapping_gene_t& x) const {
			return start < x.start;
		}
	};
	vector<overlapping_gene_t> overlapping_genes; // keeps track of the genes overlapping the current region and their start coordinates

	gene_set_t gene_set; // holds the gene_set_t that is going to be added to the annotation_index next

	for (unsigned int i = 0; i <= annotation.size(); ++i) { // traverse genome from end to beginning
		if (i != 0) { // skip first entry, because there is nothing to do yet

			// sort overlapping_genes so that higher coordinates are at the end of the vector
			// and can be removed by simply shrinking the vector
			sort(overlapping_genes.begin(), overlapping_genes.end());

			// for each combination of overlapping genes, add a gene_set_t to annotation_index
			int j;
			for (j = overlapping_genes.size()-1; j >= 0; j--) {
				if (i == annotation.size() || overlapping_genes[j].start > annotation[i].end || current_contig != annotation[i].contig) {
					if (j == overlapping_genes.size()-1 || overlapping_genes[j].start != overlapping_genes[j+1].start) {
						// add new record to annotation_index
						annotation_index[current_contig][current_end] = gene_set;
						current_end = overlapping_genes[j].start-1;
					}
					gene_set.erase(overlapping_genes[j].gene);
				} else break; // the end of overlapping_gene[j] has not been reached yet
			}

			// discard genes that we have passed
			overlapping_genes.resize(j+1);
		}

		if (i < annotation.size()) { // in the last iteration, i == annotation.size(); this iteration is needed to process the last record
			if (current_end != annotation[i].end || current_contig != annotation[i].contig) {
				// add new record to annotation index
				annotation_index[current_contig][current_end] = gene_set;
				current_end = annotation[i].end;
				current_contig = annotation[i].contig;
			}
			gene_set.insert(annotation[i].id);
			overlapping_gene_t overlapping_gene = { annotation[i].id, annotation[i].start };
			overlapping_genes.push_back(overlapping_gene);
		}
	}
	annotation_index[current_contig][current_end] = gene_set; // this adds an empty gene set at the end
}

void combine_annotations(gene_set_t& genes1, gene_set_t& genes2, gene_set_t& combined, bool make_union) {
	// when the two ends of a read map to different genes, the mapping is ambiguous
	// in this case, we try to resolve the ambiguity by taking the gene that both - start and end - overlap with
	set_intersection(genes1.begin(), genes1.end(), genes2.begin(), genes2.end(), inserter(combined, combined.begin()));
	if (combined.empty() && make_union)
		set_union(genes1.begin(), genes1.end(), genes2.begin(), genes2.end(), inserter(combined, combined.begin()));
}

void get_annotation_by_coordinate(const contig_t contig, const position_t start, const position_t end, gene_set_t& gene_set, annotation_index_t& annotation_index) {
//TODO support strand-specific libraries
	if (contig < annotation_index.size()) {
		contig_annotation_index_t::iterator result_start = annotation_index[contig].lower_bound(start);
		gene_set_t empty_set;
		if (start == end) {
			gene_set = (result_start != annotation_index[contig].end()) ? result_start->second : empty_set;
		} else {
			contig_annotation_index_t::iterator result_end = annotation_index[contig].lower_bound(end);
			combine_annotations(
				(result_start != annotation_index[contig].end()) ? result_start->second : empty_set,
				(result_end != annotation_index[contig].end()) ? result_end->second : empty_set,
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
		}
}

string dna_to_reverse_complement(string& dna) {
	string reverse_complement;
	dna_to_reverse_complement(dna, reverse_complement);
	return reverse_complement;
}

// check if a breakpoint is near an annotated splice site
bool is_breakpoint_spliced(const gene_t gene, const direction_t direction, const contig_t contig, const position_t breakpoint, annotation_index_t& exon_annotation_index) {
	gene_set_t genes_before_breakpoint;
	get_annotation_by_coordinate(contig, breakpoint-3, breakpoint-3, genes_before_breakpoint, exon_annotation_index);
	gene_set_t genes_after_breakpoint;
	get_annotation_by_coordinate(contig, breakpoint+3, breakpoint+3, genes_after_breakpoint, exon_annotation_index);
	if (direction == UPSTREAM) {
		return genes_before_breakpoint.find(gene) == genes_before_breakpoint.end() && genes_after_breakpoint.find(gene) != genes_after_breakpoint.end();
	} else { // direction == DOWNSTREAM
		return genes_before_breakpoint.find(gene) != genes_before_breakpoint.end() && genes_after_breakpoint.find(gene) == genes_after_breakpoint.end();
	}
}

