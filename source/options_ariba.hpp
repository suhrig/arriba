#ifndef _OPTIONS_ARIBA_H
#define _OPTIONS_ARIBA_H 1

#include <iostream>
#include <unordered_map>
#include <string>
#include <sstream>

using namespace std;

struct options_t {
	string chimeric_bam_file;
	string read_through_bam_file;
	string rna_bam_file;
	string genomic_breakpoints_file;
	unsigned int max_genomic_breakpoint_distance;
	string gene_annotation_file;
	string exon_annotation_file;
	string known_fusions_file;
	string output_file;
	string discarded_output_file;
	string assembly_file;
	string blacklist_file;
	string interesting_contigs;
	unsigned int homopolymer_length;
	unsigned int min_read_through_distance;
	unordered_map<string,bool> filters;
	float evalue_cutoff;
	unsigned int min_support;
	float max_mismapper_fraction;
	unsigned int min_anchor_length;
	bool print_supporting_reads;
	bool print_supporting_reads_for_discarded_fusions;
	bool print_fusion_sequence;
	bool print_fusion_sequence_for_discarded_fusions;
	bool single_end;
	float max_kmer_content;
	unsigned int fragment_length;
	string gtf_features;
};

options_t get_default_options();

void print_usage(const string& error_message = "");

options_t parse_arguments(int argc, char **argv);

#endif /* _OPTIONS_ARIBA_H */
