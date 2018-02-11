#ifndef _H_OPTIONS_EXTRACT_READS_H
#define _H_OPTIONS_EXTRACT_READS_H 1

#include <string>

using namespace std;

struct options_t {
	string gene_annotation_file;
	string rna_bam_file;
	string assembly_file;
	string read_through_file;
	string fastq_file_prefix;
	string chimeric_file;
	string gtf_features;
	unsigned int min_clipped_length;
};

options_t parse_arguments(int argc, char **argv);

#endif /* _H_OPTIONS_EXTRACT_READS_H */
