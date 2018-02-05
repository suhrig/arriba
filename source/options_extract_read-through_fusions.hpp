#ifndef _H_OPTIONS_EXTRACT_READ_THROUGH_FUSIONS_H
#define _H_OPTIONS_EXTRACT_READ_THROUGH_FUSIONS_H 1

#include <string>

using namespace std;

struct options_t {
	string gene_annotation_file;
	string input_bam_file;
	string output_bam_file;
	string gtf_features;
};

options_t parse_arguments(int argc, char **argv);

#endif /* _H_OPTIONS_EXTRACT_READ_THROUGH_FUSIONS_H */
