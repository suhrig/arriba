#include <iostream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <unistd.h>
#include <algorithm>
#include "options.hpp"
#include "options_extract_read-through_fusions.hpp"

using namespace std;

options_t get_default_options() {
	options_t options;

	options.input_bam_file = "/dev/stdin";
	options.output_bam_file = "/dev/stdout";

	return options;
}

void print_usage(const string& error_message) {
	if (error_message != "")
		cerr << "ERROR: " << error_message << endl;

	options_t default_options = get_default_options();

	cerr << endl
	     << "Ariba RNA fusion detector - extract_read-through_fusions" << endl
	     << "---------------------------------------------------------" << endl
	     << "Version: " << ARIBA_VERSION << endl << endl
	     << "This is a helper utility of Ariba The STAR RNA-Seq aligner does " << endl
	     << "not report read-through fusions in the chimeric BAM file. This program " << endl
	     << "extracts reads supporting read-through fusions from the RNA BAM file. " << endl
	     << "The output file should be passed to ariba via the parameter -r." << endl
	     << "For optimal performance extract_read-through_fusions should be run " << endl
	     << "while STAR is running (see usage)." << endl << endl
	     << "Usage: extract_read-through_fusions -g genes.bed -i rna.bam -o read_through.bam" << endl
	     << "Usage: STAR --outStd BAM [...] | tee rna.bam | extract_read-through_fusions -g genes.bed > read_through.bam" << endl << endl
	     << wrap_help("-g FILE", "BED file with gene annotation. The following columns are "
	                  "required: (1) contig, (2) gene_start, (3) gene_end, "
	                  "(4) gene_name, (5) ignored, (6) strand.")
	     << wrap_help("-i FILE", "Input file in BAM format containing alignments from STAR. "
	                  "The file need not be sorted. Default: " + default_options.input_bam_file)
	     << wrap_help("-o FILE", "Output file in BAM format containing reads which support "
	                  "read-through fusions. Default: " + default_options.output_bam_file)
	     << wrap_help("-h", "Print help and exit.")
	     << "Questions or problems may be sent to: " << HELP_CONTACT << endl;
	exit(1);
}

options_t parse_arguments(int argc, char **argv) {
	options_t options = get_default_options();
	istringstream disabled_filters;

	// parse arguments
	opterr = 0;
	int c;
	while ((c = getopt(argc, argv, "g:i:o:h")) != -1) {
		switch (c) {
			case 'g':
				options.gene_annotation_file = optarg;
				if (access(options.gene_annotation_file.c_str(), R_OK) != 0) {
					cerr << "ERROR: File '" << options.gene_annotation_file << "' not found.";
					exit(1);
				}
				break;
			case 'i':
				options.input_bam_file = optarg;
				if (access(options.input_bam_file.c_str(), R_OK) != 0) {
					cerr << "ERROR: File '" << options.input_bam_file << "' not found.";
					exit(1);
				}
				break;
			case 'o':
				options.output_bam_file = optarg;
				break;
			case '?':
				switch (optopt) {
					case 'g': case 'i': case 'o':
						print_usage(string("Option -") + ((char) optopt) + " requires an argument.");
						break;
					default:
						print_usage(string("Unknown option: -") + ((char) optopt));
						break;
				}
			case 'h':
			default:
				print_usage();
		}
	}

	// check for mandatory arguments
	if (options.gene_annotation_file.empty())
		print_usage("Missing mandatory option: -g");
	if (options.input_bam_file.empty())
		print_usage("Missing mandatory option: -i");
	if (options.output_bam_file.empty())
		print_usage("Missing mandatory option: -o");

	return options;
}

