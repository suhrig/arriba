#include <iostream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <unistd.h>
#include <algorithm>
#include "annotation.hpp"
#include "options.hpp"
#include "options_extract_read-through_fusions.hpp"

using namespace std;

options_t get_default_options() {
	options_t options;

	options.input_bam_file = "/dev/stdin";
	options.output_bam_file = "/dev/stdout";
	options.gtf_features = DEFAULT_GTF_FEATURES;

	return options;
}

void print_usage(const string& error_message) {
	if (error_message != "")
		cerr << "ERROR: " << error_message << endl;

	options_t default_options = get_default_options();

	cerr << endl
	     << "Arriba gene fusion detector - extract_read-through_fusions" << endl
	     << "----------------------------------------------------------" << endl
	     << "Version: " << ARRIBA_VERSION << endl << endl
	     << "This is a helper utility of Arriba. The STAR RNA-Seq aligner does " << endl
	     << "not report read-through fusions in the chimeric BAM file. This program " << endl
	     << "extracts reads supporting read-through fusions from the RNA BAM file. " << endl
	     << "The output file should be passed to Arriba via the parameter -r." << endl
	     << "For optimal performance extract_read-through_fusions should be run " << endl
	     << "while STAR is running (see usage)." << endl << endl
	     << "Usage: # run on existing BAM file:" << endl
	     << "       extract_read-through_fusions -g annotation.gtf -i rna.bam -o read_through.bam" << endl
	     << "Usage: # run during alignment and produce sorted BAM file:" << endl
	     << "       STAR --outStd BAM_Unsorted --outSAMtype BAM Unsorted SortedByCoordinate [...] |" << endl
	     << "       extract_read-through_fusions -g annotation.gtf > read_through.bam" << endl << endl
	     << wrap_help("-i FILE", "Input file in BAM format containing alignments from STAR. "
	                  "The file need not be sorted. Default: " + default_options.input_bam_file)
	     << wrap_help("-o FILE", "Output file in BAM format containing reads which support "
	                  "read-through fusions. Default: " + default_options.output_bam_file)
	     << wrap_help("-g FILE", "GTF file with gene annotation. The file may be gzip compressed.")
	     << wrap_help("-G GTF_FEATURES", "Comma-/space-separated list of names of GTF features.\n"
	                  "Default: " + default_options.gtf_features)
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
	while ((c = getopt(argc, argv, "i:o:g:G:1h")) != -1) {
		switch (c) {
			case 'i':
				options.input_bam_file = optarg;
				if (access(options.input_bam_file.c_str(), R_OK) != 0) {
					cerr << "ERROR: File '" << options.input_bam_file << "' not found.";
					exit(1);
				}
				break;
			case 'o':
				options.output_bam_file = optarg;
				if (!output_directory_exists(options.output_bam_file)) {
					cerr << "ERROR: Parent directory of output file '" << options.output_bam_file << "' does not exist." << endl;
					exit(1);
				}
				break;
			case 'g':
				options.gene_annotation_file = optarg;
				if (access(options.gene_annotation_file.c_str(), R_OK) != 0) {
					cerr << "ERROR: File '" << options.gene_annotation_file << "' not found.";
					exit(1);
				}
				break;
			case 'G':
				options.gtf_features = optarg;
				{
					gtf_features_t gtf_features;
					if (!parse_gtf_features(options.gtf_features, gtf_features)) {
						cerr << "ERROR: Malformed GTF features: " << options.gtf_features << endl;
						exit(1);
					}
				}
				break;
			case 'h':
				print_usage();
				break;
			default:
				switch (optopt) {
					case 'g': case 'i': case 'o':
						cerr << "ERROR: " << string("Option -") + ((char) optopt) + " requires an argument." << endl;
						exit(1);
						break;
					default:
						cerr << "ERROR: " << string("Unknown option: -") + ((char) optopt) << endl;
						exit(1);
						break;
				}
		}
	}

	// check for mandatory arguments
	if (argc == 1)
		print_usage("No arguments given.");
	if (options.input_bam_file.empty()) {
		cerr << "ERROR: Missing mandatory option: -i" << endl;
		exit(1);
	}
	if (options.output_bam_file.empty()) {
		cerr << "ERROR: Missing mandatory option: -o" << endl;
		exit(1);
	}
	if (options.gene_annotation_file.empty()) {
		cerr << "ERROR: Missing mandatory option: -g" << endl;
		exit(1);
	}

	return options;
}

