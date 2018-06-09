#include <iostream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <unistd.h>
#include <algorithm>
#include "annotation.hpp"
#include "options.hpp"
#include "options_extract_reads.hpp"

using namespace std;

options_t get_default_options() {
	options_t options;

	options.rna_bam_file = "/dev/stdin";
	options.gtf_features = DEFAULT_GTF_FEATURES;
	options.min_clipped_length = 10;

	return options;
}

void print_usage() {

	options_t default_options = get_default_options();

	cout << endl
	     << "Arriba gene fusion detector - extract_reads" << endl
	     << "-------------------------------------------" << endl
	     << "Version: " << ARRIBA_VERSION << endl
	     << endl
	     << "This is a helper utility of Arriba. It has three modes, each of " << endl
	     << "which extracts a certain type of reads from the normal alignments " << endl
	     << "of STAR (Aligned.out.bam)." << endl
	     << endl
	     << "read-through fusions mode (-r):" << endl
	     << "  The STAR RNA-Seq aligner does not report read-through fusions " << endl
	     << "  in the chimeric alignments file (Chimeric.out.sam). This program " << endl
	     << "  extracts reads supporting read-through fusions from the normal " << endl
	     << "  alignments file. The output file should be passed to Arriba via " << endl
	     << "  the parameter -r." << endl
	     << "  For optimal performance extract_reads should be run while " << endl
	     << "  STAR is running (see usage)." << endl
	     << "FastQ mode (-f):" << endl
	     << "  If STAR is run without the parameters --chimSegmentMin or " << endl
	     << "  --chimOutType SeparateSAMold, then the output file " << endl
	     << "  Chimeric.out.sam is not generated. Arriba requires this file " << endl
	     << "  for fusion detection, however. In order to avoid having to " << endl
	     << "  rerun STAR on the already aligned reads this utility can be " << endl
	     << "  used to extract only those reads from the normal alignments " << endl
	     << "  file, which are eligible for chimeric alignment, i.e., unmapped " << endl
	     << "  reads, discordant mates, and clipped reads. The reads are " << endl
	     << "  extracted in FastQ format to align them anew with STAR." << endl
	     << "chimeric alignments mode (-c):" << endl
	     << "  If STAR was run with the parameter --chimOutType WithinBAM, " << endl
	     << "  then the chimeric alignments need to be extracted to a separate " << endl
	     << "  file. This file can then be passed to Arriba via the parameter -c." << endl
	     << endl
	     << "It is possible to enable multiple modes at the same time to extract " << endl
	     << "various types of reads simultaneously." << endl
	     << endl
	     << "Usage: # extract read-through fusions from existing BAM file:" << endl
	     << "       extract_reads -g annotation.gtf -x Aligned.out.bam -r read_through.bam" << endl
	     << "Usage: # extract read-through fusions during alignment and produce sorted BAM file:" << endl
	     << "       STAR --outStd BAM_Unsorted --outSAMtype BAM Unsorted SortedByCoordinate [...] |" << endl
	     << "       extract_reads -g annotation.gtf -r read_through.bam" << endl
	     << "Usage: # extract chimeric candidates in FastQ format from existing BAM file:" << endl
	     << "       extract_reads -x Aligned.out.bam -f reads_as_fastq" << endl
	     << "Usage: # extract chimeric alignments from existing BAM file:" << endl
	     << "       extract_reads -x Aligned.out.bam -c chimeric.bam" << endl
	     << endl
	     << wrap_help("-x FILE", "Input file in SAM/BAM/CRAM format containing alignments from STAR. "
	                  "The file does not need to be sorted. Default: " + default_options.rna_bam_file)
	     << wrap_help("-r FILE", "Output file containing reads which support read-through fusions. "
	                  "The file has the same format as the input file (-x).")
	     << wrap_help("-f FILE", "Prefix of output files in FastQ format containing reads which "
	                  "are eligible for chimeric alignment, i.e., unmapped reads, discordant "
	                  "mates, and clipped reads.")
	     << wrap_help("-c FILE", "Output file containing chimeric reads, i.e., discordant mates "
	                  "and split reads with supplementary alignments. The file has the same format "
	                  "as the input file (-x).")
	     << wrap_help("-a FILE", "Assembly in FastA format. This option is only required "
	                  "when CRAM files are processed. An index (.fai) must exist.")
	     << wrap_help("-g FILE", "GTF file with gene annotation. The file may be gzip-compressed.")
	     << wrap_help("-G GTF_FEATURES", "Comma-/space-separated list of names of GTF features.\n"
	                  "Default: " + default_options.gtf_features)
	     << wrap_help("-l MIN_CLIPPED_LENGTH", "When -f is specified, only those reads are "
	                  "extracted which have a clipped segment of the given length or longer. "
	                  "This argument should match the argument to --chimSegmentMin of STAR. "
	                  "Default: " + to_string(static_cast<long long unsigned int>(default_options.min_clipped_length)))
	     << wrap_help("-h", "Print help and exit.")
	     << "For more information or help, visit: " << HELP_CONTACT << endl;
}

options_t parse_arguments(int argc, char **argv) {

	options_t options = get_default_options();

	// parse arguments
	opterr = 0;
	int c;
	while ((c = getopt(argc, argv, "x:r:f:c:a:g:G:l:h")) != -1) {
		switch (c) {
			case 'x':
				options.rna_bam_file = optarg;
				if (access(options.rna_bam_file.c_str(), R_OK) != 0) {
					cerr << "ERROR: File '" << options.rna_bam_file << "' not found.";
					exit(1);
				}
				break;
			case 'r':
				options.read_through_file = optarg;
				if (!output_directory_exists(options.read_through_file)) {
					cerr << "ERROR: Parent directory of output file '" << options.read_through_file << "' does not exist." << endl;
					exit(1);
				}
				break;
			case 'f':
				options.fastq_file_prefix = optarg;
				if (!output_directory_exists(options.fastq_file_prefix)) {
					cerr << "ERROR: Parent directory of output file '" << options.fastq_file_prefix << "' does not exist." << endl;
					exit(1);
				}
				break;
			case 'c':
				options.chimeric_file = optarg;
				if (!output_directory_exists(options.chimeric_file)) {
					cerr << "ERROR: Parent directory of output file '" << options.chimeric_file << "' does not exist." << endl;
					exit(1);
				}
				break;
			case 'a':
				options.assembly_file = optarg;
				if (access(options.assembly_file.c_str(), R_OK) != 0) {
					cerr << "ERROR: File '" << options.assembly_file << "' not found." << endl;
					exit(1);
				}
				if (access((options.assembly_file + ".fai").c_str(), R_OK) != 0) {
					cerr << "ERROR: Index for '" << options.assembly_file << "' not found." << endl;
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
			case 'l':
				if (!validate_int(optarg, options.min_clipped_length, 0)) {
					cerr << "ERROR: " << "Invalid argument to -" << ((char) c) << "." << endl;
					exit(1);
				}
				break;
			case 'h':
				print_usage();
				exit(0);
				break;
			default:
				switch (optopt) {
					case 'x': case 'r': case 'f': case 'c': case 'a': case 'g': case 'G': case 'l':
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
	if (argc == 1) {
		cerr << "ERROR: No arguments given." << endl;
		print_usage();
		exit(1);
	}
	if (options.rna_bam_file.empty()) {
		cerr << "ERROR: Missing mandatory option: -x" << endl;
		exit(1);
	}
	if (options.read_through_file.empty() && options.fastq_file_prefix.empty() && options.chimeric_file.empty()) {
		cerr << "ERROR: One of the options -r, -f, or -c must be specified." << endl;
		exit(1);
	}
	if (!options.read_through_file.empty() && options.gene_annotation_file.empty()) {
		cerr << "ERROR: Option -g must be specified, when -r is given." << endl;
		exit(1);
	}
	if (options.read_through_file.empty() && !options.gene_annotation_file.empty())
		cerr << "WARNING: Option -g has no effect, when -r is not given." << endl;

	return options;
}

