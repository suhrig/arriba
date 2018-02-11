#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "options.hpp"
#include "options_extract_reads.hpp"
#include "extract_reads_chimeric.hpp"
#include "extract_reads_fastq.hpp"
#include "extract_reads_read_through_fusion.hpp"

using namespace std;

typedef map<string,bam1_t*> buffered_bam_records_t;

int main(int argc, char **argv) {

	options_t options = parse_arguments(argc, argv);

	// open input BAM file
	BGZF* rna_bam_file = bam_open(options.rna_bam_file.c_str(), "rb");
	bam_header_t* bam_header = bam_header_read(rna_bam_file);

	// if chimeric extraction is enabled, open chimeric output file
	BGZF* chimeric_file = NULL;
	if (!options.chimeric_file.empty()) {
		chimeric_file = bam_open(options.chimeric_file.c_str(), "w");
		if (chimeric_file == 0) {
			cerr << "ERROR: failed to open output file '" << options.chimeric_file << "'." << endl;
			exit(1);
		}
		bam_header_write(chimeric_file, bam_header);
	}

	// if read-through fusion extraction is enabled, load annotation and open output file
	gene_annotation_t gene_annotation;
	gene_annotation_index_t gene_annotation_index;
	BGZF* read_through_file = NULL;
	if (!options.read_through_file.empty()) {

		// make a mapping of contig name -> contig ID
		contigs_t contigs;
		for (int32_t i = 0; i < bam_header->n_targets; ++i)
			contigs[removeChr(bam_header->target_name[i])] = i;

		// load and index annotation
		exon_annotation_t exon_annotation;
		unordered_map<string,gene_t> gene_names;
		read_annotation_gtf(options.gene_annotation_file, contigs, options.gtf_features, gene_annotation, exon_annotation, gene_names);
		make_annotation_index(gene_annotation, gene_annotation_index, contigs);

		// open output BAM file
		read_through_file = bam_open(options.read_through_file.c_str(), "w");
		if (read_through_file == 0) {
			cerr << "ERROR: failed to open output file '" << options.read_through_file << "'." << endl;
			exit(1);
		}
		bam_header_write(read_through_file, bam_header);
	}

	// if FastQ read extraction is enabled, open FastQ files
	ofstream fastq_file1, fastq_file2;
	if (!options.fastq_file_prefix.empty()) {
		string fastq_file1_path = string(options.fastq_file_prefix) + "_1.fastq";
		fastq_file1.open(fastq_file1_path.c_str());
		if (!fastq_file1.good()) {
			cerr << "ERROR: failed to open output file '" << fastq_file1_path << "'." << endl;
			exit(1);
		}
		string fastq_file2_path = string(options.fastq_file_prefix) + "_2.fastq";
		fastq_file2.open(fastq_file2_path.c_str());
		if (!fastq_file2.good()) {
			cerr << "ERROR: failed to open output file '" << fastq_file2_path << "'." << endl;
			exit(1);
		}
	}

	// read BAM records
	bam1_t* bam_record = bam_init1();
	if (bam_record == NULL) {
		cerr << "ERROR: failed to allocate memory." << endl;
		exit(1);
	}
	buffered_bam_records_t buffered_bam_records; // holds the first mate until we have found the second
	while (bam_read1(rna_bam_file, bam_record) > 0) {

		if (bam_record->core.flag & BAM_FSECONDARY) // ignore multi-mapping reads
			continue;

		// extract chimeric alignments, if requested
		if (bam_record->core.flag & BAM_FSUPPLEMENTARY) { // supplementary alignment of a split read
			if (!options.chimeric_file.empty()) {
				bam_record->core.flag ^= BAM_FSUPPLEMENTARY; // change supplementary flag to secondary flag
				bam_record->core.flag |= BAM_FSECONDARY;
				bam_write1(chimeric_file, bam_record);
			}
			continue; // supplementary alignments are written directly; only discordant mates and split reads need to be buffered (see below)
		}

		// for paired-end data we need to wait until we have received both mates
		bam1_t* previously_seen_mate = NULL;
		if (bam_record->core.flag & BAM_FPAIRED) {

			// try to insert the mate into the buffered BAM records
			// if there was already a record with the same read name, insertion will fail (->second set to false) and
			// previously_seen_mate->first will point to the mate which was already in the buffered BAM records
			pair<buffered_bam_records_t::iterator,bool> find_previously_seen_mate = buffered_bam_records.insert(pair<string,bam1_t*>((char*) bam1_qname(bam_record), bam_record));
			if (!find_previously_seen_mate.second) { // this is the second mate we have seen
				previously_seen_mate = find_previously_seen_mate.first->second;
				buffered_bam_records.erase(find_previously_seen_mate.first); // remove from lookup buffer, we don't need it anymore
			}

		}

		if ((bam_record->core.flag & BAM_FPAIRED) && previously_seen_mate == NULL) { // this is the first mate with the given read name, which we encounter
			
			bam_record = bam_init1(); // allocate memory for the next record
			if (bam_record == NULL) {
				cerr << "ERROR: failed to allocate memory." << endl;
				exit(1);
			}

		} else { // single-end data or we have already read the first mate previously

			if (!options.chimeric_file.empty())
				extract_chimeric(chimeric_file, bam_record, previously_seen_mate);

			if (!options.fastq_file_prefix.empty())
				extract_fastq(fastq_file1, fastq_file2, bam_record, previously_seen_mate, options.min_clipped_length);

			if (!options.read_through_file.empty())
				extract_read_through_fusion(read_through_file, bam_record, previously_seen_mate, gene_annotation_index);

			if (previously_seen_mate != NULL)
				bam_destroy1(previously_seen_mate);
		}

	}

	// close output files
	bam_destroy1(bam_record);
	bam_header_destroy(bam_header);
	bam_close(rna_bam_file);
	if (!options.chimeric_file.empty())
		bam_close(chimeric_file);
	if (!options.read_through_file.empty())
		bam_close(read_through_file);
	if (!options.fastq_file_prefix.empty()) {
		fastq_file1.close();
		fastq_file2.close();
	}

	return 0;
}

