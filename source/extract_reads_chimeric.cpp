#include <iostream>
#include "sam.h"
#include "extract_reads_chimeric.hpp"

using namespace std;

void write_chimeric_alignment(samFile* bam_file, bam_hdr_t* bam_header, const bam1_t* const bam_record) {
	if (sam_write1(bam_file, bam_header, bam_record) < 0) {
		cerr << "ERROR: failed to write BAM record to chimeric alignments file." << endl;
		exit(1);
	}
}

bool is_chimeric_alignment(const bam1_t* const bam_record) {
	// check if a read is a discordant mate or
	// a split read, by checking if it has a supplementary alignment (SA tag must be present)
	if (bam_record == NULL)
		return false;
	else
		return ((bam_record->core.flag & BAM_FPAIRED) && !(bam_record->core.flag & BAM_FPROPER_PAIR) && !(bam_record->core.flag & BAM_FMUNMAP) || // discordant mate
		        bam_aux_get(bam_record, "SA") != NULL); // supplementary alignment
}

bool extract_supplementary(samFile* chimeric_file, bam_hdr_t* bam_header, bam1_t* bam_record) {
	if (bam_record->core.flag & BAM_FSUPPLEMENTARY) { // supplementary alignment of a split read
		bam_record->core.flag ^= BAM_FSUPPLEMENTARY; // change supplementary flag to secondary flag
		bam_record->core.flag |= BAM_FSECONDARY;
		write_chimeric_alignment(chimeric_file, bam_header, bam_record);
		return true; // this is a supplementary alignment
	}
	return false; // this is not a supplementary alignment
}

void extract_chimeric(samFile* chimeric_file, bam_hdr_t* bam_header, const bam1_t* const read1, const bam1_t* const read2) {
	if (is_chimeric_alignment(read1) || is_chimeric_alignment(read2)) {
		if (read1 != NULL)
			write_chimeric_alignment(chimeric_file, bam_header, read1);
		if (read2 != NULL)
			write_chimeric_alignment(chimeric_file, bam_header, read2);
	}
}

