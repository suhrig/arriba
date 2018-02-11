#include <iostream>
#include "sam.h"
#include "extract_reads_chimeric.hpp"

using namespace std;

void write_chimeric_alignment(BGZF* bam_file, const bam1_t* const bam_record) {
	if (bam_write1(bam_file, bam_record) < 0) {
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

void extract_chimeric(BGZF* chimeric_file, const bam1_t* const read1, const bam1_t* const read2) {
	if (is_chimeric_alignment(read1) || is_chimeric_alignment(read2)) {
		if (read1 != NULL)
			write_chimeric_alignment(chimeric_file, read1);
		if (read2 != NULL)
			write_chimeric_alignment(chimeric_file, read2);
	}
}

