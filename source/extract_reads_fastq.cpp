#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include "sam.h"
#include "common.hpp"
#include "assembly.hpp"
#include "extract_reads_fastq.hpp"

using namespace std;

void bam_record_to_fastq(ofstream& fastq_file, const bam1_t* const bam_record) {

	// convert sequence and qualities to printable string
	string sequence, qualities;
	sequence.resize(bam_record->core.l_qseq);
	qualities.resize(bam_record->core.l_qseq);
	for (unsigned int i = 0; i < bam_record->core.l_qseq; ++i) {
		sequence[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(bam_record), i)];
		qualities[i] = 33 + bam1_qual(bam_record)[i];
	}

	// reverse complement, if the alignment is on the reverse strand
	if (bam_record->core.flag & BAM_FREVERSE) {
		sequence = dna_to_reverse_complement(sequence);
		reverse(qualities.begin(), qualities.end());
	}

	// write to FastQ file
	fastq_file << "@" << ((char*) bam1_qname(bam_record)) << endl
	           << sequence << endl
	           << "+" << endl
	           << qualities << endl;
	if (!fastq_file.good()) {
		cerr << "ERROR: failed to write to FastQ file." << endl;
		exit(1);
	}
}

bool is_read_clipped(const bam1_t* const bam_record, const uint32_t min_clipped_length) {
	bool first_cigar_op_clipped = (bam1_cigar(bam_record)[0] & 15) == BAM_CSOFT_CLIP &&
	                              (bam1_cigar(bam_record)[0] >> 4) >= min_clipped_length;
	bool last_cigar_op_clipped = (bam1_cigar(bam_record)[bam_record->core.n_cigar-1] & 15) == BAM_CSOFT_CLIP &&
	                             (bam1_cigar(bam_record)[bam_record->core.n_cigar-1] >> 4) >= min_clipped_length;
	if (bam_record->core.flag & BAM_FPAIRED) {
		// in case of paired-end data, clipping must be at end of fragment
		if (bam_record->core.flag & BAM_FREVERSE) { // read is on reverse strand
			return last_cigar_op_clipped;
		} else { // read is on forward strand
			return first_cigar_op_clipped;
		}
	} else { 
		// in case of single-end data, we accept clipping at either end
		return first_cigar_op_clipped || last_cigar_op_clipped;
	}
}

void extract_fastq(ofstream& fastq_file1, ofstream& fastq_file2, bam1_t* read1, bam1_t* read2, const uint32_t min_clipped_length) {

	// make sure the first-in-pair is read1
	if (read1->core.flag & BAM_FREAD2)
		swap(read1, read2);

	// check for clipped segment or discordant mates
	if ((read1->core.flag & BAM_FUNMAP) || // read is unmapped
	    is_read_clipped(read1, min_clipped_length) || // read is clipped
	    ((read1->core.flag & BAM_FPAIRED) && // in case of paired-end data, also check mate
	     ((read2->core.flag & BAM_FUNMAP) || // mate is unmapped
	      !(read1->core.flag & BAM_FPROPER_PAIR) || // discordant mates
	      is_read_clipped(read2, min_clipped_length)))) { // mate is clipped
		bam_record_to_fastq(fastq_file1, read1);
		if (read1->core.flag & BAM_FPAIRED)
			bam_record_to_fastq(fastq_file2, read2);
	}
}

