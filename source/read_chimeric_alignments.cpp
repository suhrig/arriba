#include <algorithm>
#include <iostream>
#include <string>
#include "sam.h"
#include "annotation.hpp"
#include "common.hpp"
#include "read_chimeric_alignments.hpp"

using namespace std;

unsigned int read_chimeric_alignments(const string& bam_file_path, chimeric_alignments_t& chimeric_alignments, contigs_t& contigs, const bool read_through) {

	chimeric_alignments_t read_through_alignments; // we only need this, when reading read-through alignments

	// open BAM file
	BGZF* bam_file = bam_open(bam_file_path.c_str(), "rb");
	bam_header_t* bam_header = bam_header_read(bam_file);

	if (read_through) {
		// check if chimeric.bam and read_through.bam have the same contigs
		for (int32_t target = 0; target < bam_header->n_targets; ++target)
			if (contigs.find(removeChr(bam_header->target_name[target])) == contigs.end()) {
				cerr << "ERROR: Input BAM files have different contigs." << endl;
				exit(1);
			}
	} else {
		// make a mapping of contig name -> contig ID
		for (int32_t target = 0; target < bam_header->n_targets; ++target)
			contigs[removeChr(bam_header->target_name[target])] = target;
	}

	// read BAM records
	bam1_t* bam_record = bam_init1();
	while (bam_read1(bam_file, bam_record) > 0) {

		alignment_t alignment;
		string name = (char*) bam1_qname(bam_record);
		alignment.supplementary = bam_record->core.flag & BAM_FSECONDARY;
		alignment.strand = (bam_record->core.flag & BAM_FREVERSE) ? REVERSE : FORWARD;
		alignment.first_in_pair = bam_record->core.flag & BAM_FREAD1;
		alignment.contig = bam_record->core.tid;
		alignment.start = bam_record->core.pos;
		alignment.end = bam_endpos(bam_record);
		alignment.cigar = bam_record;
//		if (!alignment.supplementary) {
			char sequence[bam_record->core.l_qseq+1];
			sequence[bam_record->core.l_qseq] = '\0';
			for (unsigned int i = 0; i < bam_record->core.l_qseq; ++i)
				sequence[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(bam_record), i)];
			alignment.sequence = sequence;
			// convert to uppercase
			std::transform(alignment.sequence.begin(), alignment.sequence.end(), alignment.sequence.begin(), (int (*)(int))std::toupper);
//		}
		
		if (!read_through) {
			mates_t& mates = chimeric_alignments[name];
			mates.push_back(alignment);
			mates.name = name;
		} else {
			mates_t& mates = read_through_alignments[name];
			mates.push_back(alignment);
			mates.name = name;
		}
	}

	// close BAM file
	bam_destroy1(bam_record);
	bam_header_destroy(bam_header);
	bam_close(bam_file);

	// append read_through_alignments to chimeric_alignments if the alignment is not already contained in the latter
	if (read_through)
		for (chimeric_alignments_t::iterator read_through_alignment = read_through_alignments.begin(); read_through_alignment != read_through_alignments.end();) {
			chimeric_alignments.insert(make_pair(read_through_alignment->first, read_through_alignment->second)); // only succeeds, if does not exist in chimeric_alignments
			read_through_alignment = read_through_alignments.erase(read_through_alignment);
		}

	return chimeric_alignments.size();
}

void assign_strands_from_strandedness(chimeric_alignments_t& chimeric_alignments, const strandedness_t strandedness) {
	if (strandedness != STRANDEDNESS_NO) {
		for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {
			unsigned int first_in_pair = (chimeric_alignment->second[MATE1].first_in_pair) ? MATE1 : MATE2;
			unsigned int second_in_pair = (chimeric_alignment->second[MATE1].first_in_pair) ? MATE2 : MATE1;
			chimeric_alignment->second[first_in_pair].predicted_strand = complement_strand_if(chimeric_alignment->second[first_in_pair].strand, strandedness == STRANDEDNESS_REVERSE);
			chimeric_alignment->second[first_in_pair].predicted_strand_ambiguous = false;
			chimeric_alignment->second[second_in_pair].predicted_strand = complement_strand_if(chimeric_alignment->second[first_in_pair].predicted_strand, chimeric_alignment->second[first_in_pair].strand == chimeric_alignment->second[second_in_pair].strand);
			chimeric_alignment->second[second_in_pair].predicted_strand_ambiguous = false;
			if (chimeric_alignment->second.size() == 3) { // split read
				chimeric_alignment->second[SUPPLEMENTARY].predicted_strand = complement_strand_if(chimeric_alignment->second[SPLIT_READ].predicted_strand, chimeric_alignment->second[SUPPLEMENTARY].strand != chimeric_alignment->second[SPLIT_READ].strand);
				chimeric_alignment->second[SUPPLEMENTARY].predicted_strand_ambiguous = false;
			}
		}
	}
}

