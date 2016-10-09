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
		for (int32_t i = 0; i < bam_header->n_targets; ++i)
			if (contigs.find(removeChr(bam_header->target_name[i])) == contigs.end()) {
				cerr << "ERROR: Input BAM files have different contigs." << endl;
				exit(1);
			}
	} else {
		// make a mapping of contig name -> contig ID
		for (int32_t i = 0; i < bam_header->n_targets; ++i)
			contigs[removeChr(bam_header->target_name[i])] = i;
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
		if (!alignment.supplementary) {
			char sequence[bam_record->core.l_qseq+1];
			sequence[bam_record->core.l_qseq] = '\0';
			for (unsigned int i = 0; i < bam_record->core.l_qseq; ++i)
				sequence[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(bam_record), i)];
			alignment.sequence = sequence;
			// convert to uppercase
			std::transform(alignment.sequence.begin(), alignment.sequence.end(), alignment.sequence.begin(), (int (*)(int))std::toupper);
		}
		
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
		for (chimeric_alignments_t::iterator i = read_through_alignments.begin(); i != read_through_alignments.end();) {
			chimeric_alignments.insert(make_pair(i->first, i->second)); // only succeeds, if does not exist in chimeric_alignments
			i = read_through_alignments.erase(i);
		}

	return chimeric_alignments.size();
}

