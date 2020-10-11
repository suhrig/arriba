#include "sam.h"
#include "common.hpp"
#include "filter_long_gap.hpp"

using namespace std;

unsigned int filter_long_gap(chimeric_alignments_t& chimeric_alignments) {

	// If the parameter alignIntronMax of STAR is set large (>1Mbp), then occassionally
	// STAR finds an alignment with a long gap and short matching segments, which happen to match by chance, e.g.: 12M832512N13M25S
	// Particularly ostensible deletions are prone to this, where an alignment can have multiple long gaps and
	// short matching segments, e.g.: 49M902241N14M104923N12M25S
	// In the previous example, the gap of length 902241 could be a candidate for a deletion.
	// => If we see deletions of ~1Mbp and short matching segments OR alignments with long gaps and short matching segments,
	//    then we discard the alignment.

	const int min_long_gap = 700000; // we consider gaps of this size (or longer) to be too long
	const int max_long_gap = 1500000; // let's hope nobody sets alignIntronMax greater than this
	const unsigned int short_segment = 15; // we consider aligned segments of this size (or shorter) to be too short

	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != FILTER_none)
			continue; // read has already been filtered

		// check if event is a deletion between min_long_gap and max_long_gap in size
		int size_of_deletion = 0;
		if (chimeric_alignment->second.size() == 3) { // split-read
			if (chimeric_alignment->second[SPLIT_READ].contig == chimeric_alignment->second[SUPPLEMENTARY].contig) {
				if (chimeric_alignment->second[SPLIT_READ].strand == REVERSE && chimeric_alignment->second[SUPPLEMENTARY].strand == REVERSE) {
					size_of_deletion = chimeric_alignment->second[SUPPLEMENTARY].start - chimeric_alignment->second[SPLIT_READ].end;
				} else if (chimeric_alignment->second[SPLIT_READ].strand == FORWARD && chimeric_alignment->second[SUPPLEMENTARY].strand == FORWARD) {
					size_of_deletion = chimeric_alignment->second[SPLIT_READ].start - chimeric_alignment->second[SUPPLEMENTARY].end;
				}
			}
		}

		for (mates_t::iterator mate = chimeric_alignment->second.begin(); mate != chimeric_alignment->second.end(); ++mate) {

			// look for long gap
			for (unsigned int i = 1; i < mate->cigar.size()-1; ++i) {
				if (mate->cigar.operation(i) == BAM_CREF_SKIP && ((int) mate->cigar.op_length(i) >= min_long_gap || size_of_deletion >= min_long_gap && size_of_deletion <= max_long_gap)) {

					// look for short matching segment flanking the gap on the left
					unsigned int matching_segment_left = 0;
					for (int j = i-1; j >= 0; --j) {
						switch (mate->cigar.operation(j)) {
							case BAM_CMATCH: case BAM_CDIFF: case BAM_CEQUAL:
								matching_segment_left += mate->cigar.op_length(j); // sum up length of matching segment
								break;
							case BAM_CDEL: case BAM_CINS: case BAM_CPAD:
								break; // ignore indels
							default:
								goto end_of_loop_left; // end of matching segment
						}
					}
					end_of_loop_left:

					// look for short matching segment flanking the gap on the right
					unsigned int matching_segment_right = 0;
					for (unsigned int j = i+1; j < mate->cigar.size(); ++j) {
						switch (mate->cigar.operation(j)) {
							case BAM_CMATCH: case BAM_CDIFF: case BAM_CEQUAL:
								matching_segment_right += mate->cigar.op_length(j); // sum up length of matching_segment
								break;
							case BAM_CDEL: case BAM_CINS: case BAM_CPAD:
								break; // ignore indels
							default:
								goto end_of_loop_right; // end of matching segment
						}
					}
					end_of_loop_right:

					if (matching_segment_left <= short_segment && matching_segment_right <= short_segment) {
						chimeric_alignment->second.filter = FILTER_long_gap;
						goto next_read;
					}
				}
			}
		}

		remaining++; // is skipped, when the read has been filtered

		next_read: continue;
	}

	return remaining;
}

