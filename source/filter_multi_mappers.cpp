#include "common.hpp"
#include "annotation.hpp"
#include "filter_multi_mappers.hpp"

using namespace std;

unsigned int filter_multi_mappers(chimeric_alignments_t& chimeric_alignments) {
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end();) {

		if (chimeric_alignment->second.single_end) {
			if (chimeric_alignment->second.size() == 2) {

				// use the alignment with the shorter anchor as the SUPPLEMENTARY and the longer one as the SPLIT_READ
				// and copy the split read in the MATE1 place (to simulate paired-end data)
				if (chimeric_alignment->second[MATE1].end - chimeric_alignment->second[MATE1].start > chimeric_alignment->second[MATE2].end - chimeric_alignment->second[MATE2].start) {
					chimeric_alignment->second.push_back(chimeric_alignment->second[MATE2]);
					chimeric_alignment->second[MATE2] = chimeric_alignment->second[MATE1];
				} else {
					chimeric_alignment->second.push_back(chimeric_alignment->second[MATE1]);
					chimeric_alignment->second[MATE1] = chimeric_alignment->second[MATE2];
				}

				// MATE1 and SPLIT_READ must have the sequence, SUPPLEMENTARY must not
				if (!chimeric_alignment->second[MATE1].supplementary) {
					chimeric_alignment->second[SPLIT_READ].sequence = chimeric_alignment->second[MATE1].sequence;
				} else if (!chimeric_alignment->second[SPLIT_READ].supplementary) {
					chimeric_alignment->second[MATE1].sequence = chimeric_alignment->second[SPLIT_READ].sequence;
				} else { // !chimeric_alignment->second[SUPPLEMENTARY].supplementary
					chimeric_alignment->second[MATE1].sequence = chimeric_alignment->second[SUPPLEMENTARY].sequence;
					chimeric_alignment->second[SPLIT_READ].sequence = chimeric_alignment->second[SUPPLEMENTARY].sequence;
				}
				chimeric_alignment->second[SUPPLEMENTARY].sequence.clear();

				// set supplementary flag like it would be set if we had paired-end data
				chimeric_alignment->second[SUPPLEMENTARY].supplementary = true;
				chimeric_alignment->second[MATE1].supplementary = false;
				chimeric_alignment->second[SPLIT_READ].supplementary = false;

				// set strands like they would be set if we had paired-end data
				if (chimeric_alignment->second[SPLIT_READ].sequence.length() - chimeric_alignment->second[SPLIT_READ].preclipping() - ((chimeric_alignment->second[SPLIT_READ].strand == chimeric_alignment->second[SUPPLEMENTARY].strand) ? chimeric_alignment->second[SUPPLEMENTARY].postclipping() : chimeric_alignment->second[SUPPLEMENTARY].preclipping()) <
				    chimeric_alignment->second[SPLIT_READ].sequence.length() - chimeric_alignment->second[SPLIT_READ].postclipping() - ((chimeric_alignment->second[SPLIT_READ].strand == chimeric_alignment->second[SUPPLEMENTARY].strand) ? chimeric_alignment->second[SUPPLEMENTARY].preclipping() : chimeric_alignment->second[SUPPLEMENTARY].postclipping())) {
					if (chimeric_alignment->second[SPLIT_READ].strand == FORWARD) {
						chimeric_alignment->second[MATE1].strand = complement_strand(chimeric_alignment->second[MATE1].strand);
					} else {
						chimeric_alignment->second[SPLIT_READ].strand = complement_strand(chimeric_alignment->second[SPLIT_READ].strand);
						chimeric_alignment->second[SUPPLEMENTARY].strand = complement_strand(chimeric_alignment->second[SUPPLEMENTARY].strand);
					}
				} else {
					if (chimeric_alignment->second[SPLIT_READ].strand == REVERSE) {
						chimeric_alignment->second[MATE1].strand = complement_strand(chimeric_alignment->second[MATE1].strand);
					} else {
						chimeric_alignment->second[SPLIT_READ].strand = complement_strand(chimeric_alignment->second[SPLIT_READ].strand);
						chimeric_alignment->second[SUPPLEMENTARY].strand = complement_strand(chimeric_alignment->second[SUPPLEMENTARY].strand);
					}
				}

				++chimeric_alignment;
			} else {
				chimeric_alignment = chimeric_alignments.erase(chimeric_alignment);
			}

		} else { // paired_end

			if (chimeric_alignment->second.size() == 3) { // split read

				// make sure supplementary alignment is in the SUPPLEMENTARY place
				if (chimeric_alignment->second[MATE1].supplementary) {
					swap(chimeric_alignment->second[MATE1], chimeric_alignment->second[SUPPLEMENTARY]);
				} else if (chimeric_alignment->second[MATE2].supplementary) {
					swap(chimeric_alignment->second[MATE2], chimeric_alignment->second[SUPPLEMENTARY]);
				} // else supplementary alignment is already in the correct place

				// make sure the split read is in the SPLIT_READ place
				if (chimeric_alignment->second[SPLIT_READ].first_in_pair != chimeric_alignment->second[SUPPLEMENTARY].first_in_pair)
					swap(chimeric_alignment->second[MATE1], chimeric_alignment->second[MATE2]);

				++chimeric_alignment;
			} else if (chimeric_alignment->second.size() == 2) { // discordant mate
				++chimeric_alignment;
			} else {
				// discard mates which are multi-mappers or where one mate is missing
				chimeric_alignment = chimeric_alignments.erase(chimeric_alignment);
			}

		}
	}

	return chimeric_alignments.size();
}

