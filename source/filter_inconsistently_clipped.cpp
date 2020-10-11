#include "common.hpp"
#include "filter_inconsistently_clipped.hpp"

using namespace std;

unsigned int filter_inconsistently_clipped_mates(chimeric_alignments_t& chimeric_alignments) {

	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {
		if (chimeric_alignment->second.filter != FILTER_none)
			continue; // read has already been filtered

		if (chimeric_alignment->second.size() == 3) { // these are alignments of a split read
			if ((chimeric_alignment->second[MATE1].strand == FORWARD && chimeric_alignment->second[MATE1].end > chimeric_alignment->second[SPLIT_READ].end+3) ||
			    (chimeric_alignment->second[MATE1].strand == REVERSE && chimeric_alignment->second[MATE1].start < chimeric_alignment->second[SPLIT_READ].start-3)) {
				chimeric_alignment->second.filter = FILTER_inconsistently_clipped;
				continue;
			}
		}

		++remaining;
	}

	return remaining;
}

