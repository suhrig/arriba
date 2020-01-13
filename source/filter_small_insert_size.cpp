#include <cmath>
#include "common.hpp"
#include "filter_small_insert_size.hpp"

using namespace std;

unsigned int filter_small_insert_size(chimeric_alignments_t& chimeric_alignments, const unsigned int max_overhang) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != FILTER_none)
			continue; // the read has already been filtered

		// remove chimeric alignment when insert size is too small
		if (chimeric_alignment->second.size() == 2) { // discordant mates
			if (chimeric_alignment->second[MATE1].strand != chimeric_alignment->second[MATE2].strand &&
			    chimeric_alignment->second[MATE1].contig == chimeric_alignment->second[MATE2].contig &&
			    (abs(chimeric_alignment->second[MATE1].start - chimeric_alignment->second[MATE2].start) <= max_overhang ||
			     abs(chimeric_alignment->second[MATE1].end - chimeric_alignment->second[MATE2].end) <= max_overhang)) {
				chimeric_alignment->second.filter = FILTER_small_insert_size;
				continue;
			}
		}

		// we only get here, if the chimeric alignments were not filtered
		++remaining;
	}

	return remaining;
}

