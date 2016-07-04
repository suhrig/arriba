#include "common.hpp"
#include "filter_inconsistently_clipped.hpp"

using namespace std;

unsigned int filter_inconsistently_clipped_mates(chimeric_alignments_t& chimeric_alignments) {

	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {
		if (!i->second.filters.empty())
			continue; // read has already been filtered

		if (i->second.size() == 3) { // these are alignments of a split read
			if ((i->second[MATE1].strand == FORWARD && i->second[MATE1].end > i->second[SPLIT_READ].end+3) ||
			    (i->second[MATE1].strand == REVERSE && i->second[MATE1].start < i->second[SPLIT_READ].start-3)) {
				i->second.filters.insert(FILTERS.at("inconsistently_clipped"));
				continue;
			}
		}

		++remaining;
	}

	return remaining;
}

