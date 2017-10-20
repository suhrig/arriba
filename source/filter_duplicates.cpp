#include <unordered_map>
#include <tuple>
#include "common.hpp"
#include "filter_duplicates.hpp"

using namespace std;

unsigned int filter_duplicates(chimeric_alignments_t& chimeric_alignments) {
	unsigned int remaining = 0;
	unordered_map< tuple<position_t,position_t> , unsigned int> duplicate_count;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {
		if (chimeric_alignment->second.filter != NULL)
			continue; // read has already been filtered

		unsigned int mate2 = (chimeric_alignment->second.size() == 2) ? MATE2 : SUPPLEMENTARY;

		tuple<position_t,position_t> mate_coordinates = make_tuple(
			(chimeric_alignment->second[MATE1].strand == FORWARD) ?
				chimeric_alignment->second[MATE1].start - chimeric_alignment->second[MATE1].preclipping() :
				chimeric_alignment->second[MATE1].end   + chimeric_alignment->second[MATE1].postclipping(),
			(chimeric_alignment->second[mate2].strand == FORWARD) ?
				chimeric_alignment->second[mate2].start - chimeric_alignment->second[mate2].preclipping() :
				chimeric_alignment->second[mate2].end   + chimeric_alignment->second[mate2].postclipping()
		);

		if (duplicate_count[mate_coordinates]++ > 0)
			chimeric_alignment->second.filter = FILTERS.at("duplicates");
		else
			++remaining;
	}

	return remaining;
}

