#include <unordered_map>
#include <tuple>
#include "common.hpp"
#include "filter_duplicates.hpp"

using namespace std;

unsigned int filter_duplicates(chimeric_alignments_t& chimeric_alignments) {
	unsigned int remaining = 0;
	unordered_map< tuple<position_t,position_t,position_t,position_t> , unsigned int> duplicate_count;
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {
		if (!i->second.filters.empty())
			continue; // read has already been filtered

		tuple<position_t,position_t,position_t,position_t> mate_coordinates = make_tuple(i->second[MATE1].start - i->second[MATE1].preclipping(),
												 i->second[MATE1].end   + i->second[MATE1].postclipping(),
												 i->second[MATE2].start - i->second[MATE2].preclipping(),
												 i->second[MATE2].end   + i->second[MATE2].postclipping());
		if (duplicate_count[mate_coordinates]++ > 0)
			i->second.filters.insert(FILTERS.at("duplicates"));
		else
			++remaining;
	}

	return remaining;
}

