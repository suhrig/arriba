#include <cmath>
#include "common.hpp"
#include "filter_short_anchor.hpp"

using namespace std;

unsigned int filter_short_anchor(fusions_t& fusions, unsigned int min_length) {
	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (!i->second.filters.empty())
			continue; // fusion has already been filtered

		if (abs(i->second.anchor_start1 - i->second.breakpoint1) < min_length ||
		    abs(i->second.anchor_start2 - i->second.breakpoint2) < min_length) {
			i->second.filters.insert(FILTERS.at("short_anchor"));
		} else {
			remaining++;
		}
	}
	return remaining;
}

