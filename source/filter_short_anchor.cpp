#include <cmath>
#include "common.hpp"
#include "filter_short_anchor.hpp"

using namespace std;

unsigned int filter_short_anchor(fusions_t& fusions, unsigned int min_length) {
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (!fusion->second.filters.empty())
			continue; // fusion has already been filtered

		if (!(fusion->second.spliced1 && fusion->second.spliced2) &&
		    (abs(fusion->second.anchor_start1 - fusion->second.breakpoint1) < min_length ||
		     abs(fusion->second.anchor_start2 - fusion->second.breakpoint2) < min_length)) {
			fusion->second.filters.insert(FILTERS.at("short_anchor"));
		} else {
			remaining++;
		}
	}
	return remaining;
}

