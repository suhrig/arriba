#include "common.hpp"
#include "filter_both_intronic.hpp"

using namespace std;

bool list_contains_exonic_reads(const vector<mates_t*>& read_list) {
	for (auto i = read_list.begin(); i != read_list.end(); ++i)
		if ((**i).filters.empty())
			for (mates_t::iterator mate = (**i).begin(); mate != (**i).end(); ++mate)
				if (mate->exonic)
					return true;
	return false;
}

unsigned int filter_both_intronic(fusions_t& fusions) {
	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
		if (!i->second.filters.empty())
			continue; // read has already been filtered

		if (!list_contains_exonic_reads(i->second.split_read1_list) &&
		    !list_contains_exonic_reads(i->second.split_read2_list) &&
		    !list_contains_exonic_reads(i->second.discordant_mate_list)) {
			i->second.filters.insert(FILTERS.at("intronic"));
		} else {
			++remaining;
		}
	}

	return remaining;
}

