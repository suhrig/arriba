#include "common.hpp"
#include "filter_both_intronic.hpp"

using namespace std;

bool list_contains_exonic_reads(const vector<mates_t*>& read_list) {
	for (auto i = read_list.begin(); i != read_list.end(); ++i)
		if ((**i).filter == NULL)
			for (mates_t::iterator mate = (**i).begin(); mate != (**i).end(); ++mate)
				if (mate->exonic)
					return true;
	return false;
}

unsigned int filter_both_intronic(fusions_t& fusions) {
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter != NULL)
			continue; // read has already been filtered

		if (!list_contains_exonic_reads(fusion->second.split_read1_list) &&
		    !list_contains_exonic_reads(fusion->second.split_read2_list) &&
		    !list_contains_exonic_reads(fusion->second.discordant_mate_list) &&
		    fusion->second.closest_genomic_breakpoint1 < 0) {
			fusion->second.filter = FILTERS.at("intronic");
		} else {
			++remaining;
		}
	}

	return remaining;
}

