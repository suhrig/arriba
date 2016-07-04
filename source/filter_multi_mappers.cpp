#include "common.hpp"
#include "filter_multi_mappers.hpp"

using namespace std;

unsigned int filter_multi_mappers(chimeric_alignments_t& chimeric_alignments) {

	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end();) {
		if (i->second.size() == 3) { // split read

			// make sure supplementary alignment is in the SUPPLEMENTARY_ALIGNMENT place
			if (i->second[MATE1].supplementary) {
				swap(i->second[MATE1], i->second[SUPPLEMENTARY]);
			} else if (i->second[MATE2].supplementary) {
				swap(i->second[MATE2], i->second[SUPPLEMENTARY]);
			} // else supplementary alignment is already in the correct place

			// make sure the split read is in the SPLIT_READ place
			if (i->second[SPLIT_READ].first_in_pair != i->second[SUPPLEMENTARY].first_in_pair) {
				swap(i->second[MATE1], i->second[MATE2]);
			}

			++i;
		} else if (i->second.size() == 2) { // discordant mate
			++i;
		} else {
			// discard mates which are multi-mappers or where one mate is missing
			i = chimeric_alignments.erase(i);
		}
	}

	return chimeric_alignments.size();
}

