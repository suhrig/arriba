#include "common.hpp"
#include "filter_both_intronic.hpp"

using namespace std;

unsigned int filter_both_intronic(chimeric_alignments_t& chimeric_alignments) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {
		if (!i->second.filters.empty())
			continue; // read has already been filtered

		if ((i->second.size() == 2 && !i->second[MATE1].exonic && !i->second[MATE2].exonic) ||
		    (i->second.size() == 3 && !i->second[SPLIT_READ].exonic && !i->second[SUPPLEMENTARY].exonic)) {
			i->second.filters.insert(FILTERS.at("intronic"));
		} else {
			++remaining;
		}
	}

	return remaining;
}

