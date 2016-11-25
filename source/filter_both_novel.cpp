#include "common.hpp"
#include "filter_both_novel.hpp"

using namespace std;

unsigned int filter_both_novel(fusions_t& fusions) {
	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
		if (!i->second.filters.empty())
			continue; // read has already been filtered

		if (!i->second.gene1->is_known && !i->second.gene2->is_known)
			i->second.filters.insert(FILTERS.at("novel"));
		else
			++remaining;
	}

	return remaining;
}

