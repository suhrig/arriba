#include "common.hpp"
#include "filter_non_coding_neighbors.hpp"

using namespace std;

unsigned int filter_non_coding_neighbors(fusions_t& fusions) {
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter != FILTER_none)
			continue; // read has already been filtered

		if (!fusion->second.gene1->is_protein_coding && !fusion->second.gene2->is_protein_coding &&
		    fusion->second.is_read_through())
			fusion->second.filter = FILTER_non_coding_neighbors;
		else
			++remaining;
	}

	return remaining;
}

