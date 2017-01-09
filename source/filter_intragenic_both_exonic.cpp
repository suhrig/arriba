#include "common.hpp"
#include "filter_intragenic_both_exonic.hpp"

using namespace std;

unsigned int filter_intragenic_both_exonic(fusions_t& fusions) {
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter != NULL)
			continue; // read has already been filtered

		if ((fusion->second.breakpoint_overlaps_both_genes() || fusion->second.gene1 == fusion->second.gene2) &&
		    fusion->second.exonic1 && fusion->second.exonic2 &&
		    !fusion->second.spliced1 && !fusion->second.spliced2) {
			fusion->second.filter = FILTERS.at("intragenic_exonic");
		} else {
			++remaining;
		}
	}

	return remaining;
}

