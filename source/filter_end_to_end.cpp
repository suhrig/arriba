#include "common.hpp"
#include "filter_end_to_end.hpp"

using namespace std;

unsigned int filter_end_to_end_fusions(fusions_t& fusions) {

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != NULL)
			continue; // fusion has already been filtered

		if (!fusion->second.is_read_through() &&
		    (fusion->second.spliced1 || fusion->second.spliced2)) { // spliced breakpoints are likely true
			++remaining;
			continue;
		}

		if (fusion->second.discordant_mates + fusion->second.split_reads1 == 0 || // only filter breakpoints with low support
		    fusion->second.discordant_mates + fusion->second.split_reads2 == 0 ||
		    fusion->second.split_reads1 + fusion->second.split_reads2 == 0 ||
		    fusion->second.breakpoint_overlaps_both_genes() && (fusion->second.split_reads1 == 0 || fusion->second.split_reads2 == 0)) {
			if ((fusion->second.gene1->is_dummy || (fusion->second.gene1->strand == FORWARD && fusion->second.direction1 == UPSTREAM) || (fusion->second.gene1->strand == REVERSE && fusion->second.direction1 == DOWNSTREAM)) &&
			    (fusion->second.gene2->is_dummy || (fusion->second.gene2->strand == FORWARD && fusion->second.direction2 == UPSTREAM) || (fusion->second.gene2->strand == REVERSE && fusion->second.direction2 == DOWNSTREAM))) {
				fusion->second.filter = FILTERS.at("end_to_end");
				continue;
			}
		}

		++remaining;
	}

	return remaining;
}

