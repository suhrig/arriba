#include "common.hpp"
#include "filter_min_support.hpp"

using namespace std;

// throw away fusions with few supporting reads
unsigned int filter_min_support(fusions_t& fusions, const int min_support) {
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != FILTER_none)
			continue; // fusion has already been filtered

		if (fusion->second.split_reads1 + fusion->second.split_reads2 + fusion->second.discordant_mates < min_support ||
		    fusion->second.breakpoint_overlaps_both_genes() && fusion->second.split_reads1 + fusion->second.split_reads2 < min_support)
			fusion->second.filter = FILTER_min_support;
		else
			remaining++;
	}
	return remaining;
}
