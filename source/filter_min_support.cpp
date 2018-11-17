#include "common.hpp"
#include "filter_min_support.hpp"

using namespace std;

// throw away fusions with few supporting reads
unsigned int filter_min_support(fusions_t& fusions, const int min_support) {
	unsigned int remaining = 0;
	for (fusions_t::iterator chimeric_alignment = fusions.begin(); chimeric_alignment != fusions.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != NULL)
			continue; // fusion has already been filtered

		if (chimeric_alignment->second.split_reads1 + chimeric_alignment->second.split_reads2 + chimeric_alignment->second.discordant_mates < min_support ||
		    chimeric_alignment->second.breakpoint_overlaps_both_genes() && chimeric_alignment->second.split_reads1 + chimeric_alignment->second.split_reads2 < min_support)
			chimeric_alignment->second.filter = FILTERS.at("min_support");
		else
			remaining++;
	}
	return remaining;
}
