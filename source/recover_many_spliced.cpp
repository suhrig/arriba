#include <tuple>
#include <unordered_map>
#include "common.hpp"
#include "recover_many_spliced.hpp"

using namespace std;

unsigned int recover_many_spliced(fusions_t& fusions, const unsigned int min_spliced_events) {

	// look for any spliced reads between two genes
	unordered_map< tuple<gene_t,gene_t>, set< tuple<position_t,position_t> > > spliced_fusions_by_gene_pair;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (!fusion->second.is_read_through() &&
		    (fusion->second.spliced1 || fusion->second.spliced2) &&
		    fusion->second.gene1 != fusion->second.gene2 &&
		    !fusion->second.breakpoint_overlaps_both_genes() &&
		    (fusion->second.filter == FILTER_none ||
		     fusion->second.filter == FILTER_inconsistently_clipped ||
		     fusion->second.filter == FILTER_relative_support ||
		     fusion->second.filter == FILTER_min_support ||
		     fusion->second.filter == FILTER_select_best)) {
			spliced_fusions_by_gene_pair[make_tuple(fusion->second.gene1, fusion->second.gene2)].insert(make_tuple(fusion->second.breakpoint1/10, fusion->second.breakpoint2/10)); // bin breakpoints a little to avoid counting misaligned splice sites twice
		}
	}

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter == FILTER_none) { // fusion has not been filtered, no need to recover
			remaining++;
			continue;
		}

		if (fusion->second.is_read_through() ||
		    fusion->second.gene1 == fusion->second.gene2 ||
		    fusion->second.breakpoint_overlaps_both_genes())
			continue; // don't recover events between partners which are likely to occur by chance

		if (fusion->second.filter == FILTER_inconsistently_clipped ||
                    fusion->second.filter == FILTER_relative_support ||
                    fusion->second.filter == FILTER_min_support ||
                    fusion->second.filter == FILTER_select_best) {
			if ((fusion->second.spliced1 || fusion->second.spliced2) &&
			    spliced_fusions_by_gene_pair[make_tuple(fusion->second.gene1, fusion->second.gene2)].size() >= min_spliced_events) {
				fusion->second.filter = FILTER_none;
				remaining++;
			}
		}
	}
	return remaining;
}
