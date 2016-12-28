#include <map>
#include <tuple>
#include "common.hpp"
#include "recover_isoforms.hpp"

using namespace std;

unsigned int recover_isoforms(fusions_t& fusions) {

	// make a lookup table for all the genes that have fusions that passed all filters
	map< tuple<gene_t,gene_t,direction_t,direction_t>, bool > fused_gene_pairs;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter == NULL)
			fused_gene_pairs[make_tuple(fusion->second.gene1, fusion->second.gene2, fusion->second.direction1, fusion->second.direction2)] = true;

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter == NULL) { // fusion has not been filtered, no need to recover
			remaining++;
			continue;
		}

		if (fusion->second.filter == FILTERS.at("merge_adjacent") || // don't recover alternative alignments
		    fusion->second.filter == FILTERS.at("blacklist") || // don't recover normal splice variants and artifacts
		    fusion->second.filter == FILTERS.at("end_to_end")) // don't recover alignments that happen to end at splice-sites
			continue;

		// find all splice-variants
		if (fusion->second.spliced1 && fusion->second.spliced2 &&
	            fused_gene_pairs.find(make_tuple(fusion->second.gene1, fusion->second.gene2, fusion->second.direction1, fusion->second.direction2)) != fused_gene_pairs.end()) {
			fusion->second.filter = NULL;
			remaining++;
		}

	}

	return remaining;
}
