#include <unordered_map>
#include <tuple>
#include "common.hpp"
#include "fusions.hpp"
#include "select_best.hpp"

using namespace std;

unsigned int rank_fusion(fusion_t& fusion) {
	if (fusion.split_reads1 != 0 && fusion.split_reads2 != 0 && fusion.discordant_mates != 0) {
		return 4; // fusions with split reads in both genes and discordant mates get the best rank
	} else if (fusion.split_reads1 != 0 && fusion.split_reads2 != 0) {
		return 3; // fusions with only split reads (but no discordant mates) are usually better than only discordant mates, because they are more accurate
	} else if ((fusion.split_reads1 != 0 || fusion.split_reads2 != 0) && fusion.discordant_mates != 0) {
		return 2; // the next most accurate fusion is one with discordant mates and split reads in one gene
	} else if (fusion.split_reads1 != 0 || fusion.split_reads2 != 0) {
		return 1; // the next best fusion is one with split reads in only one gene
	} else { // only discordant mates
		return 0; // the worst fusions are those with only discordant mates, because the precise breakpoint is unknown
	}
}

unsigned int select_most_supported_breakpoints(fusions_t& fusions) {

	typedef tuple<gene_t /*gene1*/, gene_t /*gene2*/, direction_t /*direction1*/, direction_t /*direction2*/> gene_pair_t;
	unordered_map< gene_pair_t, fusions_t::iterator > best_breakpoints;

	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (i->second.filter != NULL)
			continue; // fusion has already been filtered

		gene_pair_t gene_pair = make_tuple(i->second.gene1, i->second.gene2, i->second.direction1, i->second.direction2);

		// look for fusion with most support
		if (best_breakpoints.find(gene_pair) == best_breakpoints.end()) {
			best_breakpoints[gene_pair] = i; // initialize; this is the first fusion of the given gene pair which we encountered
		} else {
			fusions_t::iterator current_best = best_breakpoints[gene_pair];
			if (rank_fusion(i->second) > rank_fusion(current_best->second)) { // preferentially look for breakpoints supported by split reads
				best_breakpoints[gene_pair] = i;
			} else if (rank_fusion(i->second) == rank_fusion(current_best->second)) { // then look for the breakpoints with most supporting reads
				if (i->second.supporting_reads() > current_best->second.supporting_reads()) {
					best_breakpoints[gene_pair] = i;
				} else if (i->second.supporting_reads() == current_best->second.supporting_reads()) { // then look for the most upstream / downstream breakpoints
					if (i->second.direction1 == DOWNSTREAM && i->second.breakpoint1 > current_best->second.breakpoint1 ||
					    i->second.direction1 == UPSTREAM   && i->second.breakpoint1 < current_best->second.breakpoint1) {
						best_breakpoints[gene_pair] = i;
					} else if (i->second.direction1 == DOWNSTREAM && i->second.breakpoint1 == current_best->second.breakpoint1 ||
					           i->second.direction1 == UPSTREAM   && i->second.breakpoint1 == current_best->second.breakpoint1) {
						if (i->second.direction2 == DOWNSTREAM && i->second.breakpoint2 > current_best->second.breakpoint2 ||
					            i->second.direction2 == UPSTREAM   && i->second.breakpoint2 < current_best->second.breakpoint2)
							best_breakpoints[gene_pair] = i;
					}
				}
			}
		}

	}

	// delete all fusions but the best ones
	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (i->second.filter != NULL)
			continue; // the fusion has already been filtered

		if (i == best_breakpoints[make_tuple(i->second.gene1, i->second.gene2, i->second.direction1, i->second.direction2)])
			remaining++;
		else
			i->second.filter = FILTERS.at("select_best");

	}
	return remaining;
}

