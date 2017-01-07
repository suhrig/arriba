#include <vector>
#include "common.hpp"
#include "fusions.hpp"
#include "merge_adjacent_fusions.hpp"

using namespace std;

bool sort_fusions_by_coordinate(const fusion_t* x, const fusion_t* y) {
	if (x->contig1 != y->contig1)
		return x->contig1 < y->contig1;
	else if (x->breakpoint1 != y->breakpoint1)
		return x->breakpoint1 < y->breakpoint1;
	else if (x->contig2 != y->contig2)
		return x->contig2 < y->contig2;
	else
		return x->breakpoint2 < y->breakpoint2;
}

unsigned int merge_adjacent_fusions(fusions_t& fusions, const int max_distance) {

	vector<fusion_t*> sorted_fusions;
	sorted_fusions.reserve(fusions.size());
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter == NULL) // only merge fusions which have not been filtered before
			sorted_fusions.push_back(&fusion->second);
	sort(sorted_fusions.begin(), sorted_fusions.end(), sort_fusions_by_coordinate);

	for (auto fusion = sorted_fusions.begin(); fusion != sorted_fusions.end(); ++fusion) {

		if ((**fusion).split_reads1 + (**fusion).split_reads2 == 0)
			continue; // only merge fusions with exactly known breakpoints

		// find all adjacent breakpoints
		vector<fusion_t*> adjacent_fusions;

		// look upstream for mergeable breakpoints
		for (vector<fusion_t*>::reverse_iterator previous_fusion(fusion); previous_fusion != sorted_fusions.rend() && (**previous_fusion).contig1 == (**fusion).contig1 && (**previous_fusion).breakpoint1 >= (**fusion).breakpoint1-max_distance; ++previous_fusion) {
			if ((**previous_fusion).gene1 == (**fusion).gene1 &&
			    (**previous_fusion).gene2 == (**fusion).gene2 &&
			    (**previous_fusion).direction1 == (**fusion).direction1 &&
			    (**previous_fusion).direction2 == (**fusion).direction2 &&
			    (**previous_fusion).split_reads1 + (**previous_fusion).split_reads2 > 0 &&
			    (**previous_fusion).contig2 == (**fusion).contig2 &&
			    (**previous_fusion).breakpoint2 == (**fusion).breakpoint2 + ((**fusion).breakpoint1 - (**previous_fusion).breakpoint1) * (((**fusion).direction1 == (**fusion).direction2) ? +1 : -1)) {
				adjacent_fusions.push_back(*previous_fusion);
			}
		}

		// look downstream for mergeable breakpoints
		for (vector<fusion_t*>::iterator following_fusion(fusion); following_fusion != sorted_fusions.end() && (**following_fusion).contig1 == (**fusion).contig1 && (**following_fusion).breakpoint1 <= (**fusion).breakpoint1+max_distance; ++following_fusion) {
			if (following_fusion != fusion &&
			    (**following_fusion).gene1 == (**fusion).gene1 &&
			    (**following_fusion).gene2 == (**fusion).gene2 &&
			    (**following_fusion).direction1 == (**fusion).direction1 &&
			    (**following_fusion).direction2 == (**fusion).direction2 &&
			    (**following_fusion).split_reads1 + (**following_fusion).split_reads2 > 0 &&
			    (**following_fusion).contig2 == (**fusion).contig2 &&
			    (**following_fusion).breakpoint2 == (**fusion).breakpoint2 + ((**following_fusion).breakpoint1 - (**fusion).breakpoint1) * (((**fusion).direction1 == (**fusion).direction2) ? -1 : +1)) {
				adjacent_fusions.push_back(*following_fusion);
			}
		}

		// select the one with the most supporting alignments
		unsigned int sum_split_reads1 = 0, sum_split_reads2 = 0;
		bool fusion_has_most_support = true;
		for (unsigned int k = 0; k < adjacent_fusions.size(); ++k) {
			if ((**fusion).supporting_reads() < adjacent_fusions[k]->supporting_reads()) { // fewer supporting reads => this is not the best fusion
				fusion_has_most_support = false;
				break;
			} else {
				sum_split_reads1 += adjacent_fusions[k]->split_reads1;
				sum_split_reads2 += adjacent_fusions[k]->split_reads2;
			}
		}

		// merge, i.e. add split reads of adjacent breakpoits to best breakpoint and remove inferior ones
		if (fusion_has_most_support) {
			(**fusion).split_reads1 += sum_split_reads1;
			(**fusion).split_reads2 += sum_split_reads2;
			for (unsigned int k = 0; k < adjacent_fusions.size(); ++k)
				adjacent_fusions[k]->filter = FILTERS.at("merge_adjacent");
		}
	}

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter == NULL)
			remaining++;
	return remaining;
}

