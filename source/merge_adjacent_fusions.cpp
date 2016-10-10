#include <tuple>
#include <vector>
#include "common.hpp"
#include "options.hpp"
#include "fusions.hpp"
#include "merge_adjacent_fusions.hpp"

using namespace std;

unsigned int merge_adjacent_fusions(fusions_t& fusions, const int max_distance) {
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (!fusion->second.filters.empty())
			continue; // fusion has already been filtered

		// find all adjacent breakpoints
		vector<fusions_t::iterator> adjacent_fusions;
		for (int d1 = -max_distance; d1 <= +max_distance; d1++)
			for (int d2 = -max_distance; d2 <= +max_distance; d2++)
				if (d1 != 0 || d2 != 0) { // don't merge with itself
					fusions_t::iterator adjacent_fusion = fusions.find(make_tuple(get<0>(fusion->first), get<1>(fusion->first), get<2>(fusion->first), get<3>(fusion->first), get<4>(fusion->first)+d1, get<5>(fusion->first)+d2, get<6>(fusion->first), get<7>(fusion->first)));
					if (adjacent_fusion != fusions.end() && adjacent_fusion->second.filters.empty())
						adjacent_fusions.push_back(adjacent_fusion);
				}

		// select the one with the most supporting alignments
		unsigned int sum_split_reads1 = 0, sum_split_reads2 = 0;
		bool fusion_has_most_support = true;
		for (unsigned int k = 0; k < adjacent_fusions.size(); ++k) {
			if (adjacent_fusions[k] != fusions.end()) {
				if (fusion->second.supporting_reads() < adjacent_fusions[k]->second.supporting_reads() || // fewer supporting reads => this is not the best fusion
				    fusion->second.supporting_reads() == adjacent_fusions[k]->second.supporting_reads() && // this condition does nothing but ensure deterministic behavior
				    (fusion->second.breakpoint1 < adjacent_fusions[k]->second.breakpoint1 || // this condition does nothing but ensure deterministic behavior
				     fusion->second.breakpoint1 == adjacent_fusions[k]->second.breakpoint1 && fusion->second.breakpoint2 < adjacent_fusions[k]->second.breakpoint2)) { // this condition does nothing but ensure deterministic behavior
					fusion_has_most_support = false;
					break;
				} else {
					sum_split_reads1 += adjacent_fusions[k]->second.split_reads1;
					sum_split_reads2 += adjacent_fusions[k]->second.split_reads2;
				}
			}
		}

		// merge, i.e. add split reads of adjacent breakpoits to best breakpoint and remove inferior ones
		if (fusion_has_most_support) {
			fusion->second.split_reads1 += sum_split_reads1;
			fusion->second.split_reads2 += sum_split_reads2;
			for (unsigned int k = 0; k < adjacent_fusions.size(); ++k)
				if (adjacent_fusions[k] != fusions.end())
					adjacent_fusions[k]->second.filters.insert(FILTERS.at("merge_adjacent"));
		}
	}

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filters.empty())
			remaining++;
	return remaining;
}

