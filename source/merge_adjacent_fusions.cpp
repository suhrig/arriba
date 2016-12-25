#include <tuple>
#include <vector>
#include "common.hpp"
#include "options.hpp"
#include "fusions.hpp"
#include "merge_adjacent_fusions.hpp"

using namespace std;

unsigned int merge_adjacent_fusions(fusions_t& fusions, const int max_distance) {
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != NULL)
			continue; // fusion has already been filtered

		if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0)
			continue; // only merge fusions with exactly known breakpoints

		// find all adjacent breakpoints
		vector<fusions_t::iterator> adjacent_fusions;
		for (int distance = -max_distance; distance <= +max_distance; distance++)
			if (distance != 0) { // don't merge with itself
				fusions_t::iterator adjacent_fusion = fusions.find(
					make_tuple(
						fusion->second.gene1,
						fusion->second.gene2,
						fusion->second.contig1,
						fusion->second.contig2,
						fusion->second.breakpoint1 + distance,
						fusion->second.breakpoint2 + distance * ((fusion->second.direction1 == fusion->second.direction2) ? -1 : +1),
						fusion->second.direction1,
						fusion->second.direction2
					)
				);
				if (adjacent_fusion != fusions.end() && adjacent_fusion->second.filter == NULL &&
				    adjacent_fusion->second.split_reads1 + adjacent_fusion->second.split_reads2 > 0)
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
					adjacent_fusions[k]->second.filter = FILTERS.at("merge_adjacent");
		}
	}

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter == NULL)
			remaining++;
	return remaining;
}

