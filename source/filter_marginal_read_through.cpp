#include "common.hpp"
#include "read_stats.hpp"
#include "filter_marginal_read_through.hpp"

using namespace std;

unsigned int filter_marginal_read_through(fusions_t& fusions, const coverage_t& coverage) {

	const float margin = 0.01; // fraction of gene to be considered the margins of the gene
	const float min_vaf = 0.07; // do not remove a fusion if the supporting reads make up more than this fraction of the coverage

	unsigned int remaining = 0;
	for (auto fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter == FILTER_none && fusion->second.is_read_through()) {

			// compute relative position of breakpoint in splice donor & acceptor
			// 0 means at the end of the donor / beginning of acceptor
			// 1 means at the beginning of the donor / end of acceptor
			// values close to 1 are very frequent and an indication of read-through rather than a real fusion caused by a genomic deletion
			double position_in_donor = 1;
			double position_in_acceptor = 1;
			if (!fusion->second.gene1->is_dummy && fusion->second.gene1->strand == FORWARD && fusion->second.direction1 == DOWNSTREAM) {
				position_in_donor = 1.0 * (fusion->second.breakpoint1 - fusion->second.gene1->start) / (fusion->second.gene1->end - fusion->second.gene1->start);
			} else if (!fusion->second.gene2->is_dummy && fusion->second.gene2->strand == REVERSE && fusion->second.direction2 == UPSTREAM) {
				position_in_donor = 1.0 * (fusion->second.gene2->end - fusion->second.breakpoint2) / (fusion->second.gene2->end - fusion->second.gene2->start);
			} else if (!fusion->second.gene1->is_dummy && fusion->second.gene1->strand == REVERSE && fusion->second.direction1 == DOWNSTREAM) {
				position_in_acceptor = 1.0 * (fusion->second.breakpoint1 - fusion->second.gene1->start) / (fusion->second.gene1->end - fusion->second.gene1->start);
			} else if (!fusion->second.gene2->is_dummy && fusion->second.gene2->strand == FORWARD && fusion->second.direction2 == UPSTREAM) {
				position_in_acceptor = 1.0 * (fusion->second.gene2->end - fusion->second.breakpoint2) / (fusion->second.gene2->end - fusion->second.gene2->start);
			} else // if we get here, both breakpoints are intergenic, in which case we don't apply this filter
				continue;

			// discard event, if both breakpoints are close to the boundaries of the fused genes and the supporting reads make up only a fraction of the coverage
			int coverage1 = coverage.get_coverage(fusion->second.contig1, fusion->second.breakpoint1, (fusion->second.direction1 == UPSTREAM) ? DOWNSTREAM : UPSTREAM);
			int coverage2 = coverage.get_coverage(fusion->second.contig2, fusion->second.breakpoint2, (fusion->second.direction2 == UPSTREAM) ? DOWNSTREAM : UPSTREAM);
			if (position_in_donor > 1-margin && position_in_acceptor > 1-margin && fusion->second.supporting_reads() < min_vaf * max(coverage1, coverage2))
				fusion->second.filter = FILTER_marginal_read_through;
		}

		if (fusion->second.filter == FILTER_none)
			++remaining;
	}

	return remaining;
}

