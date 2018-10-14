#include "common.hpp"
#include "annotation.hpp"
#include "filter_intragenic_both_exonic.hpp"

using namespace std;

unsigned int filter_intragenic_both_exonic(fusions_t& fusions, const exon_annotation_index_t& exon_annotation_index, const float exonic_fraction) {
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter != NULL)
			continue; // read has already been filtered

		if ((fusion->second.breakpoint_overlaps_both_genes() || fusion->second.gene1 == fusion->second.gene2) &&
		    fusion->second.exonic1 && fusion->second.exonic2 &&
		    !(fusion->second.spliced1 && fusion->second.spliced2)) {
			// if less then 20% of the region between the breakpoints is exonic,
			// but the breakpoints are both exonic nonetheless, discard the event,
			// because it is unlikely that the breakpoints of a structural variant
			// both fall inside (the relatively small) exons rather than into an intron;
			// moreover, we discard the special case where there are no introns
			// between the breakpoints, not because this is unlikely, but because
			// this is the most frequent type of false positive
			int spliced_distance = get_spliced_distance(fusion->second.contig1, fusion->second.breakpoint1, fusion->second.breakpoint2, fusion->second.direction1, fusion->second.direction2, fusion->second.gene1, exon_annotation_index);
			int distance = fusion->second.breakpoint2 - fusion->second.breakpoint1;
			if (spliced_distance == distance || spliced_distance / distance < exonic_fraction) {
				fusion->second.filter = FILTERS.at("intragenic_exonic");
				continue;
			}
		}

		++remaining;
	}

	return remaining;
}

