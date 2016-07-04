#include "annotation.hpp"
#include "common.hpp"
#include "filter_misspliced.hpp"

using namespace std;

unsigned int filter_misspliced(fusions_t& fusions, const annotation_t& gene_annotation) {
	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
		if (!i->second.filters.empty())
			continue; // fusion has already been filtered

		// discard fusions which should be spliced, but aren't
		// (i.e., two genes are fused intronically in the right orientation, but the transcript is not spliced)
		if ((i->second.split_reads1 > 0 && i->second.split_reads2 + i->second.discordant_mates == 0 || i->second.split_reads2 > 0 && i->second.split_reads1 + i->second.discordant_mates == 0) && // only breakpoints with split reads can be spliced
		    !gene_annotation[i->second.gene1].is_dummy && !gene_annotation[i->second.gene2].is_dummy && // ignore intergenic regions (dummy genes)
		    (i->second.spliced1 && !i->second.spliced2 && !i->second.exonic2 && // breakpoint1 is spliced, but breakpoint2 isn't
		     (i->second.direction1 == UPSTREAM && gene_annotation[i->second.gene1].strand == REVERSE || i->second.direction1 == DOWNSTREAM && gene_annotation[i->second.gene1].strand == FORWARD) && // gene1 is the 5' end of the transcript
		     (i->second.direction2 == UPSTREAM && gene_annotation[i->second.gene2].strand == FORWARD || i->second.direction2 == DOWNSTREAM && gene_annotation[i->second.gene2].strand == REVERSE) || // gene2 is the 3' end of the transcript
		     i->second.spliced2 && !i->second.spliced1 && !i->second.exonic1 &&
		     (i->second.direction2 == UPSTREAM && gene_annotation[i->second.gene2].strand == REVERSE || i->second.direction2 == DOWNSTREAM && gene_annotation[i->second.gene2].strand == FORWARD) && // gene2 is the 5' end of the transcript
		     (i->second.direction1 == UPSTREAM && gene_annotation[i->second.gene1].strand == FORWARD || i->second.direction1 == DOWNSTREAM && gene_annotation[i->second.gene1].strand == REVERSE))) { // gene1 is the 3' end of the transcript
			i->second.filters.insert(FILTERS.at("misspliced"));
		} else {
			remaining++;
		}
	}
	return remaining;
}

