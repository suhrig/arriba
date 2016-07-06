#include <cmath>
#include "common.hpp"
#include "annotation.hpp"
#include "recover_both_spliced.hpp"

using namespace std;

unsigned int recover_both_spliced(fusions_t& fusions, const annotation_t& gene_annotation) {

	// look for fusions with split reads where the split read and the clipped segment are both at exon boundaries
	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (i->second.filters.empty()) { // fusion has not been filtered, no need to recover
			remaining++;
			continue;
		}

		if (!i->second.filters.empty() && // fusion has been filtered
		    i->second.filters.find(FILTERS.at("promiscuous_genes")) == i->second.filters.end() && i->second.filters.find(FILTERS.at("min_support")) == i->second.filters.end()) // reason is not low support
			continue; // we won't recover fusions which were not discarded due to low support

		if (i->second.spliced1 && i->second.spliced2 && // check if both sides are spliced
		    i->second.supporting_reads() >= 2 && // require at least two reads or else the false positive rate sky-rockets
		    (gene_annotation[i->second.gene1].strand == gene_annotation[i->second.gene2].strand && i->second.direction1 != i->second.direction2 ||
		     gene_annotation[i->second.gene1].strand != gene_annotation[i->second.gene2].strand && i->second.direction1 == i->second.direction2)) { // make sure exons are joined in sensible orientation
			i->second.filters.clear();
			remaining++;
		}
	}
	return remaining;
}
