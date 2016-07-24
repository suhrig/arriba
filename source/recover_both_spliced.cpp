#include <map>
#include "common.hpp"
#include "annotation.hpp"
#include "recover_both_spliced.hpp"

using namespace std;

unsigned int recover_both_spliced(fusions_t& fusions, const annotation_t& gene_annotation, const bool low_tumor_content) {

	// look for fusions with split reads where the split read and the clipped segment are both at exon boundaries
	map< pair<gene_t,gene_t>, unsigned int > fusions_with_both_spliced;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
		if (i->second.spliced1 && i->second.spliced2 && // check if both sides are spliced
		    (gene_annotation[i->second.gene1].strand == gene_annotation[i->second.gene2].strand && i->second.direction1 != i->second.direction2 ||
		     gene_annotation[i->second.gene1].strand != gene_annotation[i->second.gene2].strand && i->second.direction1 == i->second.direction2)) { // make sure exons are joined in sensible orientation
			fusions_with_both_spliced[make_pair(i->second.gene1, i->second.gene2)] += i->second.supporting_reads();
		}
	}

	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (i->second.filters.empty()) { // fusion has not been filtered, no need to recover
			remaining++;
			continue;
		}

		if (i->second.filters.find(FILTERS.at("promiscuous_genes")) == i->second.filters.end() &&
		    i->second.filters.find(FILTERS.at("min_support")) == i->second.filters.end())
			continue; // we won't recover fusions which were not discarded due to low support

		if (!(i->second.contig1 == i->second.contig2 && i->second.breakpoint2 - i->second.breakpoint1 < 400000 && i->second.direction1 == DOWNSTREAM && i->second.direction2 == UPSTREAM) && // ignore read-through fusions
		    fusions_with_both_spliced[make_pair(i->second.gene1, i->second.gene2)] >= 2 || low_tumor_content) { // require at least two reads or else the false positive rate sky-rockets
			i->second.filters.clear();
			remaining++;
		}
	}
	return remaining;
}
