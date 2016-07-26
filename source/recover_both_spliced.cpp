#include <map>
#include "common.hpp"
#include "annotation.hpp"
#include "recover_both_spliced.hpp"

using namespace std;

bool are_both_breakpoints_spliced(const fusion_t& fusion, const annotation_t& gene_annotation) {
	if (!fusion.filters.empty() &&
	    fusion.filters.find(FILTERS.at("promiscuous_genes")) == fusion.filters.end() &&
	    fusion.filters.find(FILTERS.at("min_support")) == fusion.filters.end())
		return false; // we won't recover fusions which were not discarded due to low support

	if (!fusion.spliced1 || !fusion.spliced2)
		return false; // both sides must be spliced

	if (gene_annotation[fusion.gene1].strand == gene_annotation[fusion.gene2].strand && fusion.direction1 == fusion.direction2 ||
	    gene_annotation[fusion.gene1].strand != gene_annotation[fusion.gene2].strand && fusion.direction1 != fusion.direction2)
		return false; // exons would be spliced in a nonsensical way

	return true;
}

unsigned int recover_both_spliced(fusions_t& fusions, const annotation_t& gene_annotation, const bool low_tumor_content) {

	// look for fusions with split reads where the split read and the clipped segment are both at exon boundaries
	map< pair<gene_t,gene_t>, unsigned int > fusions_with_both_spliced;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i)
		if (are_both_breakpoints_spliced(i->second, gene_annotation))
			fusions_with_both_spliced[make_pair(i->second.gene1, i->second.gene2)] += i->second.supporting_reads();

	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (i->second.filters.empty()) { // fusion has not been filtered, no need to recover
			remaining++;
			continue;
		}

		if (are_both_breakpoints_spliced(i->second, gene_annotation) &&
		    fusions_with_both_spliced[make_pair(i->second.gene1, i->second.gene2)] >= 2 || low_tumor_content) { // require at least two reads or else the false positive rate sky-rockets
			i->second.filters.clear();
			remaining++;
		}
	}
	return remaining;
}
