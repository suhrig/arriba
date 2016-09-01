#include <map>
#include <tuple>
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
//TODO ignore read-through fusions?

	// look for any supporting reads between two genes
	map< tuple<gene_t,gene_t,direction_t,direction_t>, unsigned int > supporting_reads_by_gene_pair;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i)
		if (i->second.filters.empty() ||
		    i->second.filters.find(FILTERS.at("intronic")) != i->second.filters.end() ||
		    i->second.filters.find(FILTERS.at("promiscuous_genes")) != i->second.filters.end() ||
		    i->second.filters.find(FILTERS.at("min_support")) != i->second.filters.end()) {
			supporting_reads_by_gene_pair[make_tuple(i->second.gene1, i->second.gene2, i->second.direction1, i->second.direction2)] += max(i->second.supporting_reads(), (unsigned int) 1);
		}

	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (i->second.filters.empty()) { // fusion has not been filtered, no need to recover
			remaining++;
			continue;
		}

		if (are_both_breakpoints_spliced(i->second, gene_annotation) &&
		    supporting_reads_by_gene_pair[make_tuple(i->second.gene1, i->second.gene2, i->second.direction1, i->second.direction2)] >= 2 || low_tumor_content) { // require at least two reads or else the false positive rate sky-rockets
			i->second.filters.clear();
			remaining++;
		}
	}
	return remaining;
}
