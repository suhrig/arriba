#include <map>
#include <tuple>
#include "common.hpp"
#include "recover_both_spliced.hpp"

using namespace std;

unsigned int recover_both_spliced(fusions_t& fusions, const bool low_tumor_content) {

	// look for any supporting reads between two genes
	map< tuple<gene_t,gene_t,direction_t,direction_t>, vector<fusion_t*> > fusions_by_gene_pair;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i)
		if (i->second.filters.empty() ||
		    i->second.filters.find(FILTERS.at("intronic")) != i->second.filters.end() ||
		    i->second.filters.find(FILTERS.at("promiscuous_genes")) != i->second.filters.end() ||
		    i->second.filters.find(FILTERS.at("min_support")) != i->second.filters.end()) {
			fusions_by_gene_pair[make_tuple(i->second.gene1, i->second.gene2, i->second.direction1, i->second.direction2)].push_back(&i->second);
		}

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filters.empty()) { // fusion has not been filtered, no need to recover
			remaining++;
			continue;
		}

		if (fusion->second.gene1 == fusion->second.gene2)
			continue; // don't recover intragenic events (this would produce too many hits)

		if (!fusion->second.filters.empty() &&
		    fusion->second.filters.find(FILTERS.at("promiscuous_genes")) == fusion->second.filters.end() &&
		    fusion->second.filters.find(FILTERS.at("min_support")) == fusion->second.filters.end())
			continue; // we won't recover fusions which were not discarded due to low support

		if (!fusion->second.spliced1 || !fusion->second.spliced2)
			continue; // both sides must be spliced

		if (fusion->second.gene1->strand == fusion->second.gene2->strand && fusion->second.direction1 == fusion->second.direction2 ||
		    fusion->second.gene1->strand != fusion->second.gene2->strand && fusion->second.direction1 != fusion->second.direction2)
			continue; // exons would be spliced in a nonsensical way

		// count all supporting reads of all fusions between the pair of genes
		auto fusions_of_given_gene_pair = fusions_by_gene_pair.find(make_tuple(fusion->second.gene1, fusion->second.gene2, fusion->second.direction1, fusion->second.direction2));
		if (fusions_of_given_gene_pair != fusions_by_gene_pair.end()) {
			unsigned int sum_of_supporting_reads = 0;
			for (auto another_fusion = fusions_of_given_gene_pair->second.begin(); another_fusion != fusions_of_given_gene_pair->second.end(); ++another_fusion) {
				// for read-through fusions, we require the distance of the other fusion to be at least as big as that of the fusion with spliced breakpoints
				if (!(**another_fusion).is_read_through() ||
				    (**another_fusion).breakpoint2 - (**another_fusion).breakpoint1 >= fusion->second.breakpoint2 - fusion->second.breakpoint1)
					sum_of_supporting_reads += max((unsigned int) 1, (**another_fusion).supporting_reads());
			}

			if (sum_of_supporting_reads >= 2 || low_tumor_content) { // require at least two reads or else the false positive rate sky-rockets
				fusion->second.filters.clear();
				remaining++;
			}
		}
	}
	return remaining;
}
