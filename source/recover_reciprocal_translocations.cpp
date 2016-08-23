#include <map>
#include <tuple>
#include "common.hpp"
#include "annotation.hpp"
#include "recover_reciprocal_translocations.hpp"

#include <iostream>

using namespace std;

unsigned int recover_reciprocal(fusions_t& fusions, const annotation_t& gene_annotation) {

	// list all fusions for all combinations of genes and direction
	map< tuple<gene_t,gene_t,direction_t,direction_t>, vector<fusions_t::iterator> > fusions_by_gene_pair;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i)
		if (i->second.gene1 != i->second.gene2 && // we cannot detect reciprocal swaps between the same gene (it would appear like splicing)
		    i->second.filters.find(FILTERS.at("promiscuous_genes")) != i->second.filters.end() || i->second.filters.find(FILTERS.at("min_support")) != i->second.filters.end()) // only recover fusions that were discarded due to low support
			fusions_by_gene_pair[make_tuple(i->second.gene1, i->second.gene2, i->second.direction1, i->second.direction2)].push_back(i);

	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (i->second.filters.empty()) { // fusion has not been filtered, no need to recover
			remaining++;
			continue;
		}

		if (i->second.filters.find(FILTERS.at("promiscuous_genes")) == i->second.filters.end() && i->second.filters.find(FILTERS.at("min_support")) == i->second.filters.end() || // only recover fusions that were discarded due to low support
		    gene_annotation[i->second.gene1].is_dummy || gene_annotation[i->second.gene2].is_dummy || // don't recover intergenic fusions
		    i->second.gene1 == i->second.gene2) // we cannot detect reciprocal swaps between the same gene (it would appear like splicing)
			continue;

		if ((i->second.direction1 == i->second.direction2 && gene_annotation[i->second.gene1].strand != gene_annotation[i->second.gene2].strand ||
		     i->second.direction1 != i->second.direction2 && gene_annotation[i->second.gene1].strand == gene_annotation[i->second.gene2].strand)) {

			// look for two fusions where the 5' gene and 3' gene are swapped
			vector<fusions_t::iterator> fusions_by_given_gene_pair = fusions_by_gene_pair[make_tuple(i->second.gene1, i->second.gene2,
			                                                                                         (i->second.direction1 == UPSTREAM) ? DOWNSTREAM : UPSTREAM,
			                                                                                         (i->second.direction2 == UPSTREAM) ? DOWNSTREAM : UPSTREAM)];
			for (auto j = fusions_by_given_gene_pair.begin(); j != fusions_by_given_gene_pair.end(); ++j) {
				if ((i->second.direction1 == DOWNSTREAM && i->second.breakpoint1 <= (**j).second.breakpoint1 || i->second.direction1 == UPSTREAM && i->second.breakpoint1 >= (**j).second.breakpoint1) &&
				    (i->second.direction2 == DOWNSTREAM && i->second.breakpoint2 <= (**j).second.breakpoint2 || i->second.direction2 == UPSTREAM && i->second.breakpoint2 >= (**j).second.breakpoint2) &&
				    (i->second.spliced1 || i->second.spliced2) &&
				    ((**j).second.spliced1 || (**j).second.spliced2)) {
					i->second.filters.clear();
					remaining++;
				}
			}
		}
	}
	return remaining;
}
