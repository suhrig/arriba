#include "common.hpp"
#include "annotation.hpp"
#include "filter_same_gene.hpp"

using namespace std;




unsigned int filter_same_gene(chimeric_alignments_t& chimeric_alignments) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {
		if (!i->second.filters.empty())
			continue; // the read has already been filtered

		unsigned int mate1, mate2;
		if (i->second.size() == 2) { // discordant mates
			mate1 = MATE1; mate2 = MATE2;
		} else { // split read
			mate1 = SPLIT_READ; mate2 = SUPPLEMENTARY;
		}

		// check if mate1 and mate2 map to the same gene(s)
		// if so, remove them
		gene_set_t common_genes;
		combine_annotations(i->second[mate1].genes, i->second[mate2].genes, common_genes, false);
		if (!common_genes.empty())
			i->second.filters.insert(FILTERS.at("same_gene"));
		else
			++remaining;
	}
	return remaining;
}

