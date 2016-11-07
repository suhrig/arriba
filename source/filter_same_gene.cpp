#include "common.hpp"
#include "annotation.hpp"
#include "filter_same_gene.hpp"

using namespace std;

unsigned int filter_same_gene(chimeric_alignments_t& chimeric_alignments, exon_annotation_index_t& exon_annotation_index) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {

		if (!i->second.filters.empty())
			continue; // the read has already been filtered

		// check if mate1 and mate2 map to the same gene
		gene_set_t common_genes;
		if (i->second.size() == 2) // discordant mate
			combine_annotations(i->second[MATE1].genes, i->second[MATE2].genes, common_genes, false);
		else // split read
			combine_annotations(i->second[MATE2].genes, i->second[SUPPLEMENTARY].genes, common_genes, false);
		if (common_genes.empty()) {
			remaining++;
			continue; // we are only interested in intragenic events here
		}

		if (i->second.size() == 2) { // discordant mates

			if (i->second[MATE1].strand == FORWARD && i->second[MATE2].strand == REVERSE && i->second[MATE1].start <= i->second[MATE2].end ||
			    i->second[MATE1].strand == REVERSE && i->second[MATE2].strand == FORWARD && i->second[MATE1].end   >= i->second[MATE2].start) {
				i->second.filters.insert(FILTERS.at("same_gene")); // normal alignment
				continue;
			}

		} else { // split read

			if (i->second[SPLIT_READ].strand == FORWARD && i->second[SUPPLEMENTARY].strand == FORWARD && i->second[SPLIT_READ].start >= i->second[SUPPLEMENTARY].end ||
			    i->second[SPLIT_READ].strand == REVERSE && i->second[SUPPLEMENTARY].strand == REVERSE && i->second[SPLIT_READ].end   <= i->second[SUPPLEMENTARY].start) {
				i->second.filters.insert(FILTERS.at("same_gene")); // normal alignment
				continue;
			}

		}

		remaining++;
	}
	return remaining;
}

