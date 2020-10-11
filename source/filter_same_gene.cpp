#include "common.hpp"
#include "annotation.hpp"
#include "filter_same_gene.hpp"

using namespace std;

unsigned int filter_same_gene(chimeric_alignments_t& chimeric_alignments, exon_annotation_index_t& exon_annotation_index) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != FILTER_none)
			continue; // the read has already been filtered

		// check if mate1 and mate2 map to the same gene
		gene_set_t common_genes;
		if (chimeric_alignment->second.size() == 2) // discordant mate
			combine_annotations(chimeric_alignment->second[MATE1].genes, chimeric_alignment->second[MATE2].genes, common_genes, false);
		else // split read
			combine_annotations(chimeric_alignment->second[MATE2].genes, chimeric_alignment->second[SUPPLEMENTARY].genes, common_genes, false);
		if (common_genes.empty()) {
			remaining++;
			continue; // we are only interested in intragenic events here
		}

		if (chimeric_alignment->second.size() == 2) { // discordant mates

			if (chimeric_alignment->second[MATE1].strand == FORWARD && chimeric_alignment->second[MATE2].strand == REVERSE && chimeric_alignment->second[MATE1].start <= chimeric_alignment->second[MATE2].end ||
			    chimeric_alignment->second[MATE1].strand == REVERSE && chimeric_alignment->second[MATE2].strand == FORWARD && chimeric_alignment->second[MATE1].end   >= chimeric_alignment->second[MATE2].start) {
				chimeric_alignment->second.filter = FILTER_same_gene; // normal alignment
				continue;
			}

		} else { // split read

			if (chimeric_alignment->second[SPLIT_READ].strand == FORWARD && chimeric_alignment->second[SUPPLEMENTARY].strand == FORWARD && chimeric_alignment->second[SPLIT_READ].start >= chimeric_alignment->second[SUPPLEMENTARY].end ||
			    chimeric_alignment->second[SPLIT_READ].strand == REVERSE && chimeric_alignment->second[SUPPLEMENTARY].strand == REVERSE && chimeric_alignment->second[SPLIT_READ].end   <= chimeric_alignment->second[SUPPLEMENTARY].start) {
				chimeric_alignment->second.filter = FILTER_same_gene; // normal alignment
				continue;
			}

		}

		remaining++;
	}
	return remaining;
}

