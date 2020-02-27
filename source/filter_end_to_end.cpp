#include "common.hpp"
#include "filter_end_to_end.hpp"

using namespace std;

// slide over gene and count bases that do not belong to any exon of the given gene (= intronic bases)
float calculate_intronic_fraction(const gene_t gene, const exon_annotation_index_t& exon_annotation_index) {
	unsigned int intronic_bases = 0;
	position_t previous_position = gene->start;
	for (auto position = exon_annotation_index[gene->contig].lower_bound(gene->start); position != exon_annotation_index[gene->contig].end() && position->first <= gene->end; ++position) {
		for (auto exon = position->second.begin(); exon != position->second.end(); ++exon) {
			if ((**exon).gene == gene) {
				if (previous_position < (**exon).start)
					intronic_bases += (**exon).start - previous_position;
				if (previous_position < (**exon).end)
					previous_position = (**exon).end + 1;
				break;
			}
		}
	}
	return ((float) intronic_bases) / (gene->end - gene->start + 1);
}

unsigned int filter_end_to_end_fusions(fusions_t& fusions, const exon_annotation_index_t& exon_annotation_index) {

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != FILTER_none)
			continue; // fusion has already been filtered

		if (!fusion->second.is_read_through() &&
		    fusion->second.gene1 != fusion->second.gene2 &&
		    (fusion->second.spliced1 || fusion->second.spliced2)) { // spliced breakpoints are likely true
			++remaining;
			continue;
		}

		if (fusion->second.discordant_mates + fusion->second.split_reads1 == 0 || // only filter breakpoints with low support
		    fusion->second.discordant_mates + fusion->second.split_reads2 == 0 ||
		    fusion->second.split_reads1 + fusion->second.split_reads2 == 0 ||
		    fusion->second.breakpoint_overlaps_both_genes() && (fusion->second.split_reads1 == 0 || fusion->second.split_reads2 == 0)) {

			if ((fusion->second.gene1->is_dummy || (fusion->second.gene1->strand == FORWARD && fusion->second.direction1 == UPSTREAM) || (fusion->second.gene1->strand == REVERSE && fusion->second.direction1 == DOWNSTREAM)) &&
			    (fusion->second.gene2->is_dummy || (fusion->second.gene2->strand == FORWARD && fusion->second.direction2 == UPSTREAM) || (fusion->second.gene2->strand == REVERSE && fusion->second.direction2 == DOWNSTREAM))) {

				// translocations involving the IG and TCR loci often only have discordant mates,
				// because STAR fails to align the supporting split reads
				// => if there are many discordant mates, we spare the event,
				//    unless there is additional evidence that raises scepticism, namely:
				//    - the breakpoints are close to each other
				//    - the breakpoints are located in exons, although the genes are mostly made up of introns
				const unsigned int many_discordant_mates = 10;
				const unsigned int min_breakpoint_distance = 1000000;
				const float max_intronic_fraction = 0.66;
				if (fusion->second.discordant_mates < many_discordant_mates ||
				    fusion->second.contig1 == fusion->second.contig2 && abs(fusion->second.breakpoint1 - fusion->second.breakpoint2) < min_breakpoint_distance ||
				    fusion->second.exonic1 && fusion->second.exonic2 &&
				    calculate_intronic_fraction(fusion->second.gene1, exon_annotation_index) > max_intronic_fraction &&
				    calculate_intronic_fraction(fusion->second.gene2, exon_annotation_index) > max_intronic_fraction) {

					fusion->second.filter = FILTER_end_to_end;
					continue;

				}
			}
		}

		++remaining;
	}

	return remaining;
}

