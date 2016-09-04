#include <string>
#include "common.hpp"
#include "annotation.hpp"
#include "filter_proximal_read_through.hpp"

using namespace std;

unsigned int filter_proximal_read_through(chimeric_alignments_t& chimeric_alignments, annotation_t& gene_annotation, const unsigned int min_distance) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {

		if (!i->second.filters.empty())
			continue; // the read has already been filtered

		if (i->second.size() == 2) { // only filter discordant mates

			// find forward and reverse mate
			alignment_t& forward_mate = (i->second[MATE1].strand == FORWARD) ? i->second[MATE1] : i->second[MATE2];
			alignment_t& reverse_mate = (i->second[MATE1].strand == FORWARD) ? i->second[MATE2] : i->second[MATE1];

			// only proper pairs can be read-through fragments
			if (forward_mate.strand != reverse_mate.strand && forward_mate.contig == reverse_mate.contig && forward_mate.end < reverse_mate.start) {

				// find boundaries of biggest gene that the mates overlap with
				position_t forward_gene_start, forward_gene_end, reverse_gene_start, reverse_gene_end;
				get_boundaries_of_biggest_gene(forward_mate.genes, gene_annotation, forward_gene_start, forward_gene_end);
				get_boundaries_of_biggest_gene(reverse_mate.genes, gene_annotation, reverse_gene_start, reverse_gene_end);

				if ((!gene_annotation[*forward_mate.genes.begin()].is_dummy && forward_mate.end   >= reverse_gene_start - min_distance) ||
				    (!gene_annotation[*reverse_mate.genes.begin()].is_dummy && reverse_mate.start <= forward_gene_end   + min_distance)) {
					// delete the chimeric alignments, if one mate maps to a dummy gene
					// and the mate is proximal (< min_distance) to the gene of the other mate
					i->second.filters.insert(FILTERS.at("read_through"));
					continue;
				}
			}

		}

		// we only get here, if the chimeric alignments were not filtered
		++remaining;
	}

	return remaining;
}

