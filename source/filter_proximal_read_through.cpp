#include <string>
#include "common.hpp"
#include "annotation.hpp"
#include "filter_proximal_read_through.hpp"

using namespace std;

unsigned int filter_proximal_read_through(chimeric_alignments_t& chimeric_alignments, const int min_distance) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != FILTER_none)
			continue; // the read has already been filtered

		// find forward and reverse mate
		alignment_t* forward_mate;
		alignment_t* reverse_mate;
		if (chimeric_alignment->second.size() == 2) { // discordant mates
			forward_mate = &((chimeric_alignment->second[MATE1].strand == FORWARD) ? chimeric_alignment->second[MATE1] : chimeric_alignment->second[MATE2]);
			reverse_mate = &((chimeric_alignment->second[MATE1].strand == FORWARD) ? chimeric_alignment->second[MATE2] : chimeric_alignment->second[MATE1]);
		} else { // split read
			forward_mate = &((chimeric_alignment->second[SPLIT_READ].strand == FORWARD) ? chimeric_alignment->second[SUPPLEMENTARY] : chimeric_alignment->second[SPLIT_READ]);
			reverse_mate = &((chimeric_alignment->second[SPLIT_READ].strand == FORWARD) ? chimeric_alignment->second[SPLIT_READ] : chimeric_alignment->second[SUPPLEMENTARY]);
		}

		// only proper pairs can be read-through fragments
		if (chimeric_alignment->second.size() == 2 && forward_mate->strand != reverse_mate->strand && forward_mate->contig == reverse_mate->contig && forward_mate->end < reverse_mate->start ||
		    chimeric_alignment->second.size() == 3 && forward_mate->strand == reverse_mate->strand && forward_mate->contig == reverse_mate->contig && forward_mate->end < reverse_mate->start) {

			// find boundaries of biggest gene that the mates overlap with
			position_t forward_gene_start, forward_gene_end, reverse_gene_start, reverse_gene_end;
			get_boundaries_of_biggest_gene(forward_mate->genes, forward_gene_start, forward_gene_end);
			get_boundaries_of_biggest_gene(reverse_mate->genes, reverse_gene_start, reverse_gene_end);

			// remove chimeric alignment when mates map too close to end of gene
			if (forward_mate->end >= reverse_gene_start - min_distance || reverse_mate->start <= forward_gene_end + min_distance) {
				chimeric_alignment->second.filter = FILTER_read_through;
				continue;
			}
		}

		// we only get here, if the chimeric alignments were not filtered
		++remaining;
	}

	return remaining;
}

