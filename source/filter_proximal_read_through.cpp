#include <string>
#include "common.hpp"
#include "annotation.hpp"
#include "filter_proximal_read_through.hpp"

using namespace std;

unsigned int filter_proximal_read_through(chimeric_alignments_t& chimeric_alignments, const unsigned int min_distance) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {

		if (!i->second.filters.empty())
			continue; // the read has already been filtered

		// find forward and reverse mate
		alignment_t* forward_mate;
		alignment_t* reverse_mate;
		if (i->second.size() == 2) { // discordant mates
			forward_mate = &((i->second[MATE1].strand == FORWARD) ? i->second[MATE1] : i->second[MATE2]);
			reverse_mate = &((i->second[MATE1].strand == FORWARD) ? i->second[MATE2] : i->second[MATE1]);
		} else { // split read
			forward_mate = &((i->second[SPLIT_READ].strand == FORWARD) ? i->second[SUPPLEMENTARY] : i->second[SPLIT_READ]);
			reverse_mate = &((i->second[SPLIT_READ].strand == FORWARD) ? i->second[SPLIT_READ] : i->second[SUPPLEMENTARY]);
		}

		// only proper pairs can be read-through fragments
		if (i->second.size() == 2 && forward_mate->strand != reverse_mate->strand && forward_mate->contig == reverse_mate->contig && forward_mate->end < reverse_mate->start ||
		    i->second.size() == 3 && forward_mate->strand == reverse_mate->strand && forward_mate->contig == reverse_mate->contig && forward_mate->end < reverse_mate->start) {

			// find boundaries of biggest gene that the mates overlap with
			position_t forward_gene_start, forward_gene_end, reverse_gene_start, reverse_gene_end;
			get_boundaries_of_biggest_gene(forward_mate->genes, forward_gene_start, forward_gene_end);
			get_boundaries_of_biggest_gene(reverse_mate->genes, reverse_gene_start, reverse_gene_end);

			// remove chimeric alignment when mates map too close to end of gene
			if (forward_mate->end >= reverse_gene_start - min_distance || reverse_mate->start <= forward_gene_end + min_distance) {
				i->second.filters.insert(FILTERS.at("read_through"));
				continue;
			}
		}

		// we only get here, if the chimeric alignments were not filtered
		++remaining;
	}

	return remaining;
}

