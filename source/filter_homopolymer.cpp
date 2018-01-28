#include "common.hpp"
#include "annotation.hpp"
#include "filter_homopolymer.hpp"

using namespace std;

bool is_split_read_spliced(const alignment_t& split_read, const exon_annotation_index_t& exon_annotation_index) {
	direction_t direction = (split_read.strand == FORWARD) ? UPSTREAM : DOWNSTREAM;
	position_t breakpoint = (split_read.strand == FORWARD) ? split_read.start : split_read.end;
	for (gene_set_t::const_iterator gene = split_read.genes.begin(); gene != split_read.genes.end(); ++gene)
		if (is_breakpoint_spliced(*gene, direction, breakpoint, exon_annotation_index))
			return true;
	return false;
}

unsigned int filter_homopolymer(chimeric_alignments_t& chimeric_alignments, const unsigned int homopolymer_length, const exon_annotation_index_t& exon_annotation_index) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {
		if (chimeric_alignment->second.filter != NULL)
			continue; // read has already been filtered

		if (chimeric_alignment->second.size() == 3) { // these are alignments of a split read

			// get sequences near breakpoint
			string sequence = "";
			if (chimeric_alignment->second[SPLIT_READ].strand == FORWARD) {
				if (chimeric_alignment->second[SPLIT_READ].preclipping() >= homopolymer_length)
					sequence += chimeric_alignment->second[SPLIT_READ].sequence.substr(chimeric_alignment->second[SPLIT_READ].preclipping() - homopolymer_length, homopolymer_length) + " ";
				if (chimeric_alignment->second[SPLIT_READ].sequence.length() - chimeric_alignment->second[SPLIT_READ].preclipping() >= homopolymer_length)
					sequence += chimeric_alignment->second[SPLIT_READ].sequence.substr(chimeric_alignment->second[SPLIT_READ].preclipping(), homopolymer_length) + " ";
			} else { // chimeric_alignment->second[SPLIT_READ].strand == REVERSE
				if (chimeric_alignment->second[SPLIT_READ].postclipping() >= homopolymer_length)
					sequence += chimeric_alignment->second[SPLIT_READ].sequence.substr(chimeric_alignment->second[SPLIT_READ].sequence.length() - chimeric_alignment->second[SPLIT_READ].postclipping(), homopolymer_length) + " ";
				if (chimeric_alignment->second[SPLIT_READ].sequence.length() - chimeric_alignment->second[SPLIT_READ].postclipping() >= homopolymer_length)
					sequence += chimeric_alignment->second[SPLIT_READ].sequence.substr(chimeric_alignment->second[SPLIT_READ].sequence.length() - chimeric_alignment->second[SPLIT_READ].postclipping() - homopolymer_length, homopolymer_length) + " ";
			}

			// check for homopolymers
			unsigned int run = 1;
			for (unsigned int c = 1; c < sequence.length(); c++) {
				if (sequence[c-1] == sequence[c]) {
					run++;
					if (run == homopolymer_length) {
						if (!is_split_read_spliced(chimeric_alignment->second[SPLIT_READ], exon_annotation_index))
							chimeric_alignment->second.filter = FILTERS.at("homopolymer");
						goto next_read;
					}
				} else {
					run = 1;
				}
			}
	
		}

		++remaining;

		next_read: NULL; // NULL is a dummy statement for the goto label
	}

	return remaining;
}

