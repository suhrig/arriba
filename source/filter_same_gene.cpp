#include <cmath>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "filter_same_gene.hpp"

using namespace std;

bool is_breakpoint_within_aligned_segment(const position_t breakpoint, const alignment_t& alignment) {
	// find out where the read aligns and check if the breakpoint is within this aligned section
	position_t reference_position = alignment.start;
	for (unsigned int cigar_element = 0; cigar_element < alignment.cigar.size(); ++cigar_element) {
		switch (alignment.cigar.operation(cigar_element)) {
			case BAM_CREF_SKIP:
			case BAM_CDEL:
				reference_position += alignment.cigar.op_length(cigar_element);
				break;
			case BAM_CMATCH:
				if (breakpoint >= reference_position && breakpoint <= reference_position + alignment.cigar.op_length(cigar_element))
					return	true;
				reference_position += alignment.cigar.op_length(cigar_element);
				break;
		}
	}
	return false;
}

unsigned int filter_same_gene(chimeric_alignments_t& chimeric_alignments) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {

		if (!i->second.filters.empty())
			continue; // the read has already been filtered

		// check if mate1 and mate2 map to the same gene(s)
		// if so, remove them
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

			position_t breakpoint1 = (i->second[MATE1].strand == FORWARD) ? i->second[MATE1].end : i->second[MATE1].start;
			position_t breakpoint2 = (i->second[MATE2].strand == FORWARD) ? i->second[MATE2].end : i->second[MATE2].start;

			if (abs(breakpoint1 - breakpoint2) <= 200 ||
			    is_breakpoint_within_aligned_segment(breakpoint1, i->second[MATE2]) ||
			    is_breakpoint_within_aligned_segment(breakpoint2, i->second[MATE1])) {
				i->second.filters.insert(FILTERS.at("same_gene"));
				continue;
			}

		} else { // split read

			if (i->second[SPLIT_READ].strand == FORWARD && i->second[SUPPLEMENTARY].strand == FORWARD && i->second[SPLIT_READ].start >= i->second[SUPPLEMENTARY].end ||
			    i->second[SPLIT_READ].strand == REVERSE && i->second[SUPPLEMENTARY].strand == REVERSE && i->second[SPLIT_READ].end   <= i->second[SUPPLEMENTARY].start) {
				i->second.filters.insert(FILTERS.at("same_gene")); // normal alignment
				continue;
			}

			position_t breakpoint_split_read = (i->second[SPLIT_READ].strand == FORWARD) ? i->second[SPLIT_READ].start : i->second[SPLIT_READ].end;
			position_t breakpoint_supplementary = (i->second[SUPPLEMENTARY].strand == FORWARD) ? i->second[SUPPLEMENTARY].end : i->second[SUPPLEMENTARY].start;
			position_t gap_start_mate1 = (i->second[MATE1].strand == FORWARD) ? i->second[MATE1].end : i->second[MATE1].start;
			position_t gap_start_split_read = (i->second[SPLIT_READ].strand == FORWARD) ? i->second[SPLIT_READ].end : i->second[SPLIT_READ].start;
			if (abs(breakpoint_supplementary - gap_start_mate1) <= 200 ||
			    abs(breakpoint_supplementary - gap_start_split_read) <= 200 ||
			    is_breakpoint_within_aligned_segment(breakpoint_split_read, i->second[SUPPLEMENTARY]) ||
			    is_breakpoint_within_aligned_segment(breakpoint_supplementary, i->second[SPLIT_READ]) ||
			    is_breakpoint_within_aligned_segment(breakpoint_supplementary, i->second[MATE1])) {
				i->second.filters.insert(FILTERS.at("same_gene"));
				continue;
			}

		}

		remaining++;
	}
	return remaining;
}

