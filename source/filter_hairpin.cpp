#include <cmath>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "filter_hairpin.hpp"

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

bool spliced_distance_too_short(const contig_t contig, const position_t breakpoint1, const position_t breakpoint2, const direction_t direction1, const direction_t direction2, const gene_set_t& common_genes, int min_distance, exon_annotation_index_t& exon_annotation_index) {
	for (gene_set_t::const_iterator common_gene = common_genes.begin(); common_gene != common_genes.end(); ++common_gene)
		if (get_spliced_distance(contig, breakpoint1, breakpoint2, direction1, direction2, *common_gene, exon_annotation_index) <= min_distance)
			return true;
	return false;
}

unsigned int filter_hairpin(chimeric_alignments_t& chimeric_alignments, exon_annotation_index_t& exon_annotation_index) {
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

			position_t breakpoint1 = (i->second[MATE1].strand == FORWARD) ? i->second[MATE1].end : i->second[MATE1].start;
			position_t breakpoint2 = (i->second[MATE2].strand == FORWARD) ? i->second[MATE2].end : i->second[MATE2].start;
			direction_t direction1 = (i->second[MATE1].strand == FORWARD) ? DOWNSTREAM : UPSTREAM;
			direction_t direction2 = (i->second[MATE2].strand == FORWARD) ? DOWNSTREAM : UPSTREAM;

			if (abs(breakpoint1 - breakpoint2) <= 200 ||
			    spliced_distance_too_short(i->second[MATE1].contig, breakpoint1, breakpoint2, direction1, direction2, common_genes, 200, exon_annotation_index) ||
			    is_breakpoint_within_aligned_segment(breakpoint1, i->second[MATE2]) ||
			    is_breakpoint_within_aligned_segment(breakpoint2, i->second[MATE1])) {
				i->second.filters.insert(FILTERS.at("hairpin"));
				continue;
			}

		} else { // split read

			position_t breakpoint_split_read = (i->second[SPLIT_READ].strand == FORWARD) ? i->second[SPLIT_READ].start : i->second[SPLIT_READ].end;
			direction_t direction_split_read = (i->second[SPLIT_READ].strand == FORWARD) ? UPSTREAM : DOWNSTREAM;
			position_t breakpoint_supplementary = (i->second[SUPPLEMENTARY].strand == FORWARD) ? i->second[SUPPLEMENTARY].end : i->second[SUPPLEMENTARY].start;
			direction_t direction_supplementary = (i->second[SUPPLEMENTARY].strand == FORWARD) ? DOWNSTREAM : UPSTREAM;
			position_t gap_start_mate1 = (i->second[MATE1].strand == FORWARD) ? i->second[MATE1].end : i->second[MATE1].start;
			direction_t direction_gap_start_mate1 = (i->second[MATE1].strand == FORWARD) ? DOWNSTREAM : UPSTREAM;
			position_t gap_start_split_read = (i->second[SPLIT_READ].strand == FORWARD) ? i->second[SPLIT_READ].end : i->second[SPLIT_READ].start;
			direction_t direction_gap_start_split_read = (i->second[SPLIT_READ].strand == FORWARD) ? DOWNSTREAM : UPSTREAM;
			position_t start_mate1 = (i->second[MATE1].strand == FORWARD) ? i->second[MATE1].start : i->second[MATE1].end;
			direction_t direction_start_mate1 = (i->second[MATE1].strand == FORWARD) ? UPSTREAM : DOWNSTREAM;
			if (abs(breakpoint_supplementary - breakpoint_split_read) <= 200 ||
			    abs(breakpoint_supplementary - gap_start_mate1) <= 200 ||
			    abs(breakpoint_supplementary - start_mate1) <= 200 ||
			    abs(breakpoint_supplementary - gap_start_split_read) <= 200 ||
			    spliced_distance_too_short(i->second[SPLIT_READ].contig, breakpoint_supplementary, breakpoint_split_read, direction_supplementary, direction_split_read, common_genes, 200, exon_annotation_index) ||
			    spliced_distance_too_short(i->second[MATE1].contig, breakpoint_supplementary, gap_start_mate1, direction_supplementary, direction_gap_start_mate1, common_genes, 200, exon_annotation_index) ||
			    spliced_distance_too_short(i->second[MATE1].contig, breakpoint_supplementary, start_mate1, direction_supplementary, direction_start_mate1, common_genes, 200, exon_annotation_index) ||
			    spliced_distance_too_short(i->second[SPLIT_READ].contig, breakpoint_supplementary, gap_start_split_read, direction_supplementary, direction_gap_start_split_read, common_genes, 200, exon_annotation_index) ||
			    is_breakpoint_within_aligned_segment(breakpoint_split_read, i->second[SUPPLEMENTARY]) ||
			    is_breakpoint_within_aligned_segment(breakpoint_supplementary, i->second[SPLIT_READ]) ||
			    is_breakpoint_within_aligned_segment(breakpoint_supplementary, i->second[MATE1])) {
				i->second.filters.insert(FILTERS.at("hairpin"));
				continue;
			}

		}

		remaining++;
	}
	return remaining;
}

