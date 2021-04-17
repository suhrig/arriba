#include "common.hpp"
#include "annotation.hpp"
#include "read_stats.hpp"
#include "recover_internal_tandem_duplication.hpp"

using namespace std;

unsigned int recover_internal_tandem_duplication(fusions_t& fusions, const chimeric_alignments_t& chimeric_alignments, const coverage_t& coverage, const exon_annotation_index_t& exon_annotation_index, const unsigned int max_itd_length) {

	unsigned int min_supporting_reads = 10; // only recover fusions with this many supporting reads
	float min_fraction_of_coverage = 0.05; // the fusion-supporting reads must make up a substantial fraction of the coverage

	// estimate duplication rate, because the calculated coverage includes duplicates
	unsigned int duplicates = 0;
	for (chimeric_alignments_t::const_iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment)
		if (chimeric_alignment->second.filter == FILTER_duplicates)
			duplicates++;
	float duplication_rate = 1.0 * duplicates / chimeric_alignments.size();

	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		// only recover candidates that were discarded for a select number of reasons
		if (fusion->second.filter != FILTER_relative_support &&
		    fusion->second.filter != FILTER_intragenic_exonic &&
		    fusion->second.filter != FILTER_hairpin &&
		    fusion->second.filter != FILTER_inconsistently_clipped &&
		    fusion->second.filter != FILTER_mismatches)
			continue;

		// is it an internal tandem duplication?
		if (fusion->second.gene1 == fusion->second.gene2 &&
		    fusion->second.exonic1 && fusion->second.exonic2 &&
		    fusion->second.direction1 == UPSTREAM && fusion->second.direction2 == DOWNSTREAM &&
		    fusion->second.gene1->is_protein_coding &&
		    ((unsigned int) fusion->second.breakpoint2 - fusion->second.breakpoint1) < max_itd_length) {

			// both breakpoints must be in the same exon and in a coding region
			exon_set_t exons;
			get_annotation_by_coordinate(fusion->second.contig1, fusion->second.breakpoint1, fusion->second.breakpoint2, exons, exon_annotation_index);
			bool is_in_coding_region = false;
			for (auto exon = exons.begin(); exon != exons.end(); ++exon)
				if ((**exon).gene == fusion->second.gene1 &&
				    (**exon).coding_region_start <= fusion->second.breakpoint1 && (**exon).coding_region_end >= fusion->second.breakpoint1 &&
				    (**exon).coding_region_start <= fusion->second.breakpoint2 && (**exon).coding_region_end >= fusion->second.breakpoint2)
					is_in_coding_region = true;
			if (!is_in_coding_region)
				continue;

			// the read support must be very convincing
			int coverage1 = coverage.get_coverage(fusion->second.contig1, fusion->second.breakpoint1, (fusion->second.direction1 == UPSTREAM) ? DOWNSTREAM : UPSTREAM);
			int coverage2 = coverage.get_coverage(fusion->second.contig2, fusion->second.breakpoint2, (fusion->second.direction2 == UPSTREAM) ? DOWNSTREAM : UPSTREAM);
			unsigned int split_reads = 0;
			for (auto split_read = fusion->second.split_read1_list.begin(); split_read != fusion->second.split_read1_list.end(); ++split_read)
				if ((**split_read).second.filter == FILTER_none || (**split_read).second.filter == FILTER_hairpin || (**split_read).second.filter == FILTER_inconsistently_clipped || (**split_read).second.filter == FILTER_mismatches)
					split_reads++;
			for (auto split_read = fusion->second.split_read2_list.begin(); split_read != fusion->second.split_read2_list.end(); ++split_read)
				if ((**split_read).second.filter == FILTER_none || (**split_read).second.filter == FILTER_hairpin || (**split_read).second.filter == FILTER_inconsistently_clipped || (**split_read).second.filter == FILTER_mismatches)
					split_reads++;

			if (split_reads > min_supporting_reads && 1.0 * split_reads/max(coverage1,coverage2)/(1-duplication_rate) > min_fraction_of_coverage) {
				// clear filters of fusion candidate and of supporting reads
				fusion->second.filter = FILTER_none;
				for (auto split_read = fusion->second.split_read1_list.begin(); split_read != fusion->second.split_read1_list.end(); ++split_read)
					if ((**split_read).second.filter == FILTER_hairpin || (**split_read).second.filter == FILTER_inconsistently_clipped || (**split_read).second.filter == FILTER_mismatches) {
						(**split_read).second.filter = FILTER_none;
						fusion->second.split_reads1++;
					}
				for (auto split_read = fusion->second.split_read2_list.begin(); split_read != fusion->second.split_read2_list.end(); ++split_read)
					if ((**split_read).second.filter == FILTER_hairpin || (**split_read).second.filter == FILTER_inconsistently_clipped || (**split_read).second.filter == FILTER_mismatches) {
						(**split_read).second.filter = FILTER_none;
						fusion->second.split_reads2++;
					}
			}

		}
	}

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter == FILTER_none)
			remaining++;
	return remaining;
}

