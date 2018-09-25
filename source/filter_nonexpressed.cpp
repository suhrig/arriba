#include "common.hpp"
#include "annotation.hpp"
#include "filter_nonexpressed.hpp"
#include "read_stats.hpp"

using namespace std;

unsigned int filter_nonexpressed(fusions_t& fusions, const coverage_t& coverage, const exon_annotation_index_t& exon_annotation_index, const int max_mate_gap) {

	// for each fusion, check if there is any expression around the breakpoint
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != NULL)
			continue; // fusion has already been filtered

		if (!fusion->second.is_read_through()) {

			if (fusion->second.split_reads1 + fusion->second.split_reads2 != 0 &&
			    fusion->second.split_reads1 + fusion->second.discordant_mates != 0 &&
			    fusion->second.split_reads2 + fusion->second.discordant_mates != 0) {
				++remaining;
				continue; // don't filter fusions with high support
			}

			if (fusion->second.exonic1 && fusion->second.exonic2 ||
			    fusion->second.spliced1 || fusion->second.spliced2) {
				++remaining;
				continue; // don't filter spliced/exonic breakpoints (they probably would not be discarded anyway, because there is plenty of coverage in exons)
			}

		}  else { // read-through

			if (fusion->second.spliced1 && fusion->second.spliced2) {
				++remaining;
				continue; // don't filter spliced breakpoints (they are more credible)
			}
		}

		position_t start, end;
		bool is_in_terminal_exon;

		// check if breakpoint1 is in a terminal exon
		exon_set_t exons;
		get_annotation_by_coordinate(fusion->second.contig1, fusion->second.breakpoint1, fusion->second.breakpoint1, exons, exon_annotation_index);
		is_in_terminal_exon = false;
		for (auto exon = exons.begin(); exon != exons.end() && !is_in_terminal_exon; ++exon)
			if ((**exon).gene == fusion->second.gene1 && ((**exon).previous_exon == NULL || (**exon).next_exon == NULL))
				is_in_terminal_exon = true;

		if (!is_in_terminal_exon) {
			// check if there is coverage around breakpoint1
			if (fusion->second.direction1 == UPSTREAM) {
				start = fusion->second.breakpoint1;
				if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0)
					start -= max_mate_gap;
				end = max(fusion->second.breakpoint1 + max_mate_gap, fusion->second.anchor_start1);
			} else {
				start = min(fusion->second.breakpoint1 - max_mate_gap, fusion->second.anchor_start1);
				end = fusion->second.breakpoint1;
				if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0)
					end += max_mate_gap;
			}
			if (fusion->second.direction1 == UPSTREAM && !coverage.fragment_starts_here(fusion->second.contig1, start, end) ||
			    fusion->second.direction1 == DOWNSTREAM && !coverage.fragment_ends_here(fusion->second.contig1, start, end)) {
				fusion->second.filter = FILTERS.at("non_expressed");
				continue;
			}
		}

		// check if breakpoint2 is in a terminal exon
		exons.clear();
		get_annotation_by_coordinate(fusion->second.contig2, fusion->second.breakpoint2, fusion->second.breakpoint2, exons, exon_annotation_index);
		is_in_terminal_exon = false;
		for (auto exon = exons.begin(); exon != exons.end() && !is_in_terminal_exon; ++exon)
			if ((**exon).gene == fusion->second.gene2 && ((**exon).previous_exon == NULL || (**exon).next_exon == NULL))
				is_in_terminal_exon = true;

		if (!is_in_terminal_exon) {
			// check if there is coverage around breakpoint2
			if (fusion->second.direction2 == UPSTREAM) {
				start = fusion->second.breakpoint2;
				if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0)
					start -= max_mate_gap;
				end = max(fusion->second.breakpoint2 + max_mate_gap, fusion->second.anchor_start2);
			} else {
				start = min(fusion->second.breakpoint2 - max_mate_gap, fusion->second.anchor_start2);
				end = fusion->second.breakpoint2;
				if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0)
					end += max_mate_gap;
			}
			if (fusion->second.direction2 == UPSTREAM && !coverage.fragment_starts_here(fusion->second.contig2, start, end) ||
			    fusion->second.direction2 == DOWNSTREAM && !coverage.fragment_ends_here(fusion->second.contig2, start, end)) {
				fusion->second.filter = FILTERS.at("non_expressed");
				continue;
			}
		}

		remaining++;
	}

	return remaining;
}

