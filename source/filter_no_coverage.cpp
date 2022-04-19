#include "common.hpp"
#include "annotation.hpp"
#include "filter_no_coverage.hpp"
#include "read_stats.hpp"

using namespace std;

unsigned int filter_no_coverage(fusions_t& fusions, const coverage_t& coverage, const exon_annotation_index_t& exon_annotation_index) {

	const int scan_range = 200; // look for coverage in this range around the breakpoint

	// for each fusion, check if there is any coverage around the breakpoint
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != FILTER_none)
			continue; // fusion has already been filtered

		if (!fusion->second.is_read_through()) {
			if (fusion->second.split_reads1 + fusion->second.split_reads2 != 0 &&
			    fusion->second.split_reads1 + fusion->second.discordant_mates != 0 &&
			    fusion->second.split_reads2 + fusion->second.discordant_mates != 0) {
				++remaining;
				continue; // don't filter fusions with high support
			}
			if (fusion->second.spliced1 || fusion->second.spliced2) {
				++remaining;
				continue; // don't filter spliced breakpoints (they are more credible)
			}
		} else { // read-through
			// most read-through fusions have at least one breakpoint at a splice-site
			// => require both breakpoints to be at splice-sites to be a little more strict
			if (fusion->second.spliced1 && fusion->second.spliced2) {
				++remaining;
				continue;
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
					start -= scan_range;
				end = max(fusion->second.breakpoint1 + scan_range, fusion->second.anchor_start1);
			} else {
				start = min(fusion->second.breakpoint1 - scan_range, fusion->second.anchor_start1);
				end = fusion->second.breakpoint1;
				if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0)
					end += scan_range;
			}
			if (fusion->second.direction1 == UPSTREAM && !coverage.fragment_starts_here(fusion->second.contig1, start, end) ||
			    fusion->second.direction1 == DOWNSTREAM && !coverage.fragment_ends_here(fusion->second.contig1, start, end)) {
				fusion->second.filter = FILTER_no_coverage;
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
					start -= scan_range;
				end = max(fusion->second.breakpoint2 + scan_range, fusion->second.anchor_start2);
			} else {
				start = min(fusion->second.breakpoint2 - scan_range, fusion->second.anchor_start2);
				end = fusion->second.breakpoint2;
				if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0)
					end += scan_range;
			}
			if (fusion->second.direction2 == UPSTREAM && !coverage.fragment_starts_here(fusion->second.contig2, start, end) ||
			    fusion->second.direction2 == DOWNSTREAM && !coverage.fragment_ends_here(fusion->second.contig2, start, end)) {
				fusion->second.filter = FILTER_no_coverage;
				continue;
			}
		}

		remaining++;
	}

	return remaining;
}

