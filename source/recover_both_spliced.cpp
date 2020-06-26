#include <map>
#include <tuple>
#include <unordered_map>
#include <vector>
#include "common.hpp"
#include "annotation.hpp"
#include "filter_in_vitro.hpp"
#include "read_stats.hpp"
#include "recover_both_spliced.hpp"

using namespace std;

direction_t opposite_direction(const direction_t direction) {
	return ((direction == DOWNSTREAM) ? UPSTREAM : DOWNSTREAM);
}

unsigned int count_supporting_reads(fusion_t& fusion, unordered_map<gene_t,unsigned int>& read_count_by_gene, const exon_annotation_index_t& exon_annotation_index, const coverage_t& coverage, const unsigned int high_expression_threshold, const int max_exon_size, const unsigned int max_coverage) {

	// check if one of the fusion genes is highly expressed => risk of in vitro-generated artifact
	if (read_count_by_gene[fusion.gene1] > high_expression_threshold || read_count_by_gene[fusion.gene2] > high_expression_threshold) {
		if (fusion.both_breakpoints_spliced() && fusion.discordant_mates <= fusion.split_reads1 + fusion.split_reads2) // only count the most reliable events
			return 1; // ignore the real number of supporting reads and count the event as 1 read
		else
			return 0;
	}

	if (!fusion.both_breakpoints_spliced()) {

		// ignore reads in high-coverage regions; they are likely libprep-mediated artifacts
		unsigned int coverage1 = coverage.get_coverage(fusion.contig1, fusion.breakpoint1, (fusion.direction1 == UPSTREAM) ? DOWNSTREAM : UPSTREAM);
		unsigned int coverage2 = coverage.get_coverage(fusion.contig2, fusion.breakpoint2, (fusion.direction2 == UPSTREAM) ? DOWNSTREAM : UPSTREAM);
		if (coverage1 + coverage2 > fusion.supporting_reads() * max_coverage)
			return 0;

		// count intra-exonic breakpoints only if the exons are not too big
		// otherwise it's probably an libprep-mediated artifact caused by fragments sticking together due to hybridization
		exon_set_t exons;
		get_annotation_by_coordinate(fusion.contig1, fusion.breakpoint1, fusion.breakpoint1, exons, exon_annotation_index);
		for (auto exon = exons.begin(); exon != exons.end(); ++exon)
			if ((**exon).end + 1 - (**exon).start > max_exon_size)
				return 0;
		exons.clear();
		get_annotation_by_coordinate(fusion.contig2, fusion.breakpoint2, fusion.breakpoint2, exons, exon_annotation_index);
		for (auto exon = exons.begin(); exon != exons.end(); ++exon)
			if ((**exon).end + 1 - (**exon).start > max_exon_size)
				return 0;
	}

	// count reads, but ignore multi-mapping reads
	unsigned int multimappers = 0;
	unsigned int unique_mappers = 0;
	for (auto chimeric_alignment = fusion.split_read1_list.begin(); chimeric_alignment != fusion.split_read1_list.end(); ++chimeric_alignment)
		if ((**chimeric_alignment).second.multimapper)
			multimappers++;
		else if ((**chimeric_alignment).second.filter == FILTER_none)
			unique_mappers++;
	for (auto chimeric_alignment = fusion.split_read2_list.begin(); chimeric_alignment != fusion.split_read2_list.end(); ++chimeric_alignment)
		if ((**chimeric_alignment).second.multimapper)
			multimappers++;
		else if ((**chimeric_alignment).second.filter == FILTER_none)
			unique_mappers++;
	for (auto chimeric_alignment = fusion.discordant_mate_list.begin(); chimeric_alignment != fusion.discordant_mate_list.end(); ++chimeric_alignment)
		if ((**chimeric_alignment).second.multimapper)
			multimappers++;
		else if ((**chimeric_alignment).second.filter == FILTER_none)
			unique_mappers++;

	if (multimappers >= 0.5 * (fusion.split_read1_list.size() + fusion.split_read2_list.size() + fusion.discordant_mate_list.size()))
		return 0; // ignore fusions that are supported predominantly by multi-mapping reads, there are too many

	if (unique_mappers == 0)
		return 1; // all supporting reads were discarded by filters => nevertheless count the event as 1 read

	return unique_mappers;
}

unsigned int recover_both_spliced(fusions_t& fusions, const chimeric_alignments_t& chimeric_alignments, const exon_annotation_index_t& exon_annotation_index, const coverage_t& coverage, const unsigned int max_fusions_to_recover, const float high_expression_quantile, const int max_exon_size, const unsigned int max_coverage) {

	// find highly expressed genes; we apply more stringent criteria to them, because they produce lots of artifacts
	unordered_map<gene_t,unsigned int> read_count_by_gene;
	unsigned int high_expression_threshold;
	find_top_expressed_genes(chimeric_alignments, high_expression_quantile, read_count_by_gene, high_expression_threshold);

	// look for any supporting reads between two genes
	map< tuple<gene_t,gene_t,direction_t,direction_t>, vector<fusion_t*> > fusions_by_gene_pair;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter != FILTER_merge_adjacent)
			if (fusion->second.filter == FILTER_none ||
			    fusion->second.filter == FILTER_in_vitro ||
			    fusion->second.filter == FILTER_intronic ||
			    fusion->second.filter == FILTER_relative_support ||
			    fusion->second.filter == FILTER_min_support ||
			    fusion->second.filter == FILTER_inconsistently_clipped && fusion->second.both_breakpoints_spliced())
				if (count_supporting_reads(fusion->second, read_count_by_gene, exon_annotation_index, coverage, high_expression_threshold, max_exon_size, max_coverage) > 0)
					fusions_by_gene_pair[make_tuple(fusion->second.gene1, fusion->second.gene2, (direction_t) fusion->second.direction1, (direction_t) fusion->second.direction2)].push_back(&fusion->second);

	unsigned int remaining = 0;
	const char MODE_COUNTING = 0;
	const char MODE_RECOVER = 1;
	map<unsigned int /*supporting reads*/, unsigned int /*number of recovered fusions*/> recovered_fusions_by_supporting_reads;
	unsigned int min_supporting_reads = 1; // is adjusted dynamically depending on <max_fusions_to_recover>
	for (char mode = MODE_COUNTING; mode <= MODE_RECOVER; ++mode) {
		for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

			if (fusion->second.filter == FILTER_none) { // fusion has not been filtered, no need to recover
				if (mode == MODE_RECOVER)
					remaining++;
				continue;
			}

			if (!fusion->second.both_breakpoints_spliced())
				continue; // only recover fusions with two spliced breakpoints (the purpose of this filter)

			if (fusion->second.gene1 == fusion->second.gene2 || fusion->second.breakpoint_overlaps_both_genes())
				continue; // don't recover intragenic events (this would produce too many hits)

			if (fusion->second.is_read_through())
				continue; // again, there would be too many hits

			if (fusion->second.filter != FILTER_none &&
			    fusion->second.filter != FILTER_relative_support &&
			    fusion->second.filter != FILTER_min_support &&
			    fusion->second.filter != FILTER_in_vitro)
				continue; // we won't recover fusions which were not discarded due to low support

			// count all supporting reads of all fusions between the pair of genes
			unsigned int sum_of_supporting_reads = 0;

			// look for other reads with the same orientation
			auto fusions_of_given_gene_pair = fusions_by_gene_pair.find(make_tuple(fusion->second.gene1, fusion->second.gene2, (direction_t) fusion->second.direction1, (direction_t) fusion->second.direction2));
			if (fusions_of_given_gene_pair != fusions_by_gene_pair.end())
				for (auto another_fusion = fusions_of_given_gene_pair->second.begin(); another_fusion != fusions_of_given_gene_pair->second.end(); ++another_fusion)
					sum_of_supporting_reads += count_supporting_reads(**another_fusion, read_count_by_gene, exon_annotation_index, coverage, high_expression_threshold, max_exon_size, max_coverage);

			// consider reciprocal translocations
			auto reciprocal_fusions_of_given_gene_pair = fusions_by_gene_pair.find(make_tuple(fusion->second.gene1, fusion->second.gene2, opposite_direction(fusion->second.direction1), opposite_direction(fusion->second.direction2)));
			if (reciprocal_fusions_of_given_gene_pair != fusions_by_gene_pair.end())
				for (auto another_fusion = reciprocal_fusions_of_given_gene_pair->second.begin(); another_fusion != reciprocal_fusions_of_given_gene_pair->second.end(); ++another_fusion)
					if (!(**another_fusion).is_read_through())
						// if the fusion is supported by 2 events of which all breakpoints are spliced, we don't question it
						// if the other fusion is not spliced, then we check if the reciprocal fusions support a common genomic breakpoint
						if ((**another_fusion).both_breakpoints_spliced() ||
						    (((fusion->second.direction1 == DOWNSTREAM) != /*xor*/ (fusion->second.breakpoint1 > (**another_fusion).breakpoint1)) &&
						     ((fusion->second.direction2 == DOWNSTREAM) != /*xor*/ (fusion->second.breakpoint2 > (**another_fusion).breakpoint2))))
							sum_of_supporting_reads += count_supporting_reads(**another_fusion, read_count_by_gene, exon_annotation_index, coverage, high_expression_threshold, max_exon_size, max_coverage);

			if (sum_of_supporting_reads >= 2) { // require at least two reads or else the false positive rate sky-rockets
				if (mode == MODE_RECOVER) { // we are in recover mode => actually recover the fusion by clearing the filters
					unsigned int add_one_read_for_proximal_events = (fusion->second.contig1 == fusion->second.contig2 && abs(fusion->second.breakpoint1 - fusion->second.breakpoint2) < 1000000) ? 1 : 0;
					if (fusion->second.supporting_reads() >= min_supporting_reads + add_one_read_for_proximal_events) {
						fusion->second.filter = FILTER_none;
						remaining++;
					}
				} else { // mode == MODE_COUNTING
					// We are in dry-run mode, where we only count how many fusions would be recovered.
					// If too many would be recovered, we increase the minimum required supporting reads.
					recovered_fusions_by_supporting_reads[fusion->second.supporting_reads()]++;
				}
			}
		}

		if (mode == MODE_COUNTING) {
			// The map "recovered_fusions_by_supporting_reads" now holds for each number of supporting reads
			// the number of fusions that would be recovered (e.g., X fusions with Y supporting reads would be recovered).
			// In some samples we observe an extraordinary number of recovered fusions with less than 3 supporting reads.
			// We don't want to recover those, since most of them are probably false positives.
			// Instead, we increase the minimum number of required supporting reads, e.g., if we observe a gazillion
			// of recovered fusions with 2 supporting reads or less, we only recover ones with more than 2, unless
			// we also see a gazillion of those, in which case we only recover ones with more than 3, and so on.
			// => Find the minimum number of supporting reads, where a moderate number of fusions will be recovered.
			unsigned int would_be_recovered = 0;
			for (auto supporting_reads = recovered_fusions_by_supporting_reads.rbegin(); supporting_reads != recovered_fusions_by_supporting_reads.rend(); ++supporting_reads) {
				would_be_recovered += supporting_reads->second;
				if (would_be_recovered >= max_fusions_to_recover) {
					min_supporting_reads = supporting_reads->first + 1;
					break;
				}
			}
		}
	}
	return remaining;
}
