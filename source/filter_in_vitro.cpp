#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "read_stats.hpp"
#include "filter_in_vitro.hpp"

using namespace std;

// convenience wrapper to lookup a value in an unordered_map or return a default value, if the key does not exist
template <class K, class V> V find_or_default(const unordered_map<K,V>& m, const K& k, const V default_value) {
	auto i = m.find(k);
	return (i == m.end()) ? default_value : i->second;
}

// check if there is a gene which overlaps the given breakpoint and is expressed at a higher level than <highest_expressed_gene>
gene_t find_higher_expressed_gene(const contig_t contig, const position_t breakpoint, const gene_annotation_index_t& gene_annotation_index, const unordered_map<gene_t,unsigned int>& expression_by_gene, gene_t highest_expressed_gene) {
	unsigned int highest_expression = find_or_default(expression_by_gene, highest_expressed_gene, (unsigned int) 0);
	gene_set_t genes_overlapping_breakpoint;
	get_annotation_by_coordinate(contig, breakpoint, breakpoint, genes_overlapping_breakpoint, gene_annotation_index);
	for (gene_set_t::iterator gene = genes_overlapping_breakpoint.begin(); gene != genes_overlapping_breakpoint.end(); ++gene) {
		unsigned int expression = find_or_default(expression_by_gene, *gene, (unsigned int) 0);
		if (expression > highest_expression) {
			highest_expression = expression;
			highest_expressed_gene = *gene;
		}
	}
	return highest_expressed_gene;
}

// make helper class to sort genes by the number of chimeric reads
struct sort_genes_by_reads_t {
	unordered_map<gene_t,unsigned int>* read_count;
	bool operator()(const gene_t x, const gene_t y) const {
		unsigned int reads_x = read_count->at(x);
		unsigned int reads_y = read_count->at(y);
		if (reads_x != reads_y) {
			return reads_x < reads_y;
		} else {
			return x->id < y->id; // ensures deterministic behavior in case of ties
		}
	}
};

void find_top_expressed_genes(const chimeric_alignments_t& chimeric_alignments, const float high_expression_quantile, unordered_map<gene_t,unsigned int>& read_count_by_gene, unsigned int& high_expression_threshold) {

	// RT fusions are mostly observed between highly expressed genes
	// we use the number of chimeric reads as a proxy to measure expression
	for (chimeric_alignments_t::const_iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {
		for (gene_set_t::const_iterator gene = chimeric_alignment->second[MATE1].genes.begin(); gene != chimeric_alignment->second[MATE1].genes.end(); ++gene)
			read_count_by_gene[*gene]++;
		unsigned int mate2 = (chimeric_alignment->second.size() == 2) ? MATE2 : SUPPLEMENTARY;
		for (gene_set_t::const_iterator gene = chimeric_alignment->second[mate2].genes.begin(); gene != chimeric_alignment->second[mate2].genes.end(); ++gene)
			read_count_by_gene[*gene]++;
	}

	// (partially) sort genes by reads using nth_element() to calculate quantiles
	high_expression_threshold = 0;
	if (read_count_by_gene.size() > 0) {

		// put genes in vector for sorting
		vector<gene_t> genes_sorted_by_reads;
		genes_sorted_by_reads.reserve(read_count_by_gene.size());
		for (auto gene = read_count_by_gene.begin(); gene != read_count_by_gene.end(); ++gene)
			genes_sorted_by_reads.push_back(gene->first);

		// calculate position of quantile
		unsigned int quantile = static_cast<int>(floor(high_expression_quantile * genes_sorted_by_reads.size()));
		if (quantile >= genes_sorted_by_reads.size())
			quantile = genes_sorted_by_reads.size() - 1;

		// partially sort using nth_element
		sort_genes_by_reads_t sort_genes_by_reads;
		sort_genes_by_reads.read_count = &read_count_by_gene;
		nth_element(genes_sorted_by_reads.begin(), genes_sorted_by_reads.begin()+quantile, genes_sorted_by_reads.end(), sort_genes_by_reads);

		// extract quantile
		high_expression_threshold = read_count_by_gene[genes_sorted_by_reads[quantile]];
	}
}

unsigned int filter_in_vitro(fusions_t& fusions, const chimeric_alignments_t& chimeric_alignments, const float high_expression_quantile, const gene_annotation_index_t& gene_annotation_index, const coverage_t& coverage) {

	// older version of STAR occasionally clipped discordant mates for no good reason,
	// which appeared as though the mate overlaps a breakpoint
	// => only consider discordant mates with this many clipped bases (or more) to overlap the breakpoint
	const unsigned int min_clipped_length = 3;
	// genes fused during reverse transcription (RT)
	// often have multiple breakpoints within exons (rather than at splice-sites)
	// => we consider the following to be many exonic breakpoints
	const unsigned int max_exonic_breakpoints_by_gene_pair = 8;

	// count the number of breakpoints within exons for each gene pair
	unordered_map< tuple<gene_t/*gene1*/,gene_t/*gene2*/>, unsigned int > exonic_breakpoints_by_gene_pair;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.gene1 != fusion->second.gene2 && // it is perfectly normal to have many breakpoints within the same gene (hairpin fusions)
		    !fusion->second.spliced1 && !fusion->second.spliced2 && // breakpoints at splice sites are almost exclusively a result of splicing and thus, no RT-mediated fusions
		    fusion->second.exonic1 && fusion->second.exonic2 && // RT fusions only contain spliced transcripts, so we ignore intronic/intergenic breakpoints
		    fusion->second.split_read1_list.size() + fusion->second.split_read2_list.size() > 0 && // require a split read for exact location of the breakpoint
		    fusion->second.filter != FILTER_merge_adjacent && // slightly varying alignments may lead to adjacent breakpoints, we should not count them as separate breakpoints
		    fusion->second.filter != FILTER_uninteresting_contigs) { // skip uninteresting contigs to save some runtime/memory
			exonic_breakpoints_by_gene_pair[make_tuple(fusion->second.gene1, fusion->second.gene2)]++;
			exonic_breakpoints_by_gene_pair[make_tuple(fusion->second.gene2, fusion->second.gene1)]++;
		}
	}

	// RT fusions are mostly observed between highly expressed genes
	unordered_map<gene_t,unsigned int> read_count_by_gene;
	unsigned int high_expression_threshold;
	find_top_expressed_genes(chimeric_alignments, high_expression_quantile, read_count_by_gene, high_expression_threshold);

	// remove fusions with a lot of characteristics of RT-mediated fusions
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != FILTER_none && // fusion has already been filtered
		    // also tag filtered events, if they are spliced to prevent the filters 'spliced' and 'many_spliced' from recovering them
		    !((fusion->second.spliced1 || fusion->second.spliced2) && (fusion->second.filter == FILTER_relative_support || fusion->second.filter == FILTER_min_support || fusion->second.filter == FILTER_homopolymer)))
			continue;

		// RT-mediated fusions often have both breakpoints within exons,
		// but sometimes one (or even both) of them can be in introns or even at splice-sites
		// the latter only happens with very highly expressed genes
		// => count the number of potentially RT-mediated breakpoints and
		//    penalize events with one/two RT-mediated breakpoints more
		float potential_rt_breakpoints = 0;
		if (!fusion->second.exonic1) {
			potential_rt_breakpoints += 0.5; // intron => unclear if it is a RT artifact (not likely, but possible)
		} else if (!fusion->second.spliced1) {
			potential_rt_breakpoints += 1; // breakpoint inside exon => very likely RT artifact
		} // else: splice-site => probably from splicing and not RT-mediated
		if (!fusion->second.exonic2) {
			potential_rt_breakpoints += 0.5; // intron => unclear if it is a RT artifact (not likely, but possible)
		} else if (!fusion->second.spliced2) {
			potential_rt_breakpoints += 1; // breakpoint inside exon => very likely RT artifact
		} // else: splice-site => probably from splicing and not RT-mediated

		// RT artifacts usually have no/few split reads
		// then again, some true events have no split reads, because of alignment issues (e.g., many inserted non-template bases)
		// we check if some of the discordant mates overlap the breakpoint and are clipped there
		// if so, they are counted as split reads
		unsigned int clipped_discordant_mates1 = 0;
		unsigned int clipped_discordant_mates2 = 0;
		for (auto discordant_mates = fusion->second.discordant_mate_list.begin(); discordant_mates != fusion->second.discordant_mate_list.end(); ++discordant_mates) {
			if ((**discordant_mates).second.filter == FILTER_none) {
				for (mates_t::iterator mate = (**discordant_mates).second.begin(); mate != (**discordant_mates).second.end(); ++mate) {
					if (mate->strand == FORWARD && mate->postclipping() >= min_clipped_length) {
						if (mate->contig == fusion->second.contig1 && mate->end == fusion->second.breakpoint1) {
							clipped_discordant_mates1++;
						} else if (mate->contig == fusion->second.contig2 && mate->end == fusion->second.breakpoint2) {
							clipped_discordant_mates2++;
						}
					} else if (mate->strand == REVERSE && mate->preclipping() >= min_clipped_length) {
						if (mate->contig == fusion->second.contig1 && mate->start == fusion->second.breakpoint1) {
							clipped_discordant_mates1++;
						} else if (mate->contig == fusion->second.contig2 && mate->start == fusion->second.breakpoint2) {
							clipped_discordant_mates2++;
						}
					}
				}
			}
		}

		// add the clipped discordant mates to the split reads
		// we take the min of both counts, because RT artifacts often have clipped reads on only one side
		// and because we would potentially count the same fragment twice, if both discordant mates overlap both breakpoints
		unsigned int total_split_reads = min(clipped_discordant_mates1, clipped_discordant_mates2) + fusion->second.split_reads1 + fusion->second.split_reads2;

		// find the highest expressed genes overlapping the breakpoint
		gene_t gene1 = find_higher_expressed_gene(fusion->second.contig1, fusion->second.breakpoint1, gene_annotation_index, read_count_by_gene, fusion->second.gene1);
		gene_t gene2 = find_higher_expressed_gene(fusion->second.contig2, fusion->second.breakpoint2, gene_annotation_index, read_count_by_gene, fusion->second.gene2);

		// use the number of chimeric reads as a proxy to measure expression
		unsigned int gene1_expression = find_or_default(read_count_by_gene, gene1, (unsigned int) 0);
		unsigned int gene2_expression = find_or_default(read_count_by_gene, gene2, (unsigned int) 0);

		// find out how many non-spliced, non-intronic breakpoints there are between the fused genes (or genes overlapping the breakpoints)
		unsigned int exonic_breakpoints = max(
			find_or_default(exonic_breakpoints_by_gene_pair, make_tuple(gene1, gene2), (unsigned int) 0),
			find_or_default(exonic_breakpoints_by_gene_pair, make_tuple(fusion->second.gene1, fusion->second.gene2), (unsigned int) 0)
		);

		int coverage1 = coverage.get_coverage(fusion->second.contig1, fusion->second.breakpoint1, (fusion->second.direction1 == UPSTREAM) ? DOWNSTREAM : UPSTREAM);
		int coverage2 = coverage.get_coverage(fusion->second.contig2, fusion->second.breakpoint2, (fusion->second.direction2 == UPSTREAM) ? DOWNSTREAM : UPSTREAM);

		// check if the event exhibits the characteristics of a RT-mediated fusion
		if (
			// RT artifacts often have very few split reads
			// => check if the event has <=2 split reads or fewer than 1 supporting read per 10,000 chimeric reads
			total_split_reads <= 2 + 0.0001 * (gene1_expression + gene2_expression) &&
			// RT artifacts often have many discordant mates, but only few split reads
			// => the number of split reads and discordant mates must be unbalanced (unless the total number of reads is very low)
			(total_split_reads * 2 <= fusion->second.discordant_mates || total_split_reads <= 2) &&
			// RT artifacts occur between highly expressed genes
			// => the sum of the chimeric reads of the fused genes must be in the top percentile
			gene1_expression + gene2_expression > high_expression_threshold &&
			// make an exception for fusions that make up a substantial fraction of the overall expression of a gene and have a spliced breakpoint
			!(fusion->second.supporting_reads() >= 10 && ((int) fusion->second.supporting_reads() * 4) >= max(coverage1, coverage2) &&
			  coverage1 > (int) fusion->second.supporting_reads() && coverage2 > (int) fusion->second.supporting_reads() &&
			  (fusion->second.spliced1 || fusion->second.spliced2) && ((fusion->second.spliced1 || !fusion->second.exonic1) && (fusion->second.spliced2 || !fusion->second.exonic2))) &&
			// RT artifacts typically have breakpoints inside exons, otherwise they might be true events
			// => if we see an intronic/spliced breakpoint, only discard the event, if the genes are expressed at extreme levels
			(
				// both breakpoints in exons or one in exon and the other in intron => discard event
				potential_rt_breakpoints > 1 ||
				// both breakpoints in introns or one breakpoint at a splice-site => only discard event, if one of the genes is in the top quantile expression
				potential_rt_breakpoints > 0 && (gene1_expression > high_expression_threshold || gene2_expression > high_expression_threshold) ||
				// both breakpoints at splice-sites => only discard even, if both genes are in the top quantile expression
				gene1_expression > 2 * high_expression_threshold || gene2_expression > 2 * high_expression_threshold || gene1_expression > high_expression_threshold && gene2_expression > high_expression_threshold ||
				// there are many exonic breakpoints between the fused genes => probably the genes have sequence homology, which often yields spliced events
				exonic_breakpoints > max_exonic_breakpoints_by_gene_pair ||
				// the fusion is spliced but only supported by very few reads => ignore that the event is spliced, because it is probably trans-spliced
				fusion->second.supporting_reads() <= 1
			)
		   )
			fusion->second.filter = FILTER_in_vitro;

	}

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter == FILTER_none)
			++remaining;
	return remaining;
}

