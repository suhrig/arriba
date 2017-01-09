#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include "common.hpp"
#include "annotation.hpp"
#include "filter_pcr_fusions.hpp"

using namespace std;















unsigned int filter_pcr_fusions(fusions_t& fusions, const float max_pcr_fusion_score, const unsigned int max_exonic_breakpoints, const unsigned int max_partners_with_many_exonic_breakpoints, const unsigned int min_split_reads, const unsigned int max_exonic_breakpoints_for_gene_pair) {

	unordered_map< gene_t,unsigned int > exonic_breakpoint_count; // count the number of fusions with exonic (non-spliced) breakpoints for each gene
	unordered_map< gene_t,unsigned int > partners_with_many_exonic_breakpoints; // count the number of gene partners which have many exonic breakpoints for each gene
	unordered_map< tuple<gene_t/*gene1*/, gene_t/*gene2*/>, unsigned int > exonic_breakpoints_by_gene_pair; // count the number of exonic breakpoints for each gene pair
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.gene1 != fusion->second.gene2 && // there are often many breakpoints within the same gene (probably hairpin fusions)
		    !fusion->second.spliced1 && !fusion->second.spliced2 && // breakpoints at splice sites are almost exclusively a result of splicing and thus, no PCR-mediated fusions
		    (!fusion->second.overlap_duplicate1 || !fusion->second.overlap_duplicate2) && // don't count breakpoints twice, if they are in overlapping genes
		    fusion->second.exonic1 && fusion->second.exonic2 && // PCR fusions only contain spliced transcripts, so we ignore intronic/intergenic breakpoints
		    fusion->second.filter != FILTERS.at("uninteresting_contigs") && // ignore fusions with mitochondrial DNA and decoy sequences, there are usually tons of those
		    fusion->second.filter != FILTERS.at("merge_adjacent")) { // slightly varying alignments may lead to adjacent breakpoints, we should not count them as separate breakpoints

			if (fusion->second.split_read1_list.size() + fusion->second.split_read2_list.size()) { // this also counts discarded reads
				if (!fusion->second.overlap_duplicate2)
					exonic_breakpoint_count[fusion->second.gene1]++;
				if (!fusion->second.overlap_duplicate1)
					exonic_breakpoint_count[fusion->second.gene2]++;
				if (++exonic_breakpoints_by_gene_pair[make_tuple(fusion->second.gene1, fusion->second.gene2)] == max_exonic_breakpoints) {
					partners_with_many_exonic_breakpoints[fusion->second.gene1]++;
					partners_with_many_exonic_breakpoints[fusion->second.gene2]++;
				}
			}
		}
	}

	// calculate score which reflects the likelihood of PCR fusions for each gene
	unordered_map< gene_t,float > pcr_fusion_scores;
	for (auto partners = partners_with_many_exonic_breakpoints.begin(); partners != partners_with_many_exonic_breakpoints.end(); ++partners)
		pcr_fusion_scores[partners->first] = partners->second * log10(exonic_breakpoint_count[partners->first]); // we take the logarithm so that the factors have about equal order of magnitude / weight

	// remove fusions with a lot of evidence for being PCR-mediated
	map< pair<gene_t,gene_t>, bool > filtered_gene_pairs;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != NULL)
			continue; // fusion has already been filtered

		//TODO how about we don't care if it is spliced in the PCR amplified gene and only care about if it is spliced in the other gene?
		if ((fusion->second.split_reads1 + fusion->second.split_reads2 < fusion->second.discordant_mates || fusion->second.split_reads1 + fusion->second.split_reads2 <= min_split_reads) &&
		    (!fusion->second.spliced1 && !fusion->second.spliced2 || fusion->second.split_reads1 + fusion->second.split_reads2 <= min_split_reads || pcr_fusion_scores[fusion->second.gene1] >= max_pcr_fusion_score && pcr_fusion_scores[fusion->second.gene2] >= max_pcr_fusion_score) && // discard fusions, which are not spliced / have few reads / are between genes with high PCR fusion scores
		    (fusion->second.exonic1 || fusion->second.exonic2) && // one of the breakpoints must be exonic (in theory, both should be, but there are too many exceptions)
		    (pcr_fusion_scores[fusion->second.gene1] + pcr_fusion_scores[fusion->second.gene2] >= max_pcr_fusion_score && // PCR fusion score must be above threshold
		     (partners_with_many_exonic_breakpoints[fusion->second.gene1] >= max_partners_with_many_exonic_breakpoints || partners_with_many_exonic_breakpoints[fusion->second.gene2] >= max_partners_with_many_exonic_breakpoints) || // one of the partners must have many other partners with many exonic breakpoints
		      exonic_breakpoints_by_gene_pair[make_tuple(fusion->second.gene1, fusion->second.gene2)] > max_exonic_breakpoints_for_gene_pair)) { // the gene pair has many breakpoints between each other
			filtered_gene_pairs[make_pair(fusion->second.gene1, fusion->second.gene2)];
			fusion->second.filter = FILTERS.at("pcr_fusions");
		}
	}

	// filter fusions which overlap with genes removed by this filter
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter != NULL)
			continue; // fusion has already been filtered

		if ((fusion->second.split_reads1 + fusion->second.split_reads2 < fusion->second.discordant_mates || fusion->second.split_reads1 + fusion->second.split_reads2 <= min_split_reads) &&
		    (!fusion->second.spliced1 && !fusion->second.spliced2 || fusion->second.split_reads1 + fusion->second.split_reads2 <= min_split_reads) && // discard fusions, which are not spliced / have few reads
		    (fusion->second.exonic1 || fusion->second.exonic2)) { // one of the breakpoints must be exonic (in theory, both should be, but there are too many exceptions)

			bool breakpoints_overlap_filtered_gene_pair = false;
			for (auto filtered_gene_pair = filtered_gene_pairs.begin(); filtered_gene_pair != filtered_gene_pairs.end(); ++filtered_gene_pair) {
				if (fusion->second.contig1 == filtered_gene_pair->first.first->contig && fusion->second.breakpoint1 >= filtered_gene_pair->first.first->start && fusion->second.breakpoint1 <= filtered_gene_pair->first.first->end &&
				    fusion->second.contig2 == filtered_gene_pair->first.second->contig && fusion->second.breakpoint2 >= filtered_gene_pair->first.second->start && fusion->second.breakpoint2 <= filtered_gene_pair->first.second->end) {
					breakpoints_overlap_filtered_gene_pair = true;
					break;
				}
			}

			if (breakpoints_overlap_filtered_gene_pair)
				fusion->second.filter = FILTERS.at("pcr_fusions");
		}

	}

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter == NULL)
			++remaining;
	return remaining;
}

