#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <tuple>
#include <string>
#include <set>
#include <unordered_map>
#include <vector>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "filter_relative_support.hpp"

using namespace std;

void estimate_expected_fusions(fusions_t& fusions, const unsigned long int mapped_reads, const exon_annotation_index_t& exon_annotation_index) {

	// find all fusion partners for each gene
	unordered_map< gene_t,gene_set_t > fusion_partners;
	unordered_map< tuple<gene_t,position_t,position_t>,char > overlap_duplicates;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter == FILTER_none && fusion->second.gene1 != fusion->second.gene2) {
			if (!overlap_duplicates[make_tuple(fusion->second.gene2, fusion->second.breakpoint1, fusion->second.breakpoint2)]++)
				fusion_partners[fusion->second.gene2].insert(fusion->second.gene1);
			if (!overlap_duplicates[make_tuple(fusion->second.gene1, fusion->second.breakpoint1, fusion->second.breakpoint2)]++)
				fusion_partners[fusion->second.gene1].insert(fusion->second.gene2);
		}
	}
	overlap_duplicates.clear();

	// count the number of fusion partners for each gene
	// fusions with genes that have more fusion partners are ignored
	unordered_map<gene_t,int> fusion_partner_count;
	for (auto fusion_partner1 = fusion_partners.begin(); fusion_partner1 != fusion_partners.end(); ++fusion_partner1) {
		for (auto fusion_partner2 = fusion_partner1->second.begin(); fusion_partner2 != fusion_partner1->second.end(); ++fusion_partner2) {
			if (fusion_partner1->second.size() >= fusion_partners[*fusion_partner2].size()) {
				fusion_partner_count[fusion_partner1->first]++;
			}
		}
	}

	// estimate the fraction of breakpoints by location (splice-site vs. exon vs. intron)
	// non-spliced breakpoints get a penalty based on how much more frequent they are than spliced breakpoints
	unsigned int spliced_breakpoints = 0;
	unsigned int exonic_breakpoints = 0;
	unsigned int intronic_breakpoints = 0;
	unsigned int exonic_intronic_breakpoints = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter == FILTER_none &&
		    (fusion->second.contig1 != fusion->second.contig2 || fusion->second.breakpoint2 - fusion->second.breakpoint1 > 500000) && // ignore proximity artifacts
		    fusion->second.supporting_reads() >= 2 && fusion->second.split_reads1 + fusion->second.split_reads2 > 0 && // require at least 2 reads, because most events with 1 read are artifacts
		    !fusion->second.gene1->is_dummy && !fusion->second.gene2->is_dummy) {
			if (fusion->second.spliced1 || fusion->second.spliced2)
				spliced_breakpoints++;
			else if (fusion->second.exonic1 && fusion->second.exonic2)
				exonic_breakpoints++;
			else if (!fusion->second.exonic1 && !fusion->second.exonic2)
				intronic_breakpoints++;
			else
				exonic_intronic_breakpoints++;
		}
	}
	// use some reasonable default values if there are not enough data points to estimate the fractions accurately
	if (spliced_breakpoints + exonic_breakpoints + intronic_breakpoints + exonic_intronic_breakpoints < 100 ||
	    spliced_breakpoints == 0 || exonic_breakpoints == 0 || intronic_breakpoints == 0 || exonic_intronic_breakpoints == 0) {
		spliced_breakpoints = 10;
		exonic_breakpoints = 65;
		intronic_breakpoints = 10;
		exonic_intronic_breakpoints = 15;
	}

	// penalize intragenic events according to event type (inversion vs. duplication),
	// because some libraries produce a huge amount of artifacts of on of these two types of events:
	// stranded libraries produce many inversions/unstranded libraries produce many duplications
	unsigned int intragenic_duplications = 0;
	unsigned int intragenic_inversions = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter == FILTER_none && fusion->second.gene1 == fusion->second.gene2 && fusion->second.split_reads1 + fusion->second.split_reads2 >= 2) {
			if (fusion->second.direction1 == UPSTREAM && fusion->second.direction2 == DOWNSTREAM)
				intragenic_duplications++;
			else if (fusion->second.direction1 == fusion->second.direction2)
				intragenic_inversions++;
		}
	}
	// use reasonable defaut values, if sample size is too small
	if (intragenic_inversions + intragenic_duplications < 100) {
		intragenic_inversions = 1;
		intragenic_duplications = 1;
	}

	// some samples have an extraordinary number of intragenic events
	// if this is the case, we penalize intragenic events proportionately
	// consider only spliced events to compute the ratio, otherwise we would penalize TCR- and IG-rearranged tumors too much
	unsigned int spliced_events_in_same_gene = 0;
	unsigned int spliced_events_in_different_genes = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.spliced1 && fusion->second.spliced2) {
			if (fusion->second.gene1 == fusion->second.gene2)
				spliced_events_in_same_gene++;
			else
				spliced_events_in_different_genes++;
		}
	}
	// use reasonable defaut values, if sample size is too small
	if (spliced_events_in_same_gene + spliced_events_in_different_genes < 100) {
		spliced_events_in_same_gene = 0; // effectively disables penalty
		spliced_events_in_different_genes = 100;
	}

	// for each fusion, check if the observed number of supporting reads cannot be explained by random chance,
	// i.e. if the number is higher than expected given the number of fusion partners in both genes
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		// pick the gene with the most fusion partners
		float max_fusion_partners = max(
			10000.0 / fusion->second.gene1->exonic_length * max(fusion_partner_count[fusion->second.gene1]-1, 1),
			10000.0 / fusion->second.gene2->exonic_length * max(fusion_partner_count[fusion->second.gene2]-1, 1)
		);

		// calculate expected number of fusions (e-value)

		// the more reads there are in the rna.bam file, the more likely we find fusions supported by just a few reads (2-4)
		// the likelihood increases linearly, therefore we scale up the e-value proportionately to the number of mapped reads
		// for every 20 million reads, the scaling factor increases by 1 (this is an empirically determined value)
		fusion->second.evalue = max_fusion_partners * max(1.0, mapped_reads / 20000000.0 * pow(0.02, fusion->second.supporting_reads()-2));

		// intergenic and intragenic fusions are scored differently, because they have different frequencies
		if (fusion->second.is_intragenic()) {

			// events get a bonus based on their type
			// the bonus is proportionate to the frequency of the type
			// we multiply by 2.0 so that the overall effect is neutral, because there are two types (duplication vs. inversion)
			fusion->second.evalue *= 2.0 / (intragenic_duplications + intragenic_inversions);
			if (fusion->second.direction1 == UPSTREAM && fusion->second.direction2 == DOWNSTREAM)
				fusion->second.evalue *= intragenic_duplications;
			else if (fusion->second.direction1 == fusion->second.direction2)
				fusion->second.evalue *= intragenic_inversions;

			// the more fusion partners a gene has, the less likely a fusion is true (hence we multiply the e-value by max_fusion_partners)
			// but the likehood of a false positive decreases near-polynomially with the number of supporting reads
			if (fusion->second.supporting_reads() >= 1) {
				fusion->second.evalue *= pow(fusion->second.supporting_reads()-0.42, -2.11) * pow(10, -1.11);
				int spliced_distance = get_spliced_distance(fusion->second.contig1, fusion->second.breakpoint1, fusion->second.breakpoint2, fusion->second.gene1, exon_annotation_index);
				if (spliced_distance < 1000) {
					fusion->second.evalue *= pow(max(400, spliced_distance)/1000.0, -2);
					if (spliced_distance < 400)
						fusion->second.evalue *= pow(max(1, spliced_distance)/400.0, -4.58);
				}
			}

			// penalize intragenic events, if there are excessively many, i.e.
			// when the ratio of intragenic to intergenic events exceeds 0.25 (a value determined empirically from good-quality samples)
			fusion->second.evalue *= max(1.0, spliced_events_in_same_gene / 0.25 / spliced_events_in_different_genes);

		} else { // intergenic event

			if (fusion->second.supporting_reads() >= 1) {
				// the more fusion partners a gene has, the less likely a fusion is true (hence we multiply the e-value by max_fusion_partners)
				// but the likehood of a false positive decreases near-polynomially with the number of supporting reads
				fusion->second.evalue *= pow(fusion->second.supporting_reads()-0.73, -2.28) * pow(10, -1.75);

				if (fusion->second.is_read_through()) { // penalize read-through fusions
					fusion->second.evalue *= pow(max(1, fusion->second.breakpoint2 - fusion->second.breakpoint1)/400000.0, -0.63);
				} else if (fusion->second.contig1 == fusion->second.contig2 && fusion->second.breakpoint2 - fusion->second.breakpoint1 < 400000) { // penalize proximal events
					fusion->second.evalue *= pow(max(1, fusion->second.breakpoint2 - fusion->second.breakpoint1)/400000.0, -1.53);
				}
			}

		}

		// events get a bonus based on their location
		// the bonus is proportionate to the frequency of the events
		// we multiply by 4.0 so that the overall effect is neutral, because there are four possible locations (splice-site vs. intron vs. exon vs. mixed)
		// we always take max(spliced_breakpoints, ...), because spliced breakpoints should be the rarest or else the estimates are probably faulty
		fusion->second.evalue *= 4.0 / (spliced_breakpoints + exonic_breakpoints + intronic_breakpoints + exonic_intronic_breakpoints);
		if (fusion->second.spliced1 || fusion->second.spliced2)
			fusion->second.evalue *= spliced_breakpoints;
		else if (fusion->second.exonic1 && fusion->second.exonic2)
			fusion->second.evalue *= max(spliced_breakpoints, exonic_breakpoints);
		else if (!fusion->second.exonic1 && !fusion->second.exonic2)
			fusion->second.evalue *= max(spliced_breakpoints, intronic_breakpoints);
		else
			fusion->second.evalue *= max(spliced_breakpoints, exonic_intronic_breakpoints);
	}
}

unsigned int filter_relative_support(fusions_t& fusions, const float evalue_cutoff) {
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter != FILTER_none)
			continue; // fusion has already been filtered

		// throw away fusions which are expected to occur by random chance
		if (fusion->second.evalue < evalue_cutoff && // only keep fusions with good e-value
		    !(fusion->second.is_intragenic() && fusion->second.split_reads1 + fusion->second.split_reads2 == 0)) { // but ignore intragenic fusions only supported by discordant mates
			remaining++;
		} else {
			fusion->second.filter = FILTER_relative_support;
		}
	}
	return remaining;
}
