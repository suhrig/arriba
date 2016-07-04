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
#include "filter_promiscuous_genes.hpp"

using namespace std;

unsigned long int count_mapped_reads(const string& bam_file_path, const vector<bool>& interesting_contigs) {
	unsigned long int result = 0;
	bam_index_t* bam_index = bam_index_load(bam_file_path.c_str());
	for (unsigned int i = 0; i < interesting_contigs.size(); ++i) {
        	unsigned long int mapped, unmapped;
		hts_idx_get_stat(bam_index, i, &mapped, &unmapped);
		if (interesting_contigs[i]) // only count reads on interesting contigs
			result += mapped;
	}
	return result;
}

void estimate_expected_fusions(fusions_t& fusions, const annotation_t& gene_annotation, const unsigned long int mapped_reads) {

	// find all fusion partners for each gene
	vector< set<gene_t> > fusion_partners(gene_annotation.size());
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
		if (i->second.filters.empty()) {
			if (!i->second.overlap_duplicate1)
				fusion_partners[i->second.gene2].insert(i->second.gene1);
			if (!i->second.overlap_duplicate2)
				fusion_partners[i->second.gene1].insert(i->second.gene2);
		}
	}

	// count the number of fusion partners for each gene
	// fusions with genes that have more fusion partners are ignored
	vector<int> fusion_partner_count(gene_annotation.size());
	for (gene_t i = 0; i < gene_annotation.size(); ++i) {
		for (set<gene_t>::iterator j = fusion_partners[i].begin(); j != fusion_partners[i].end(); ++j) {
			if (fusion_partners[i].size() >= fusion_partners[*j].size()) {
				fusion_partner_count[i]++;
			}
		}
	}

	// for each fusion, check if the observed number of supporting reads cannot be explained by random chance,
	// i.e. if the number is higher than expected given the number of fusion partners in both genes
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		// pick the gene with the most fusion partners
		float max_fusion_partners = max(
			10000.0 / gene_annotation[i->second.gene1].exonic_length * max(fusion_partner_count[i->second.gene1]-1, 0),
			10000.0 / gene_annotation[i->second.gene2].exonic_length * max(fusion_partner_count[i->second.gene2]-1, 0)
		);
 
		unsigned int supporting_reads = i->second.split_reads1 + i->second.split_reads2 + i->second.discordant_mates;

		// calculate expected number of fusions (e-value)

		// the more reads there are in the rna.bam file, the more likely we find fusions supported by just a few reads (2-4)
		// the likelihood increases linearly, therefore we scale up the e-value proportionately to the number of mapped reads
		// for every 20 million reads, the scaling factor increases by 1 (this is an empirically determined value)
		i->second.evalue = max_fusion_partners * max(1.0, mapped_reads / 20000000.0 * pow(0.02, supporting_reads-2));

		// the more fusion partners a gene has, the less likely a fusion is true (hence we multiply the e-value by max_fusion_partners)
		// but the likehood of a false positive decreases near-exponentially with the number of supporting reads (hence the mutiply by 0.2*0.4^(supporting_reads-2) )
		if (supporting_reads > 1)
			i->second.evalue *= 0.02 * pow(0.4, supporting_reads-2);
	}

}

unsigned int filter_promiscuous_genes(fusions_t& fusions, const float evalue_cutoff) {
	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
		if (!i->second.filters.empty())
			continue; // fusion has already been filtered

		// throw away fusions which are expected to occur by random chance
		if (i->second.evalue < evalue_cutoff) {
			remaining++;
		} else {
			i->second.filters.insert(FILTERS.at("promiscuous_genes"));
		}
	}
	return remaining;
}
