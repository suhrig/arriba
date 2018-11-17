#include <cmath>
#include <string>
#include "common.hpp"
#include "annotation.hpp"
#include "assembly.hpp"
#include "filter_mismappers.hpp"
#include "filter_homologs.hpp"

using namespace std;

bool is_homolog(const gene_t gene1, const gene_t gene2, const kmer_indices_t& kmer_indices, const char kmer_length, const assembly_t& assembly, const float max_identity_fraction) {

	// we look for kmers of length <kmer_length> + <extended_kmer_length> that are present in both genes
	const char extended_kmer_length = 8;

	// looking for homology only makes sense between different genes
	if (gene1 == gene2)
		return false;

	// find the smaller of the two genes
	gene_t small_gene = gene1;
	gene_t big_gene = gene2;
	if (small_gene->length() > big_gene->length())
		swap(small_gene, big_gene);

	// genes must not overlap, otherwise there is for sure going to be sequence similarity
	if (small_gene->contig == big_gene->contig &&
	    (small_gene->start >= big_gene->start && small_gene->start <= big_gene->end ||
	     small_gene->end   >= big_gene->start && small_gene->end   <= big_gene->end))
		return false;

	// retrieve sequence of smaller gene
	string small_gene_sequence = assembly.at(small_gene->contig).substr(small_gene->start, small_gene->length());
	if (small_gene->strand != big_gene->strand)
		small_gene_sequence = dna_to_reverse_complement(small_gene_sequence);

	// count how many k-mers of the small gene can be found in the big gene
	unsigned int matching_kmers = 0;
	for (string::size_type pos = 0; pos + kmer_length < small_gene_sequence.size(); pos += kmer_length) {
		kmer_index_t::const_iterator kmer_hits = kmer_indices[big_gene->contig].find(kmer_to_int(small_gene_sequence, pos, kmer_length));
		if (kmer_hits != kmer_indices[big_gene->contig].end()) {
			for (auto kmer_hit = lower_bound(kmer_hits->second.begin(), kmer_hits->second.end(), big_gene->start); kmer_hit != kmer_hits->second.end() && *kmer_hit <= big_gene->end; ++kmer_hit) {
				if (small_gene->contig != big_gene->contig || *kmer_hit < small_gene->start || *kmer_hit > small_gene->end) {
					if (assembly.at(big_gene->contig).substr(*kmer_hit+kmer_length, extended_kmer_length) == small_gene_sequence.substr(pos+kmer_length, extended_kmer_length)) {
						matching_kmers++;
						break;
					}
				}
			}
		}
	}

	return (matching_kmers * kmer_length >= small_gene->length() * max_identity_fraction);
}

unsigned int filter_homologs(fusions_t& fusions, const kmer_indices_t& kmer_indices, const char kmer_length, const assembly_t& assembly, const float max_identity_fraction) {

	// discard fusion, if gene1 and gene2 are homologs
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != NULL)
			continue;

		if (is_homolog(fusion->second.gene1, fusion->second.gene2, kmer_indices, kmer_length, assembly, max_identity_fraction)) {

			fusion->second.filter = FILTERS.at("homologs");
			continue;

		} else {

			remaining++;

			// sometimes a fusion between geneA and geneB leads to two predictions between
			// geneA and geneB as well as between geneA and a homolog of geneB due to mismapping reads
			// => look for other fusions concerning geneA and check if the fusion partners are homologs;
			//    if so, keep the one with more supporting reads or lower e-value
			for (fusions_t::iterator other_fusion = fusions.begin(); other_fusion != fusions.end(); ++other_fusion) {

				if (other_fusion->second.filter != NULL)
					continue;

				// check if geneA of fusion == geneA of other fusion
				// to determine which genes need to be checked for homology (geneB and geneC)
				gene_t homolog1, homolog2;
				if (fusion->second.gene1 == other_fusion->second.gene1 && fusion->second.breakpoint2 != other_fusion->second.breakpoint2) {
					homolog1 = fusion->second.gene2;
					homolog2 = other_fusion->second.gene2;
				} else if (fusion->second.gene1 == other_fusion->second.gene2 && fusion->second.breakpoint2 != other_fusion->second.breakpoint1) {
					homolog1 = fusion->second.gene2;
					homolog2 = other_fusion->second.gene1;
				} else if (fusion->second.gene2 == other_fusion->second.gene1 && fusion->second.breakpoint1 != other_fusion->second.breakpoint2) {
					homolog1 = fusion->second.gene1;
					homolog2 = other_fusion->second.gene2;
				} else if (fusion->second.gene2 == other_fusion->second.gene2 && fusion->second.breakpoint1 != other_fusion->second.breakpoint1) {
					homolog1 = fusion->second.gene1;
					homolog2 = other_fusion->second.gene1;
				} else
					continue; // the given fusions have no genes in common

				// find out which fusion has better alignments
				unsigned int anchor1 = (fusion->second.split_reads1 > 0) + (fusion->second.split_reads2 > 0) + (fusion->second.discordant_mates > 0);
				unsigned int anchor2 = (other_fusion->second.split_reads1 > 0) + (other_fusion->second.split_reads2 > 0) + (other_fusion->second.discordant_mates > 0);

				// other event must have poorer alignments or fewer reads or a worse e-value for us to consider its supporting reads to be mismappers
				if (anchor1 > anchor2 ||
				    anchor1 == anchor2 && fusion->second.supporting_reads() > other_fusion->second.supporting_reads() ||
				    fusion->second.supporting_reads() == other_fusion->second.supporting_reads() && fusion->second.evalue <= other_fusion->second.evalue) {

					// check if the fusion partners geneB and geneC are homologs
					if (is_homolog(homolog1, homolog2, kmer_indices, kmer_length, assembly, max_identity_fraction))
						other_fusion->second.filter = FILTERS.at("homologs");

				}
			}
		}
	}

	return remaining;
}

