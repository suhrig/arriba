#include <cmath>
#include <cstring>
#include <list>
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
	for (string::size_type pos = 0; pos + 2*kmer_length < small_gene_sequence.size(); pos += kmer_length) {

		if (matching_kmers * kmer_length + (small_gene_sequence.size() - pos) < small_gene->length() * max_identity_fraction)
			return false; // abort early, if there is no way we can possibly reach max_identity_fraction

		kmer_index_t::const_iterator kmer_hits = kmer_indices[big_gene->contig].find(kmer_to_int(small_gene_sequence, pos, kmer_length));
		if (kmer_hits != kmer_indices[big_gene->contig].end()) {
			for (auto kmer_hit = lower_bound(kmer_hits->second.begin(), kmer_hits->second.end(), big_gene->start); kmer_hit != kmer_hits->second.end() && *kmer_hit <= big_gene->end; ++kmer_hit) {
				if (small_gene->contig != big_gene->contig || *kmer_hit < small_gene->start || *kmer_hit > small_gene->end) {
					if (strncmp(assembly.at(big_gene->contig).c_str()+*kmer_hit+kmer_length, small_gene_sequence.c_str()+pos+kmer_length, extended_kmer_length) == 0) {
						matching_kmers++;
						if (matching_kmers * kmer_length >= small_gene->length() * max_identity_fraction)
							return true;
						break;
					}
				}
			}
		}
	}

	// if we get here, there weren't enough identical k-mers
	return false;
}

unsigned int filter_homologs(fusions_t& fusions, const kmer_indices_t& kmer_indices, const char kmer_length, const assembly_t& assembly, const float max_identity_fraction) {

	// select non-discarded fusions for better speed,
	// we need to iterate over them many times
	list<fusion_t*> remaining_fusions;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter == FILTER_none)
			remaining_fusions.push_front(&fusion->second);

	// discard fusion, if gene1 and gene2 are homologs
	for (auto fusion = remaining_fusions.begin(); fusion != remaining_fusions.end(); ++fusion) {

		if ((**fusion).filter != FILTER_none)
			continue;

		if (is_homolog((**fusion).gene1, (**fusion).gene2, kmer_indices, kmer_length, assembly, max_identity_fraction)) {

			(**fusion).filter = FILTER_homologs;

		} else {

			// sometimes a fusion between geneA and geneB leads to two predictions between
			// geneA and geneB as well as between geneA and a homolog of geneB due to mismapping reads
			// => look for other fusions concerning geneA and check if the fusion partners are homologs;
			//    if so, keep the one with more supporting reads or lower e-value
			for (auto other_fusion = next(fusion); other_fusion != remaining_fusions.end(); ++other_fusion) {

				if ((**other_fusion).filter != FILTER_none)
					continue;

				// check if geneA of fusion == geneA of other fusion
				// to determine which genes need to be checked for homology (geneB and geneC)
				gene_t homolog1, homolog2;
				if ((**fusion).gene1 == (**other_fusion).gene1 && (**fusion).breakpoint2 != (**other_fusion).breakpoint2) {
					homolog1 = (**fusion).gene2;
					homolog2 = (**other_fusion).gene2;
				} else if ((**fusion).gene1 == (**other_fusion).gene2 && (**fusion).breakpoint2 != (**other_fusion).breakpoint1) {
					homolog1 = (**fusion).gene2;
					homolog2 = (**other_fusion).gene1;
				} else if ((**fusion).gene2 == (**other_fusion).gene1 && (**fusion).breakpoint1 != (**other_fusion).breakpoint2) {
					homolog1 = (**fusion).gene1;
					homolog2 = (**other_fusion).gene2;
				} else if ((**fusion).gene2 == (**other_fusion).gene2 && (**fusion).breakpoint1 != (**other_fusion).breakpoint1) {
					homolog1 = (**fusion).gene1;
					homolog2 = (**other_fusion).gene1;
				} else
					continue; // the given fusions have no genes in common

				// find out which fusion has better alignments
				unsigned int anchor1 = ((**fusion).split_reads1 > 0) + ((**fusion).split_reads2 > 0) + ((**fusion).discordant_mates > 0);
				unsigned int anchor2 = ((**other_fusion).split_reads1 > 0) + ((**other_fusion).split_reads2 > 0) + ((**other_fusion).discordant_mates > 0);

				// check if the fusion partners geneB and geneC are homologs
				if (is_homolog(homolog1, homolog2, kmer_indices, kmer_length, assembly, max_identity_fraction)) {

					// other event must have poorer alignments or fewer reads or a worse e-value for us to consider its supporting reads to be mismappers
					if (anchor1 > anchor2 ||
					    anchor1 == anchor2 && (**fusion).supporting_reads() > (**other_fusion).supporting_reads() ||
					    anchor1 == anchor2 && (**fusion).supporting_reads() == (**other_fusion).supporting_reads() && (**fusion).evalue <= (**other_fusion).evalue) {
						(**other_fusion).filter = FILTER_homologs;
					} else {
						(**fusion).filter = FILTER_homologs;
						break;
					}
				}
			}
		}
	}

	// count fusions remaining after filtering
	unsigned int remaining = 0;
	for (auto fusion = remaining_fusions.begin(); fusion != remaining_fusions.end(); ++fusion)
		if ((**fusion).filter == FILTER_none)
			++remaining;
	return remaining;
}

