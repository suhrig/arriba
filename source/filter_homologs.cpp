#include <string>
#include "common.hpp"
#include "annotation.hpp"
#include "assembly.hpp"
#include "filter_mismappers.hpp"
#include "filter_homologs.hpp"

using namespace std;

unsigned int filter_homologs(fusions_t& fusions, const kmer_indices_t& kmer_indices, const char kmer_length, const assembly_t& assembly, const float max_identity_fraction) {

	// we look for kmers of length <kmer_length> + <extended_kmer_length> that are present in both genes
	const char extended_kmer_length = 8;

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != NULL)
			continue;

		if (fusion->second.gene1 == fusion->second.gene2) {
			remaining++;
			continue; // looking for homology only makes sense between different genes
		}

		gene_t small_gene = (fusion->second.gene1->length() > fusion->second.gene2->length()) ? fusion->second.gene2 : fusion->second.gene1;
		gene_t big_gene = (fusion->second.gene1->length() > fusion->second.gene2->length()) ? fusion->second.gene1 : fusion->second.gene2;
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

		if (matching_kmers * kmer_length >= small_gene->length() * max_identity_fraction) {
			fusion->second.filter = FILTERS.at("homologs");
		} else {
			remaining++;
		}
	}

	return remaining;
}

