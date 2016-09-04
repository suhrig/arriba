#include <cmath>
#include <set>
#include <string>
#include <unordered_map>
#include "common.hpp"
#include "annotation.hpp"
#include "filter_mismappers.hpp"

using namespace std;

typedef unordered_map< string, vector<string::size_type> > kmer_index_t;

bool align(int score, const string& read_sequence, string::size_type read_pos, const string& gene_sequence, string::size_type gene_pos, kmer_index_t& kmer_index, const unsigned int kmer_length, const int min_score, unsigned int segment_count) {

	const string::size_type max_extend_left = read_pos; // we must not extend the alignment beyond this left boundary

	for (/* read_pos comes from parameters */; read_pos + kmer_length < read_sequence.length() && read_pos + min_score < read_sequence.length() + score; ++read_pos, score--) {

		auto kmer_hits = kmer_index.find(read_sequence.substr(read_pos, kmer_length));
		if (kmer_hits == kmer_index.end())
			continue; // kmer not found in gene sequence

		for (auto kmer_hit = kmer_hits->second.begin(); kmer_hit != kmer_hits->second.end(); ++kmer_hit) {

			if (*kmer_hit < gene_pos)
				continue;

			int extended_score = score + kmer_length;
			if (extended_score >= min_score)
				return true;

			// extend match locally to the left
			// (if this is the first segment and we are not near the left end of the read / gene)
			string::size_type extended_left = *kmer_hit; // keeps track of how far we can extend to the left
			if (read_pos > max_extend_left + 2 && gene_pos > 2) {
				string::size_type extended_read_pos = read_pos - 1;
				string::size_type extended_gene_pos = *kmer_hit - 1;
				unsigned int mismatch_count = 0;
				while (extended_read_pos >= max_extend_left && extended_gene_pos >= 0) {

					if (read_sequence[extended_read_pos] == gene_sequence[extended_gene_pos]) {

						// read and gene sequences match => increase score and go the next base
						extended_score++;
						if (extended_score >= min_score)
							return true;
						extended_left--;

					} else { // there is a mismatch

						// allow one mismatch
						mismatch_count++;
						if (mismatch_count > 1)
							break;
						extended_score--;
						extended_left--;
					}
					extended_read_pos--;
					extended_gene_pos--;
				}
			}

			// extend match locally to the right
			{
				string::size_type extended_read_pos = read_pos + kmer_length;
				string::size_type extended_gene_pos = *kmer_hit + kmer_length;
				unsigned int mismatch_count = 0;
				while (extended_read_pos < read_sequence.length() && extended_gene_pos < gene_sequence.length()) {

					if (read_sequence[extended_read_pos] == gene_sequence[extended_gene_pos]) {

						// read and gene sequences match => increase score and go the next base
						extended_score++;
						if (extended_score >= min_score)
							return true;

					} else { // there is a mismatch

						// allow one mismatch
						mismatch_count++;
						if (mismatch_count > 1) {
							if (read_sequence.length() >= 30 && (segment_count == 1 || extended_left <= gene_pos+1)) { // when there is more than one mismatch, try to align the other part of the read (but only allow one intron/indel)
								if (align(extended_score, read_sequence, extended_read_pos, gene_sequence, extended_gene_pos, kmer_index, kmer_length, min_score, segment_count+1)) {
									return true;
								} else {
									break;
								}
							} else {
								break;
							}
						}
						extended_score--;
					}

					extended_read_pos++;
					extended_gene_pos++;
				}
			}
		}
	}

	// we only get here, if the read could not be aligned
	return false;
}

bool align_both_strands(const string& read_sequence, vector<kmer_index_t>& kmer_indices, gene_set_t& genes, annotation_t& gene_annotation, const unsigned int kmer_length, const int max_score, const float min_align_percent) {
	int min_score = min_align_percent * read_sequence.size();
	if (min_score > max_score)
		min_score = max_score;
	for (gene_set_t::iterator i = genes.begin(); i != genes.end(); ++i) {
		if (!gene_annotation[*i].sequence.empty()) {
			if (align(0, read_sequence, 0, gene_annotation[*i].sequence, 0, kmer_indices[*i], kmer_length, min_score, 1)) { // align on forward strand
				return true;
			} else { // align on reverse strand
				string reverse_complement;
				string original = read_sequence;
				dna_to_reverse_complement(original, reverse_complement);
				return align(0, reverse_complement, 0, gene_annotation[*i].sequence, 0, kmer_indices[*i], kmer_length, min_score, 1);
			}
		}
	}
	return false;
}

void count_mismappers(vector<mates_t*>& chimeric_alignments_list, unsigned int& mismappers, unsigned int& total_reads, unsigned int& supporting_reads) {
	for (auto i = chimeric_alignments_list.begin(); i != chimeric_alignments_list.end(); ++i) {
		if ((**i).filters.empty()) {
			total_reads++;
		} else if ((**i).filters.find(FILTERS.at("mismappers")) != (**i).filters.end()) {
			total_reads++;
			mismappers++;
			if (supporting_reads > 0)
				supporting_reads--;
		}
	}
}

unsigned int filter_mismappers(fusions_t& fusions, annotation_t& gene_annotation, float max_mismapper_fraction) {

	const unsigned int kmer_length = 8;
	const int max_score = 30; // maximum amount of nucleotides which need to align to consider alignment good
	float min_align_percent = 0.8; // allow ~1 mismatch for every six matches

	// make kmer indices from gene sequences
	vector<kmer_index_t> kmer_indices(gene_annotation.size()); // contains a kmer index for every gene sequence
	for (unsigned int gene = 0; gene < gene_annotation.size(); ++gene) {
		if (!gene_annotation[gene].sequence.empty()) {
			// store positions of kmers in hash
			for (string::size_type pos = 0; pos + kmer_length < gene_annotation[gene].sequence.length(); pos++)
				kmer_indices[gene][gene_annotation[gene].sequence.substr(pos, kmer_length)].push_back(pos);
		}
	}

	// align discordnat mate / clipped segment in gene of origin
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (!fusion->second.filters.empty())
			continue; // fusion has already been filtered

		if (fusion->second.gene1 == fusion->second.gene2)
			continue; // re-aligning the read only makes sense between different genes

		// re-align split reads
		vector<mates_t*> all_split_reads = fusion->second.split_read1_list;
		all_split_reads.insert(all_split_reads.end(), fusion->second.split_read2_list.begin(), fusion->second.split_read2_list.end());
		for (auto i = all_split_reads.begin(); i != all_split_reads.end(); ++i) {

			if ((**i).filters.empty()) { // read has not yet been filtered

				alignment_t& split_read = (**i)[SPLIT_READ]; // introduce alias for cleaner code
				if (split_read.strand == FORWARD) {
					if (align_both_strands(split_read.sequence.substr(0, split_read.preclipping()), kmer_indices, split_read.genes, gene_annotation, kmer_length, max_score, min_align_percent) || // clipped segment aligns to donor
					    align_both_strands((**i)[MATE1].sequence, kmer_indices, (**i)[SUPPLEMENTARY].genes, gene_annotation, kmer_length, max_score, min_align_percent)) { // non-spliced mate aligns to acceptor
						(**i).filters.insert(FILTERS.at("mismappers"));
					}
				} else { // split_read.strand == REVERSE
					if (align_both_strands(split_read.sequence.substr(split_read.sequence.length() - split_read.postclipping()), kmer_indices, split_read.genes, gene_annotation, kmer_length, max_score, min_align_percent) || // clipped segment aligns to donor
					    align_both_strands((**i)[MATE1].sequence, kmer_indices, (**i)[SUPPLEMENTARY].genes, gene_annotation, kmer_length, max_score, min_align_percent)) { // non-spliced mate aligns to acceptor
						(**i).filters.insert(FILTERS.at("mismappers"));
					}
				}
			}

		}

		// re-align discordant mates
		for (auto i = fusion->second.discordant_mate_list.begin(); i != fusion->second.discordant_mate_list.end(); ++i)
			if ((**i).filters.empty()) // read has not yet been filtered
				if ((**i).size() == 2) // discordant mates
					if (align_both_strands((**i)[MATE1].sequence, kmer_indices, (**i)[MATE2].genes, gene_annotation, kmer_length, max_score, min_align_percent) ||
					    align_both_strands((**i)[MATE2].sequence, kmer_indices, (**i)[MATE1].genes, gene_annotation, kmer_length, max_score, min_align_percent))
						(**i).filters.insert(FILTERS.at("mismappers"));

	}

	// discard all fusions with more than XX% mismappers
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (!fusion->second.filters.empty())
			continue; // fusion has already been filtered

		unsigned int total_reads = 0;
		unsigned int mismappers = 0;
		count_mismappers(fusion->second.split_read1_list, mismappers, total_reads, fusion->second.split_reads1);
		count_mismappers(fusion->second.split_read2_list, mismappers, total_reads, fusion->second.split_reads2);
		count_mismappers(fusion->second.discordant_mate_list, mismappers, total_reads, fusion->second.discordant_mates);

		// remove fusions with mostly mismappers
		if (mismappers > 0 && mismappers >= floor(max_mismapper_fraction * total_reads))
			fusion->second.filters.insert(FILTERS.at("mismappers"));
		else
			remaining++;

	}

	return remaining;
}


