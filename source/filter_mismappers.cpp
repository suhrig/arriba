#include <cmath>
#include <set>
#include <string>
#include <unordered_map>
#include "common.hpp"
#include "annotation.hpp"
#include "filter_mismappers.hpp"

using namespace std;

typedef unsigned int kmer_as_int_t;
typedef unordered_map< kmer_as_int_t, vector<string::size_type> > kmer_index_t;
typedef vector<kmer_index_t> kmer_indices_t;

kmer_as_int_t kmer_to_int(const string& kmer, const string::size_type position, const char kmer_length) {
	unsigned int result = 0;
	for (char base = 0; base < kmer_length; ++base) {
		result = result<<2;
		switch (kmer.c_str()[position + base]) {
			case 'T': result += 0; break;
			case 'G': result += 1; break;
			case 'C': result += 2; break;
			default:  result += 3; break;
		}
	}
	return result;
}

bool align(int score, const string& read_sequence, string::size_type read_pos, const string& gene_sequence, const string::size_type gene_pos, const position_t gene_start, const position_t gene_end, const kmer_index_t& kmer_index, const char kmer_length, const int min_score, int max_deletions) {

	unsigned int skipped_bases = 0;

	for (/* read_pos comes from parameters */;
	     read_pos + kmer_length < read_sequence.length() && // don't run over end of read
	     read_pos + min_score <= read_sequence.length() + score + 2*kmer_length; // give up, when we can impossibly get above min_score, because we are near the end of the read
	                                                                             // 2*kmer_length takes into account that the score can improve, if we can extend to the left (up to kmer_length)
	     read_pos++, score--, skipped_bases++) { // if a base cannot be aligned, go to the next, but give -1 penalty and increase the number of skipped bases

		auto kmer_hits = kmer_index.find(kmer_to_int(read_sequence, read_pos, kmer_length));
		if (kmer_hits == kmer_index.end())
			continue; // kmer not found on given contig

		for (auto kmer_hit = lower_bound(kmer_hits->second.begin(), kmer_hits->second.end(), gene_start + gene_pos); kmer_hit != kmer_hits->second.end() && *kmer_hit < gene_end; ++kmer_hit) {

			int extended_score = score + kmer_length;
			if (read_pos == skipped_bases) // so far, all bases at the beginning of the read have been skipped
				extended_score += skipped_bases; // this effectively removes any penalties on leading mismatches (as in local alignment)
			if (extended_score >= min_score)
				return true;

			// extend match locally to the left
			{
				int extended_read_pos = read_pos - 1;
				int extended_gene_pos = *kmer_hit - gene_start - 1;
				unsigned int mismatch_count = 0;
				while (extended_read_pos >= read_pos - skipped_bases && // only align yet unaligned bases
				       extended_gene_pos >= 0) { // don't go beyond start of gene

					if (read_sequence[extended_read_pos] == gene_sequence[extended_gene_pos]) {

						// read and gene sequences match => increase score and go the next base
						if (read_pos == skipped_bases)
							extended_score += 1; // increase by 1, because leading skipped bases are not penalized in the main for-loop
						else
							extended_score += 2; // increase by 2, because intermittent skipped bases are already penalized in the main for-loop
						if (extended_score >= min_score)
							return true;

					} else { // there is a mismatch

						// allow one mismatch
						mismatch_count++;
						if (mismatch_count > 1)
							break;
						//extended_score--; // no need to penalize the score, because this is already done in the main for-loop
					}
					extended_read_pos--;
					extended_gene_pos--;
				}
			}

			// extend match locally to the right
			{
				int extended_read_pos = read_pos + kmer_length;
				int extended_gene_pos = *kmer_hit - gene_start + kmer_length;
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
						if (mismatch_count == 1) { // when there is more than one mismatch, do another k-mer lookup
							if (max_deletions > 0 && read_sequence.length() >= 30 && // do not allow too many deletions/introns and only if the read is reasonably long
							    align(extended_score, read_sequence, extended_read_pos, gene_sequence, extended_gene_pos, gene_start, gene_end, kmer_index, kmer_length, min_score, max_deletions-1)) {
								return true;
							}
						}
						extended_score--; // penalize mismatch
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

bool align_both_strands(const string& read_sequence, const position_t alignment_start, const position_t alignment_end, kmer_indices_t& kmer_indices, gene_set_t& genes, const char kmer_length, const float min_align_percent, int min_score) {
	min_score = min(min_score, (int) (min_align_percent * read_sequence.size() + 0.5));
	for (gene_set_t::iterator gene = genes.begin(); gene != genes.end(); ++gene) {

		if ((**gene).sequence.empty())
			continue;

		// in the case of intragenic events or overlapping genes,
		// both, the donor AND the acceptor gene overlap the breakpoint
		// => we do no alignment, because this would always discard the read
		if (alignment_start >= (**gene).start && alignment_start <= (**gene).end ||
		    alignment_end   >= (**gene).start && alignment_end   <= (**gene).end)
			continue;

		if (align(0, read_sequence, 0, (**gene).sequence, 0, (**gene).start, (**gene).end, kmer_indices[(**gene).contig], kmer_length, min_score, 1)) { // align on forward strand

			return true;
		} else { // align on reverse strand
			string reverse_complement;
			string original = read_sequence;
			dna_to_reverse_complement(original, reverse_complement);
			if (align(0, reverse_complement, 0, (**gene).sequence, 0, (**gene).start, (**gene).end, kmer_indices[(**gene).contig], kmer_length, min_score, 1))
				return true;
		}
	}
	return false;
}

void count_mismappers(vector<mates_t*>& chimeric_alignments_list, unsigned int& mismappers, unsigned int& total_reads, unsigned int& supporting_reads) {
	for (auto chimeric_alignment = chimeric_alignments_list.begin(); chimeric_alignment != chimeric_alignments_list.end(); ++chimeric_alignment) {
		if ((**chimeric_alignment).filter == NULL) {
			total_reads++;
		} else if ((**chimeric_alignment).filter == FILTERS.at("mismappers")) {
			total_reads++;
			mismappers++;
			if (supporting_reads > 0)
				supporting_reads--;
		}
	}
}

unsigned int filter_mismappers(fusions_t& fusions, const gene_annotation_t& gene_annotation, const contigs_t& contigs, const float max_mismapper_fraction) {

	const char kmer_length = 8; // must not be bigger than 16 or else conversion to int will fail
	const float min_align_percent = 0.8; // allow ~1 mismatch for every 10 matches
	const int min_score = 35; // consider this score or higher a match (even if less than min_align_percent match)

	// make kmer indices from gene sequences
	kmer_indices_t kmer_indices(contigs.size()); // contains a kmer index for every contig
	for (gene_annotation_t::const_iterator gene = gene_annotation.begin(); gene != gene_annotation.end(); ++gene) {
		if (!gene->sequence.empty()) {
			// store positions of kmers in hash
			for (string::size_type pos = 0; pos + kmer_length < gene->sequence.length(); pos++)
				if (gene->sequence[pos] != 'N') // don't index masked regions, as long stretches of N's cause alignment to take forever
					kmer_indices[gene->contig][kmer_to_int(gene->sequence, pos, kmer_length)].push_back(pos + gene->start);
		}
	}
	for (kmer_indices_t::iterator kmer_index = kmer_indices.begin(); kmer_index != kmer_indices.end(); ++kmer_index)
		for (kmer_index_t::iterator kmer_hits = kmer_index->begin(); kmer_hits != kmer_index->end(); ++kmer_hits)
			sort(kmer_hits->second.begin(), kmer_hits->second.end());

	// align discordnat mate / clipped segment in gene of origin
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.gene1 == fusion->second.gene2)
			continue; // re-aligning the read only makes sense between different genes

		if (fusion->second.filter != NULL)
			continue;

		// re-align split reads
		vector<mates_t*> all_split_reads;
		all_split_reads.insert(all_split_reads.end(), fusion->second.split_read1_list.begin(), fusion->second.split_read1_list.end());
		all_split_reads.insert(all_split_reads.end(), fusion->second.split_read2_list.begin(), fusion->second.split_read2_list.end());
		for (auto chimeric_alignment = all_split_reads.begin(); chimeric_alignment != all_split_reads.end(); ++chimeric_alignment) {

			if ((**chimeric_alignment).filter != NULL)
				continue; // read has already been filtered

			// introduce aliases for cleaner code
			alignment_t& split_read = (**chimeric_alignment)[SPLIT_READ];
			alignment_t& supplementary = (**chimeric_alignment)[SUPPLEMENTARY];
			alignment_t& mate1 = (**chimeric_alignment)[MATE1];

			if (split_read.strand == FORWARD) {
				if (align_both_strands(split_read.sequence.substr(0, split_read.preclipping()), supplementary.start, supplementary.end, kmer_indices, split_read.genes, kmer_length, min_align_percent, min_score) || // clipped segment aligns to donor
				    align_both_strands(mate1.sequence.substr(mate1.preclipping()), mate1.start, mate1.end, kmer_indices, supplementary.genes, kmer_length, min_align_percent, min_score)) { // non-spliced mate aligns to acceptor
					(**chimeric_alignment).filter = FILTERS.at("mismappers");
				}
			} else { // split_read.strand == REVERSE
				if (align_both_strands(split_read.sequence.substr(split_read.sequence.length() - split_read.postclipping()), supplementary.start, supplementary.end, kmer_indices, split_read.genes, kmer_length, min_align_percent, min_score) || // clipped segment aligns to donor
				    align_both_strands(mate1.sequence.substr(0, mate1.sequence.length() - mate1.postclipping()), mate1.start, mate1.end, kmer_indices, supplementary.genes, kmer_length, min_align_percent, min_score)) { // non-spliced mate aligns to acceptor
					(**chimeric_alignment).filter = FILTERS.at("mismappers");
				}
			}
		}

		// re-align discordant mates
		for (auto chimeric_alignment = fusion->second.discordant_mate_list.begin(); chimeric_alignment != fusion->second.discordant_mate_list.end(); ++chimeric_alignment) {
			if ((**chimeric_alignment).filter != NULL)
				continue; // read has already been filtered

			if ((**chimeric_alignment).size() == 2) { // discordant mates

				// introduce aliases for cleaner code
				alignment_t& mate1 = (**chimeric_alignment)[MATE1];
				alignment_t& mate2 = (**chimeric_alignment)[MATE2];

				if (align_both_strands(mate1.sequence, mate1.start, mate1.end, kmer_indices, mate2.genes, kmer_length, min_align_percent, min_score) ||
				    align_both_strands(mate2.sequence, mate2.start, mate2.end, kmer_indices, mate1.genes, kmer_length, min_align_percent, min_score)) {
					(**chimeric_alignment).filter = FILTERS.at("mismappers");
				}
			}
		}

	}

	// discard all fusions with more than XX% mismappers
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != NULL)
			continue; // fusion has already been filtered

		unsigned int total_reads = 0;
		unsigned int mismappers = 0;
		count_mismappers(fusion->second.split_read1_list, mismappers, total_reads, fusion->second.split_reads1);
		count_mismappers(fusion->second.split_read2_list, mismappers, total_reads, fusion->second.split_reads2);
		count_mismappers(fusion->second.discordant_mate_list, mismappers, total_reads, fusion->second.discordant_mates);

		// remove fusions with mostly mismappers
		if (mismappers > 0 && mismappers >= floor(max_mismapper_fraction * total_reads))
			fusion->second.filter = FILTERS.at("mismappers");
		else
			remaining++;

	}

	return remaining;
}

