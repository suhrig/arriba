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

bool align(int score, const string& read_sequence, string::size_type read_pos, const string& gene_sequence, const string::size_type gene_pos, const position_t gene_start, const position_t gene_end, const kmer_index_t& kmer_index, const char kmer_length, const int min_score, unsigned int segment_count) {

	const string::size_type max_extend_left = read_pos; // we must not extend the alignment beyond this left boundary

	for (/* read_pos comes from parameters */; read_pos + kmer_length < read_sequence.length() && read_pos + min_score < read_sequence.length() + score; ++read_pos, score--) {

		auto kmer_hits = kmer_index.find(kmer_to_int(read_sequence, read_pos, kmer_length));
		if (kmer_hits == kmer_index.end())
			continue; // kmer not found on given contig

		for (auto kmer_hit = lower_bound(kmer_hits->second.begin(), kmer_hits->second.end(), gene_start + gene_pos); kmer_hit != kmer_hits->second.end() && *kmer_hit < gene_end; ++kmer_hit) {

			int extended_score = score + kmer_length;
			if (extended_score >= min_score)
				return true;

			// extend match locally to the left
			// (if this is the first segment and we are not near the left end of the read / gene)
			string::size_type extended_left = *kmer_hit - gene_start; // keeps track of how far we can extend to the left
			if (read_pos > max_extend_left + 2 && gene_pos > 2) {
				string::size_type extended_read_pos = read_pos - 1;
				string::size_type extended_gene_pos = *kmer_hit - gene_start - 1;
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
				string::size_type extended_gene_pos = *kmer_hit - gene_start + kmer_length;
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
								if (align(extended_score, read_sequence, extended_read_pos, gene_sequence, extended_gene_pos, gene_start, gene_end, kmer_index, kmer_length, min_score, segment_count+1)) {
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

bool align_both_strands(const string& read_sequence, const position_t alignment_start, const position_t alignment_end, kmer_indices_t& kmer_indices, gene_set_t& genes, const char kmer_length, const int max_score, const float min_align_percent) {
	int min_score = min_align_percent * read_sequence.size();
	if (min_score > max_score)
		min_score = max_score;
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
	for (auto i = chimeric_alignments_list.begin(); i != chimeric_alignments_list.end(); ++i) {
		if ((**i).filter == NULL) {
			total_reads++;
		} else if ((**i).filter == FILTERS.at("mismappers")) {
			total_reads++;
			mismappers++;
			if (supporting_reads > 0)
				supporting_reads--;
		}
	}
}

unsigned int filter_mismappers(fusions_t& fusions, const gene_annotation_t& gene_annotation, const contigs_t& contigs, const float max_mismapper_fraction) {

	const char kmer_length = 8; // must not be bigger than 32 or else conversion to int will fail
	const int max_score = 30; // maximum amount of nucleotides which need to align to consider alignment good
	const float min_align_percent = 0.8; // allow ~1 mismatch for every six matches

	// make kmer indices from gene sequences
	kmer_indices_t kmer_indices(contigs.size()); // contains a kmer index for every contig
	for (gene_annotation_t::const_iterator gene = gene_annotation.begin(); gene != gene_annotation.end(); ++gene) {
		if (!gene->sequence.empty()) {
			// store positions of kmers in hash
			for (string::size_type pos = 0; pos + kmer_length < gene->sequence.length(); pos++)
				kmer_indices[gene->contig][kmer_to_int(gene->sequence, pos, kmer_length)].push_back(pos + gene->start);
		}
	}
	for (kmer_indices_t::iterator kmer_index = kmer_indices.begin(); kmer_index != kmer_indices.end(); ++kmer_index)
		for (kmer_index_t::iterator kmer_hits = kmer_index->begin(); kmer_hits != kmer_index->end(); ++kmer_hits)
			sort(kmer_hits->second.begin(), kmer_hits->second.end());

	// align discordnat mate / clipped segment in gene of origin
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != NULL)
			continue; // fusion has already been filtered

		if (fusion->second.gene1 == fusion->second.gene2)
			continue; // re-aligning the read only makes sense between different genes

		// re-align split reads
		vector<mates_t*> all_split_reads;
		all_split_reads.insert(all_split_reads.end(), fusion->second.split_read1_list.begin(), fusion->second.split_read1_list.end());
		all_split_reads.insert(all_split_reads.end(), fusion->second.split_read2_list.begin(), fusion->second.split_read2_list.end());
		for (auto i = all_split_reads.begin(); i != all_split_reads.end(); ++i) {

			if ((**i).filter != NULL)
				continue; // read has already been filtered

			// introduce aliases for cleaner code
			alignment_t& split_read = (**i)[SPLIT_READ];
			alignment_t& supplementary = (**i)[SUPPLEMENTARY];
			alignment_t& mate1 = (**i)[MATE1];

			if (split_read.strand == FORWARD) {
				if (align_both_strands(split_read.sequence.substr(0, split_read.preclipping()), supplementary.start, supplementary.end, kmer_indices, split_read.genes, kmer_length, max_score, min_align_percent) || // clipped segment aligns to donor
				    align_both_strands(mate1.sequence.substr(mate1.preclipping()), mate1.start, mate1.end, kmer_indices, supplementary.genes, kmer_length, max_score, min_align_percent)) { // non-spliced mate aligns to acceptor
					(**i).filter = FILTERS.at("mismappers");
				}
			} else { // split_read.strand == REVERSE
				if (align_both_strands(split_read.sequence.substr(split_read.sequence.length() - split_read.postclipping()), supplementary.start, supplementary.end, kmer_indices, split_read.genes, kmer_length, max_score, min_align_percent) || // clipped segment aligns to donor
				    align_both_strands(mate1.sequence.substr(0, mate1.sequence.length() - mate1.postclipping()), mate1.start, mate1.end, kmer_indices, supplementary.genes, kmer_length, max_score, min_align_percent)) { // non-spliced mate aligns to acceptor
					(**i).filter = FILTERS.at("mismappers");
				}
			}
		}

		// re-align discordant mates
		for (auto i = fusion->second.discordant_mate_list.begin(); i != fusion->second.discordant_mate_list.end(); ++i) {
			if ((**i).filter != NULL)
				continue; // read has already been filtered

			if ((**i).size() == 2) { // discordant mates

				// introduce aliases for cleaner code
				alignment_t& mate1 = (**i)[MATE1];
				alignment_t& mate2 = (**i)[MATE2];

				if (align_both_strands(mate1.sequence, mate1.start, mate1.end, kmer_indices, mate2.genes, kmer_length, max_score, min_align_percent) ||
				    align_both_strands(mate2.sequence, mate2.start, mate2.end, kmer_indices, mate1.genes, kmer_length, max_score, min_align_percent))
					(**i).filter = FILTERS.at("mismappers");
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


