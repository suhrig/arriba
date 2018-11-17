#include <cstring>
#include "sam.h"
#include "common.hpp"
#include "filter_low_entropy.hpp"

using namespace std;

unsigned int filter_low_entropy(chimeric_alignments_t& chimeric_alignments, const unsigned int kmer_length, const float kmer_content) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != NULL)
			continue; // read has already been filtered

		// look for recurrent k-mers in read sequence
		// if there are too many, discard the reads
		for (unsigned int mate = MATE1; mate <= MATE2; ++mate) {
			if (chimeric_alignment->second[mate].sequence.length() >= kmer_length) {

				// find out which part of the read aligns to the genome (is not clipped)
				unsigned int mate1_start, mate1_end, mate2_start, mate2_end;
				mate1_start = (chimeric_alignment->second[mate].cigar.operation(0) == BAM_CSOFT_CLIP) ? chimeric_alignment->second[mate].cigar.op_length(0) : 0;
				mate1_end = chimeric_alignment->second[mate].sequence.length();
				if (chimeric_alignment->second[mate].cigar.operation(chimeric_alignment->second[mate].cigar.size()-1) == BAM_CSOFT_CLIP)
					mate1_end -= chimeric_alignment->second[mate].cigar.op_length(chimeric_alignment->second[mate].cigar.size()-1);
				if (chimeric_alignment->second.size() == 3 && mate == SPLIT_READ) { // split read
					mate2_start = (chimeric_alignment->second[SUPPLEMENTARY].cigar.operation(0) == BAM_CSOFT_CLIP) ? chimeric_alignment->second[SUPPLEMENTARY].cigar.op_length(0) : 0;
					mate2_end = chimeric_alignment->second[SPLIT_READ].sequence.length();
					if (chimeric_alignment->second[SUPPLEMENTARY].cigar.operation(chimeric_alignment->second[SUPPLEMENTARY].cigar.size()-1) == BAM_CSOFT_CLIP)
						mate2_end -= chimeric_alignment->second[SUPPLEMENTARY].cigar.op_length(chimeric_alignment->second[SUPPLEMENTARY].cigar.size()-1);
					if (chimeric_alignment->second[SUPPLEMENTARY].strand != chimeric_alignment->second[SPLIT_READ].strand) {
						mate2_start = chimeric_alignment->second[SPLIT_READ].sequence.length() - mate2_start;
						mate2_end = chimeric_alignment->second[SPLIT_READ].sequence.length() - mate2_end;
						swap(mate2_start, mate2_end);
					}
				} else { // discordant mates
					mate2_start = mate1_start;
					mate2_end = mate1_end;
				}

				// extract all possible k-mers from read
				const char* sequence = chimeric_alignment->second[mate].sequence.data();
				for (unsigned int kmer_pos = 0; kmer_pos < chimeric_alignment->second[mate].sequence.length() - kmer_length; kmer_pos++) {

					// count the number of occurrences of the given k-mer
					unsigned int kmer_count = 1;
					unsigned int kmer_count_aligned1 = (kmer_pos >= mate1_start && kmer_pos + kmer_length < mate1_end) ? 1 : 0;
					unsigned int kmer_count_aligned2 = (kmer_pos >= mate2_start && kmer_pos + kmer_length < mate2_end) ? 1 : 0;
					unsigned int max_kmer_count = chimeric_alignment->second[mate].sequence.length() * kmer_content / kmer_length + 0.5;
					unsigned int max_kmer_count_aligned1 = (mate1_end - mate1_start) * kmer_content / kmer_length + 0.5;
					unsigned int max_kmer_count_aligned2 = (mate2_end - mate2_start) * kmer_content / kmer_length + 0.5;
					for (unsigned int pos = kmer_pos + kmer_length; pos < chimeric_alignment->second[mate].sequence.length() - kmer_length;) {
						if (strncmp(sequence + kmer_pos, sequence + pos, kmer_length) == 0) {
							++kmer_count;
							if (pos+2 >= mate1_start && pos < mate1_end)
								++kmer_count_aligned1;
							if (pos+2 >= mate2_start && pos < mate2_end)
								++kmer_count_aligned2;

							if (kmer_count >= max_kmer_count || kmer_count_aligned1 >= max_kmer_count_aligned1 || kmer_count_aligned2 >= max_kmer_count_aligned2) {
								// a big fraction of the read consists of repetitive k-mers => remove it
								chimeric_alignment->second.filter = FILTERS.at("low_entropy");
								goto next_read;
							}

							pos += kmer_length;
						} else {
							pos++;
						}
					}
				}
			}
		}

		++remaining;

		next_read: continue;
	}

	return remaining;
}

