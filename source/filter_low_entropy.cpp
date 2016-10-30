#include <cstring>
#include "sam.h"
#include "common.hpp"
#include "filter_low_entropy.hpp"

using namespace std;

unsigned int filter_low_entropy(chimeric_alignments_t& chimeric_alignments, const unsigned int kmer_length, const float kmer_content) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {

		if (!i->second.filters.empty())
			continue; // read has already been filtered

		// look for recurrent k-mers in read sequence
		// if there are too many, discard the reads
		for (unsigned int mate = MATE1; mate <= MATE2; ++mate) {
			if (i->second[mate].sequence.length() >= kmer_length) {

				// find out which part of the read aligns to the genome (is not clipped)
				unsigned int mate1_start, mate1_end, mate2_start, mate2_end;
				mate1_start = (i->second[mate].cigar.operation(0) == BAM_CSOFT_CLIP) ? i->second[mate].cigar.op_length(0) : 0;
				mate1_end = i->second[mate].sequence.length();
				if (i->second[mate].cigar.operation(i->second[mate].cigar.size()-1) == BAM_CSOFT_CLIP)
					mate1_end -= i->second[mate].cigar.op_length(i->second[mate].cigar.size()-1);
				if (i->second.size() == 3 && mate == SPLIT_READ) { // split read
					mate2_start = (i->second[SUPPLEMENTARY].cigar.operation(0) == BAM_CSOFT_CLIP) ? i->second[SUPPLEMENTARY].cigar.op_length(0) : 0;
					mate2_end = i->second[SPLIT_READ].sequence.length();
					if (i->second[SUPPLEMENTARY].cigar.operation(i->second[SUPPLEMENTARY].cigar.size()-1) == BAM_CSOFT_CLIP)
						mate2_end -= i->second[SUPPLEMENTARY].cigar.op_length(i->second[SUPPLEMENTARY].cigar.size()-1);
					if (i->second[SUPPLEMENTARY].strand != i->second[SPLIT_READ].strand) {
						mate2_start = i->second[SPLIT_READ].sequence.length() - mate2_start;
						mate2_end = i->second[SPLIT_READ].sequence.length() - mate2_end;
						swap(mate2_start, mate2_end);
					}
				} else { // discordant mates
					mate2_start = mate1_start;
					mate2_end = mate1_end;
				}

				// extract all possible k-mers from read
				const char* sequence = i->second[mate].sequence.data();
				for (unsigned int kmer_pos = 0; kmer_pos < i->second[mate].sequence.length() - kmer_length; kmer_pos++) {

					// count the number of occurrences of the given k-mer
					unsigned int kmer_count = 1;
					unsigned int kmer_count_aligned1 = (kmer_pos >= mate1_start && kmer_pos + kmer_length < mate1_end) ? 1 : 0;
					unsigned int kmer_count_aligned2 = (kmer_pos >= mate2_start && kmer_pos + kmer_length < mate2_end) ? 1 : 0;
					unsigned int max_kmer_count = i->second[mate].sequence.length() * kmer_content / kmer_length + 0.5;
					unsigned int max_kmer_count_aligned1 = (mate1_end - mate1_start) * kmer_content / kmer_length + 0.5;
					unsigned int max_kmer_count_aligned2 = (mate2_end - mate2_start) * kmer_content / kmer_length + 0.5;
					for (unsigned int pos = kmer_pos + kmer_length; pos < i->second[mate].sequence.length() - kmer_length;) {
						if (strncmp(sequence + kmer_pos, sequence + pos, kmer_length) == 0) {
							++kmer_count;
							if (pos+2 >= mate1_start && pos < mate1_end)
								++kmer_count_aligned1;
							if (pos+2 >= mate2_start && pos < mate2_end)
								++kmer_count_aligned2;

							if (kmer_count >= max_kmer_count || kmer_count_aligned1 >= max_kmer_count_aligned1 || kmer_count_aligned2 >= max_kmer_count_aligned2) {
								// a big fraction of the read consists of repetitive k-mers => remove it
								i->second.filters.insert(FILTERS.at("low_entropy"));
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

		next_read: NULL; // NULL is a dummy statement for the goto label
	}

	return remaining;
}

