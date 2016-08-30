#include <cstring>
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

				// we won't look for k-mers beyond the position <kmer_end>
				unsigned int kmer_end = min(
					(int) (i->second[mate].sequence.length() - kmer_length), // the remaining segment of the read is shorter than <kmer_length>
					(int) (i->second[mate].sequence.length() * (1-kmer_content)) + 1 // the remaining fraction of the read is shorter than <kmer_content>
				);

				// extract all possible k-mers from read
				for (unsigned int kmer_pos = 0; kmer_pos < kmer_end; kmer_pos++) {
					const char* sequence = i->second[mate].sequence.data();
					char kmer[kmer_length];
					strncpy(kmer, sequence + kmer_pos, kmer_length);

					// count the number of occurrences of the given k-mer
					unsigned int kmer_count = 1;
					unsigned int max_kmer_count = i->second[mate].sequence.length() * kmer_content / kmer_length + 0.5;
					for (unsigned int pos = kmer_pos + kmer_length; pos < i->second[mate].sequence.length() - kmer_length;) {
						if (strncmp(kmer, sequence + pos, kmer_length) == 0) {
							pos += kmer_length;
							kmer_count++;

							if (kmer_count >= max_kmer_count) {
								// a big fraction of the read consists of repetitive k-mers => remove it
								i->second.filters.insert(FILTERS.at("low_entropy"));
								goto next_read;
							}
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

