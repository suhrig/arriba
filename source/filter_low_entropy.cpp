#include <cmath>
#include "sam.h"
#include "common.hpp"
#include "filter_low_entropy.hpp"
#include "filter_mismappers.hpp"

using namespace std;

unsigned int filter_low_entropy(chimeric_alignments_t& chimeric_alignments, const unsigned int kmer_length, const float kmer_content, const unsigned int max_itd_length) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		// all alignments that look like internal tandem duplications are checked for low entropy,
		// even if they have already been removed by previous filters, because low entropy regions
		// give rise to artifactual ITD alignments and the ITD filter would recover them, unless
		// they are marked as artifacts by the low_entropy filter
		bool is_internal_tandem_duplication = chimeric_alignment->second.size() == 3 && // split read
		                                      chimeric_alignment->second[SPLIT_READ].strand == chimeric_alignment->second[SUPPLEMENTARY].strand &&
		                                      chimeric_alignment->second[SPLIT_READ].contig == chimeric_alignment->second[SUPPLEMENTARY].contig &&
		                                      (
		                                      	chimeric_alignment->second[SPLIT_READ].strand == FORWARD &&
		                                      	chimeric_alignment->second[SPLIT_READ].start < chimeric_alignment->second[SUPPLEMENTARY].end &&
		                                      	chimeric_alignment->second[SPLIT_READ].start + ((int) max_itd_length) >= chimeric_alignment->second[SUPPLEMENTARY].end ||
		                                      	chimeric_alignment->second[SPLIT_READ].strand == REVERSE &&
		                                      	chimeric_alignment->second[SPLIT_READ].end > chimeric_alignment->second[SUPPLEMENTARY].start &&
		                                      	chimeric_alignment->second[SPLIT_READ].end <= chimeric_alignment->second[SUPPLEMENTARY].start + ((int) max_itd_length)
		                                      ); // alignments are oriented like a duplication

		if (!is_internal_tandem_duplication || chimeric_alignment->second.filter == FILTER_duplicates)
                	if (chimeric_alignment->second.filter != FILTER_none)
	                        continue; // read has already been filtered

		// look for recurrent k-mers in read sequence
		// if there are too many, discard the reads
		for (unsigned int mate = MATE1; mate <= MATE2; ++mate) {
			if (chimeric_alignment->second[mate].sequence.length() >= kmer_length) {

				// find out which part of the read aligns to the genome (is not clipped),
				// because k-mer content is computed for the whole read AND for the aligned segments individually
				unsigned int aligned_start1, aligned_end1, aligned_start2, aligned_end2;
				aligned_start1 = (chimeric_alignment->second[mate].cigar.operation(0) == BAM_CSOFT_CLIP) ? chimeric_alignment->second[mate].cigar.op_length(0) : 0;
				aligned_end1 = chimeric_alignment->second[mate].sequence.length();
				if (chimeric_alignment->second[mate].cigar.operation(chimeric_alignment->second[mate].cigar.size()-1) == BAM_CSOFT_CLIP)
					aligned_end1 -= chimeric_alignment->second[mate].cigar.op_length(chimeric_alignment->second[mate].cigar.size()-1);
				if (chimeric_alignment->second.size() == 3 && mate == SPLIT_READ) { // split read
					aligned_start2 = (chimeric_alignment->second[SUPPLEMENTARY].cigar.operation(0) == BAM_CSOFT_CLIP) ? chimeric_alignment->second[SUPPLEMENTARY].cigar.op_length(0) : 0;
					aligned_end2 = chimeric_alignment->second[SPLIT_READ].sequence.length();
					if (chimeric_alignment->second[SUPPLEMENTARY].cigar.operation(chimeric_alignment->second[SUPPLEMENTARY].cigar.size()-1) == BAM_CSOFT_CLIP)
						aligned_end2 -= chimeric_alignment->second[SUPPLEMENTARY].cigar.op_length(chimeric_alignment->second[SUPPLEMENTARY].cigar.size()-1);
					if (chimeric_alignment->second[SUPPLEMENTARY].strand != chimeric_alignment->second[SPLIT_READ].strand) {
						aligned_start2 = chimeric_alignment->second[SPLIT_READ].sequence.length() - aligned_start2;
						aligned_end2 = chimeric_alignment->second[SPLIT_READ].sequence.length() - aligned_end2;
						swap(aligned_start2, aligned_end2);
					}
				} else { // discordant mates
					aligned_start2 = aligned_start1;
					aligned_end2 = aligned_end1;
				}

				// create counters to keep track of the number of occurrences of every possible k-mer,
				// i.e., every possible combination of A, T, C, and G in a sequence of length <kmer_length>
				vector<unsigned int> kmer_count(pow(4, kmer_length));
				vector<unsigned int> kmer_count_aligned1(kmer_count.size());
				vector<unsigned int> kmer_count_aligned2(kmer_count.size());

				// determine thresholds that we consider "too many" identical k-mers in the same read
				unsigned int max_kmer_count = chimeric_alignment->second[mate].sequence.length() * kmer_content / kmer_length + 0.5;
				unsigned int max_kmer_count_aligned1 = (aligned_end1 - aligned_start1) * kmer_content / kmer_length + 0.5;
				unsigned int max_kmer_count_aligned2 = (aligned_end2 - aligned_start2) * kmer_content / kmer_length + 0.5;

				// when k-mers overlap, we should count them only once
				// this vector keeps track of the last position where a k-mer was found
				// new instances of k-mers are only counted, if they appear after the last k-mer
				vector<string::size_type> previous_kmer_pos(kmer_count.size());

				// count all different k-mers for each read
				for (string::size_type kmer_pos = 0; kmer_pos < chimeric_alignment->second[mate].sequence.length() - kmer_length; kmer_pos++) {

					kmer_as_int_t kmer_as_int = kmer_to_int(chimeric_alignment->second[mate].sequence, kmer_pos, kmer_length);

					// only count the k-mer if it does not overlap with a k-mer with identical sequence
					if (previous_kmer_pos[kmer_as_int] <= kmer_pos) {
						previous_kmer_pos[kmer_as_int] = kmer_pos + kmer_length;

						// update stats of given k-mer
						++kmer_count[kmer_as_int];
						if (kmer_pos+1 >= aligned_start1 && kmer_pos < aligned_end1) // k-mer is in aligned segment of mate1
							++kmer_count_aligned1[kmer_as_int];
						if (kmer_pos+1 >= aligned_start2 && kmer_pos < aligned_end2) // k-mer is in aligned segment of mate2
							++kmer_count_aligned2[kmer_as_int];

						// check if we crossed the k-mer count threshold
						if (kmer_count[kmer_as_int] >= max_kmer_count ||
						    kmer_count_aligned1[kmer_as_int] >= max_kmer_count_aligned1 ||
						    kmer_count_aligned2[kmer_as_int] >= max_kmer_count_aligned2) {
							chimeric_alignment->second.filter = FILTER_low_entropy;
							goto next_read;
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

