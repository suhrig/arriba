#ifndef _FILTER_LOW_ENTROPY_H
#define _FILTER_LOW_ENTROPY_H 1

#include "common.hpp"

using namespace std;

unsigned int filter_low_entropy(chimeric_alignments_t& chimeric_alignments, const unsigned int kmer_length, const float kmer_content, const unsigned int max_itd_length);

#endif /* _FILTER_LOW_ENTROPY_H */
