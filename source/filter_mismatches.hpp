#ifndef _FILTER_MISMATCHES_H
#define _FILTER_MISMATCHES_H 1

#include "common.hpp"

using namespace std;

unsigned int filter_mismatches(chimeric_alignments_t& chimeric_alignments, const assembly_t& assembly, const float mismatch_probability, const float pvalue_cutoff);

float estimate_mismatch_probability(const chimeric_alignments_t& chimeric_alignments, const assembly_t& assembly);

#endif /* _FILTER_MISMATCHES_H */
