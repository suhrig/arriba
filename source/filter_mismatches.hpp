#ifndef _FILTER_MISMATCHES_H
#define _FILTER_MISMATCHES_H 1

#include "common.hpp"

using namespace std;

unsigned int filter_mismatches(chimeric_alignments_t& chimeric_alignments, const assembly_t& assembly, const contigs_t& contigs, const string& interesting_contigs, const float mismatch_probability, const float pvalue_cutoff);

#endif /* _FILTER_MISMATCHES_H */
