#ifndef _FILTER_PCR_FUSIONS_H
#define _FILTER_PCR_FUSIONS_H 1

#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned int filter_pcr_fusions(fusions_t& fusions, const float max_pcr_fusion_score, const unsigned int max_exonic_breakpoints, const unsigned int max_partners_with_many_exonic_breakpoints, const unsigned int min_split_reads);

#endif /* _FILTER_PCR_FUSIONS_H */
