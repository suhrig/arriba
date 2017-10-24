#ifndef _FILTER_RELATIVE_SUPPORT_H
#define _FILTER_RELATIVE_SUPPORT_H 1

#include "common.hpp"
#include "annotation.hpp"

using namespace std;

void estimate_expected_fusions(fusions_t& fusions, const unsigned long int mapped_reads);

unsigned int filter_relative_support(fusions_t& fusions, const float evalue_cutoff);

#endif /* _FILTER_RELATIVE_SUPPORT_H */
