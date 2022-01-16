#ifndef FILTER_INTRAGENIC_BOTH_EXONIC_H
#define FILTER_INTRAGENIC_BOTH_EXONIC_H 1

#include "common.hpp"

using namespace std;

unsigned int filter_intragenic_both_exonic(fusions_t& fusions, const exon_annotation_index_t& exon_annotation_index, const float exonic_fraction);

#endif /* FILTER_INTRAGENIC_BOTH_EXONIC_H */
