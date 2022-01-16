#ifndef FILTER_RELATIVE_SUPPORT_H
#define FILTER_RELATIVE_SUPPORT_H 1

#include "common.hpp"
#include "annotation.hpp"

using namespace std;

void estimate_expected_fusions(fusions_t& fusions, const unsigned long int mapped_reads, const exon_annotation_index_t& exon_annotation_index);

unsigned int filter_relative_support(fusions_t& fusions, const float evalue_cutoff);

#endif /* FILTER_RELATIVE_SUPPORT_H */
