#ifndef _FILTER_PCR_FUSIONS_H
#define _FILTER_PCR_FUSIONS_H 1

#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned int filter_pcr_fusions(fusions_t& fusions, const chimeric_alignments_t& chimeric_alignments, const float high_expression_quantile, const gene_annotation_index_t& gene_annotation_index);

#endif /* _FILTER_PCR_FUSIONS_H */
