#ifndef _FILTER_IN_VITRO_H
#define _FILTER_IN_VITRO_H 1

#include "common.hpp"

using namespace std;

unsigned int filter_in_vitro(fusions_t& fusions, const chimeric_alignments_t& chimeric_alignments, const float high_expression_quantile, const gene_annotation_index_t& gene_annotation_index);

#endif /* _FILTER_IN_VITRO_H */
