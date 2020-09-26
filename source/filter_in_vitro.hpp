#ifndef _FILTER_IN_VITRO_H
#define _FILTER_IN_VITRO_H 1

#include "common.hpp"
#include "read_stats.hpp"

using namespace std;

void find_top_expressed_genes(const chimeric_alignments_t& chimeric_alignments, const float high_expression_quantile, unordered_map<gene_t,unsigned int>& read_count_by_gene, unsigned int& high_expression_threshold);

unsigned int filter_in_vitro(fusions_t& fusions, const chimeric_alignments_t& chimeric_alignments, const float high_expression_quantile, const gene_annotation_index_t& gene_annotation_index, const coverage_t& coverage);

#endif /* _FILTER_IN_VITRO_H */
