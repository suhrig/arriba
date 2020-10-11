#ifndef _RECOVER_BOTH_SPLICED_H
#define _RECOVER_BOTH_SPLICED_H 1

#include "common.hpp"
#include "read_stats.hpp"

using namespace std;

unsigned int recover_both_spliced(fusions_t& fusions, const chimeric_alignments_t& chimeric_alignments, const exon_annotation_index_t& exon_annotation_index, const coverage_t& coverage, const unsigned int max_fusions_to_recover, const float high_expression_quantile, const int max_exon_size, const unsigned int max_coverage);

#endif /* _RECOVER_BOTH_SPLICED_H */
