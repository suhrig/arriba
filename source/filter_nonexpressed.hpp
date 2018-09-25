#ifndef _FILTER_NONEXPRESSED_H
#define _FILTER_NONEXPRESSED_H 1

#include "common.hpp"
#include "annotation.hpp"
#include "read_stats.hpp"

using namespace std;

unsigned int filter_nonexpressed(fusions_t& fusions, const coverage_t& coverage, const exon_annotation_index_t& exon_annotation_index, const int max_mate_gap);

#endif /* _FILTER_NONEXPRESSED_H */
