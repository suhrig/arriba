#ifndef FILTER_NO_COVERAGE_H
#define FILTER_NO_COVERAGE_H 1

#include "common.hpp"
#include "annotation.hpp"
#include "read_stats.hpp"

using namespace std;

unsigned int filter_no_coverage(fusions_t& fusions, const coverage_t& coverage, const exon_annotation_index_t& exon_annotation_index);

#endif /* FILTER_NO_COVERAGE_H */
