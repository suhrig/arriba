#ifndef _FILTER_LOW_COVERAGE_VIRAL_CONTIGS_H
#define _FILTER_LOW_COVERAGE_VIRAL_CONTIGS_H 1

#include <string>
#include "common.hpp"
#include "read_stats.hpp"

using namespace std;

unsigned int filter_low_coverage_viral_contigs(chimeric_alignments_t& chimeric_alignments, const coverage_t& coverage, const contigs_t& contigs, const string& viral_contigs, const float min_covered_fraction);

#endif /* _FILTER_LOW_COVERAGE_VIRAL_CONTIGS_H */
