#ifndef FILTER_LOW_COVERAGE_VIRAL_CONTIGS_H
#define FILTER_LOW_COVERAGE_VIRAL_CONTIGS_H 1

#include <vector>
#include "common.hpp"
#include "read_stats.hpp"

using namespace std;

unsigned int filter_low_coverage_viral_contigs(chimeric_alignments_t& chimeric_alignments, const coverage_t& coverage, const vector<bool>& viral_contigs, const float min_covered_fraction, const float min_covered_bases);

#endif /* FILTER_LOW_COVERAGE_VIRAL_CONTIGS_H */
