#ifndef FILTER_GENOMIC_SUPPORT_H
#define FILTER_GENOMIC_SUPPORT_H 1

#include <string>
#include <vector>
#include "common.hpp"
#include "read_stats.hpp"

using namespace std;

unsigned int mark_genomic_support(fusions_t& fusions, const string& genomic_breakpoints_file_path, const contigs_t& contigs, const int max_distance, const int max_itd_length);

void assign_confidence(fusions_t& fusions, const coverage_t& coverage);

unsigned int filter_no_genomic_support(fusions_t& fusions, const vector<bool>& viral_contigs);

unsigned int recover_genomic_support(fusions_t& fusions);

#endif /* FILTER_GENOMIC_SUPPORT_H */
