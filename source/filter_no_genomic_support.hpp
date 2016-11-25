#ifndef _FILTER_NO_GENOMIC_SUPPORT_H
#define _FILTER_NO_GENOMIC_SUPPORT_H 1

#include <string>
#include "common.hpp"

using namespace std;

unsigned int mark_genomic_support(fusions_t& fusions, const string& genomic_breakpoints_file_path, const contigs_t& contigs, const int max_distance);

unsigned int filter_no_genomic_support(fusions_t& fusions, const float evalue_cutoff);

#endif /* _FILTER_NO_GENOMIC_SUPPORT_H */
