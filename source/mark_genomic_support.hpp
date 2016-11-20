#ifndef _MARK_GENOMIC_SUPPORT_H
#define _MARK_GENOMIC_SUPPORT_H 1

#include <string>
#include "common.hpp"

using namespace std;

unsigned int mark_genomic_support(fusions_t& fusions, const string& genomic_breakpoints_file_path, const contigs_t& contigs, const int max_distance);

#endif /* _MARK_GENOMIC_SUPPORT_H */
