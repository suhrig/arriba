#ifndef _FILTER_BOTH_INTRONIC_H
#define _FILTER_BOTH_INTRONIC_H 1

#include <string>
#include "common.hpp"

using namespace std;

unsigned int filter_both_intronic(fusions_t& fusions, const contigs_t& contigs, const string& viral_contigs);

#endif /* _FILTER_BOTH_INTRONIC_H */
