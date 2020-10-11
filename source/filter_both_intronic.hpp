#ifndef _FILTER_BOTH_INTRONIC_H
#define _FILTER_BOTH_INTRONIC_H 1

#include <vector>
#include "common.hpp"

using namespace std;

unsigned int filter_both_intronic(fusions_t& fusions, const vector<bool>& viral_contigs);

#endif /* _FILTER_BOTH_INTRONIC_H */
