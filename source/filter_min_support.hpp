#ifndef FILTER_MIN_SUPPORT_H
#define FILTER_MIN_SUPPORT_H 1

#include "common.hpp"

using namespace std;

// throw away fusions with few supporting reads
unsigned int filter_min_support(fusions_t& fusions, const int min_support);

#endif /* FILTER_MIN_SUPPORT_H */

