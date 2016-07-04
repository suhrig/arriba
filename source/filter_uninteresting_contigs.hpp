#ifndef _FILTER_UNINTERESTING_CONTIGS_H
#define _FILTER_UNINTERESTING_CONTIGS_H 1

#include <vector>
#include "common.hpp"

using namespace std;

unsigned int filter_uninteresting_contigs(chimeric_alignments_t& chimeric_alignments, vector<bool>& interesting_contigs);

#endif /* _FILTER_UNINTERESTING_CONTIGS_H */
