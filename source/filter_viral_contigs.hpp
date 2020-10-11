#ifndef _FILTER_VIRAL_CONTIGS_H
#define _FILTER_VIRAL_CONTIGS_H 1

#include <vector>
#include "common.hpp"

using namespace std;

unsigned int filter_viral_contigs(chimeric_alignments_t& chimeric_alignments, const vector<bool>& viral_contigs);

#endif /* _FILTER_VIRAL_CONTIGS_H */
