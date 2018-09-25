#ifndef _FILTER_UNINTERESTING_CONTIGS_H
#define _FILTER_UNINTERESTING_CONTIGS_H 1

#include "common.hpp"

using namespace std;

unsigned int filter_uninteresting_contigs(chimeric_alignments_t& chimeric_alignments, const contigs_t& contigs, const contigs_t& interesting_contigs);

#endif /* _FILTER_UNINTERESTING_CONTIGS_H */
