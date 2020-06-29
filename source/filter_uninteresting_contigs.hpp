#ifndef _FILTER_UNINTERESTING_CONTIGS_H
#define _FILTER_UNINTERESTING_CONTIGS_H 1

#include <string>
#include "common.hpp"

using namespace std;

unsigned int filter_uninteresting_contigs(chimeric_alignments_t& chimeric_alignments, const contigs_t& contigs, const string& interesting_contigs);

#endif /* _FILTER_UNINTERESTING_CONTIGS_H */
