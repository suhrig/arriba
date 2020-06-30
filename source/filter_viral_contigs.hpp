#ifndef _FILTER_VIRAL_CONTIGS_H
#define _FILTER_VIRAL_CONTIGS_H 1

#include <string>
#include "common.hpp"

using namespace std;

unsigned int filter_viral_contigs(chimeric_alignments_t& chimeric_alignments, const contigs_t& contigs, const string& viral_contigs);

#endif /* _FILTER_VIRAL_CONTIGS_H */
