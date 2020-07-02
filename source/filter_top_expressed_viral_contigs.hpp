#ifndef _FILTER_TOP_EXPRESSED_VIRAL_CONTIGS_H
#define _FILTER_TOP_EXPRESSED_VIRAL_CONTIGS_H 1

#include <string>
#include <vector>
#include "common.hpp"

using namespace std;

unsigned int filter_top_expressed_viral_contigs(chimeric_alignments_t& chimeric_alignments, unsigned int top_count, const contigs_t& contigs, const string& viral_contigs, const vector<unsigned long int>& mapped_viral_reads_by_contig, const assembly_t& assembly);

#endif /* _FILTER_TOP_EXPRESSED_VIRAL_CONTIGS_H */
