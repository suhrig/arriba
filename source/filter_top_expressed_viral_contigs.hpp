#ifndef FILTER_TOP_EXPRESSED_VIRAL_CONTIGS_H
#define FILTER_TOP_EXPRESSED_VIRAL_CONTIGS_H 1

#include <vector>
#include "common.hpp"

using namespace std;

unsigned int filter_top_expressed_viral_contigs(chimeric_alignments_t& chimeric_alignments, unsigned int top_count, const vector<bool>& viral_contigs, const vector<bool>& interesting_contigs, const vector<unsigned long int>& mapped_viral_reads_by_contig, const assembly_t& assembly);

#endif /* FILTER_TOP_EXPRESSED_VIRAL_CONTIGS_H */
