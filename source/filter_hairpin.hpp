#ifndef FILTER_HAIRPIN_H
#define FILTER_HAIRPIN_H 1

#include "common.hpp"

using namespace std;

unsigned int filter_hairpin(chimeric_alignments_t& chimeric_alignments, exon_annotation_index_t& exon_annotation_index, const int max_mate_gap);

#endif /* FILTER_HAIRPIN_H */
