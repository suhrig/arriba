#ifndef _FILTER_HOMOPOLYMER_H
#define _FILTER_HOMOPOLYMER_H 1

#include "common.hpp"

using namespace std;

unsigned int filter_homopolymer(chimeric_alignments_t& chimeric_alignments, const unsigned int homopolymer_length, const exon_annotation_index_t& exon_annotation_index);

#endif /* _FILTER_HOMOPOLYMER_H */
