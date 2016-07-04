#ifndef _FIND_FUSIONS_H
#define _FIND_FUSIONS_H 1

#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned int find_fusions(chimeric_alignments_t& chimeric_alignments, fusions_t& fusions, annotation_index_t& exon_annotation_index);

#endif /* _FIND_FUSIONS_H */
