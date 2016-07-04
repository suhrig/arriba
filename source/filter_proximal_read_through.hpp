#ifndef _FILTER_PROXIMAL_READ_THROUGH_H
#define _FILTER_PROXIMAL_READ_THROUGH_H 1

#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned int filter_proximal_read_through(chimeric_alignments_t& chimeric_alignments, annotation_t& gene_annotation, const unsigned int min_distance);

#endif /* _FILTER_PROXIMAL_READ_THROUGH_H */

