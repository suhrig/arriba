#ifndef _FIND_GENOMIC_BREAKPOINTS_H
#define _FIND_GENOMIC_BREAKPOINTS_H 1

#include "annotation.hpp"
#include "common.hpp"

using namespace std;

unsigned int find_genomic_breakpoints(fusions_t& fusions, annotation_t& gene_annotation);

#endif /* _FIND_GENOMIC_BREAKPOINTS_H */
