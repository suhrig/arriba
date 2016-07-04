#ifndef _FILTER_MISMAPPER_H
#define _FILTER_MISMAPPER_H 1

#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned int filter_mismappers(fusions_t& fusions, annotation_t& gene_annotation, float max_mismapper_fraction);

#endif /* _FILTER_MISMAPPERS_H */
