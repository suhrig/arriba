#ifndef _FILTER_MISMAPPER_H
#define _FILTER_MISMAPPER_H 1

#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned int filter_mismappers(fusions_t& fusions, const assembly_t& assembly, gene_annotation_t& gene_annotation, const exon_annotation_index_t& exon_annotation_index, const contigs_t& contigs, const float max_mismapper_fraction);

#endif /* _FILTER_MISMAPPERS_H */
