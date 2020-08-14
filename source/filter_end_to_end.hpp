#ifndef _FILTER_END_TO_END_H
#define _FILTER_END_TO_END_H 1

#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned int filter_end_to_end_fusions(fusions_t& fusions, const exon_annotation_index_t& exon_annotation_index, const contigs_t& contigs, const string& viral_contigs);

#endif /* _FILTER_END_TO_END_H */
