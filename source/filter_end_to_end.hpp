#ifndef FILTER_END_TO_END_H
#define FILTER_END_TO_END_H 1

#include <vector>
#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned int filter_end_to_end_fusions(fusions_t& fusions, const exon_annotation_index_t& exon_annotation_index, const vector<bool>& viral_contigs);

#endif /* FILTER_END_TO_END_H */
