#ifndef _EXTRACT_READS_READ_THROUGH_FUSION_H
#define _EXTRACT_READS_READ_THROUGH_FUSION_H 1

#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"

using namespace std;

void extract_read_through_fusion(BGZF* read_through_fusions_file, bam1_t* forward_mate, bam1_t* reverse_mate, const gene_annotation_index_t& gene_annotation_index);

#endif /* _EXTRACT_READS_READ_THROUGH_FUSION_H */
