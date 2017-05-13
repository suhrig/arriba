#ifndef _FILTER_NONEXPRESSED_H
#define _FILTER_NONEXPRESSED_H 1

#include <string>
#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned int filter_nonexpressed(fusions_t& fusions, const string& bam_file_path, const chimeric_alignments_t& chimeric_alignments, const exon_annotation_index_t& exon_annotation_index, const int max_mate_gap);

#endif /* _FILTER_NONEXPRESSED_H */
