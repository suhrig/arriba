#ifndef _READ_STATS_H
#define _READ_STATS_H 1

#include <string>
#include <vector>
#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned long int count_mapped_reads(const string& bam_file_path, const vector<bool>& interesting_contigs);

bool estimate_mate_gap_distribution(const chimeric_alignments_t& chimeric_alignments, float& mate_gap_mean, float& mate_gap_stddev, const gene_annotation_index_t& gene_annotation_index, const exon_annotation_index_t& exon_annotation_index);

strandedness_t detect_strandedness(const chimeric_alignments_t& chimeric_alignments, const gene_annotation_index_t& gene_annotation_index);

#endif /* _READ_STATS_H */
