#ifndef READ_CHIMERIC_ALIGNMENTS_H
#define READ_CHIMERIC_ALIGNMENTS_H 1

#include <string>
#include <vector>
#include "common.hpp"
#include "read_stats.hpp"

using namespace std;

unsigned int read_chimeric_alignments(const string& bam_file_path, const assembly_t& assembly, const string& assembly_file_path, chimeric_alignments_t& chimeric_alignments, unsigned long int& mapped_reads, vector<unsigned long int>& mapped_viral_reads_by_contig, coverage_t& coverage, contigs_t& contigs, vector<string>& original_contig_names, const string& interesting_contigs, const string& viral_contigs, const gene_annotation_index_t& gene_annotation_index, const bool separate_chimeric_bam_file, const bool is_rna_bam_file, const bool external_duplicate_marking, const unsigned int max_itd_length, int threads);

void assign_strands_from_strandedness(chimeric_alignments_t& chimeric_alignments, const strandedness_t strandedness);

unsigned int mark_multimappers(chimeric_alignments_t& chimeric_alignments);

#endif /* READ_CHIMERIC_ALIGNMENTS_H */

