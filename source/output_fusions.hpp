#ifndef _OUTPUT_FUSIONS_H
#define _OUTPUT_FUSIONS_H 1

#include <vector>
#include <string>
#include "annotation.hpp"
#include "read_stats.hpp"

using namespace std;

void write_fusions_to_file(fusions_t& fusions, const string& output_file, const coverage_t& coverage, const assembly_t& assembly, gene_annotation_index_t& gene_annotation_index, exon_annotation_index_t& exon_annotation_index, vector<string> contigs_by_id, const bool print_supporting_reads, const bool print_fusion_sequence, const bool print_peptide_sequence, const bool write_discarded_fusions);

#endif /* _OUTPUT_FUSIONS_H */
