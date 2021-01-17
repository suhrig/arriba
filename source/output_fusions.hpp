#ifndef _OUTPUT_FUSIONS_H
#define _OUTPUT_FUSIONS_H 1

#include <string>
#include <vector>
#include "annotation.hpp"
#include "annotate_tags.hpp"
#include "annotate_protein_domains.hpp"
#include "read_stats.hpp"

using namespace std;

void write_fusions_to_file(fusions_t& fusions, const string& output_file, const coverage_t& coverage, const assembly_t& assembly, gene_annotation_index_t& gene_annotation_index, exon_annotation_index_t& exon_annotation_index, vector<string> original_contig_names, const tags_t& tags, const protein_domain_annotation_index_t& protein_domain_annotation_index, const int max_mate_gap, const unsigned max_itd_length, const bool print_extra_info, const bool fill_sequence_gaps, const bool write_discarded_fusions);

#endif /* _OUTPUT_FUSIONS_H */
