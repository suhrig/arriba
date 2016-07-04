#ifndef _FILTER_PROMISCUOUS_GENES_H
#define _FILTER_PROMISCUOUS_GENES_H 1

#include <string>
#include <unordered_map>
#include <vector>
#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned long int count_mapped_reads(const string& bam_file_path, const vector<bool>& interesting_contigs);

void estimate_expected_fusions(fusions_t& fusions, const annotation_t& gene_annotation, const unsigned long int mapped_reads);

unsigned int filter_promiscuous_genes(fusions_t& fusions, const float evalue_cutoff);

#endif /* _FILTER_PROMISCUOUS_GENES_H */
