#ifndef RECOVER_KNOWN_FUSIONS_H
#define RECOVER_KNOWN_FUSIONS_H 1

#include <string>
#include <unordered_map>
#include "common.hpp"
#include "annotation.hpp"
#include "read_stats.hpp"

using namespace std;

unsigned int recover_known_fusions(fusions_t& fusions, const string& known_fusions_file_path, const contigs_t& contigs, const unordered_map<string,gene_t>& genes, const coverage_t& coverage, const int max_mate_gap);

#endif /* RECOVER_KNOWN_FUSIONS_H */
