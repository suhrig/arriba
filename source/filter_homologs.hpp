#ifndef FILTER_HOMOLOGS_H
#define FILTER_HOMOLOGS_H 1

#include "common.hpp"
#include "assembly.hpp"
#include "filter_mismappers.hpp"

using namespace std;

unsigned int filter_homologs(fusions_t& fusions, const kmer_indices_t& kmer_indices, const char kmer_length, const assembly_t& assembly, const float max_identity_fraction);

#endif /* FILTER_HOMOLOGS_H */
