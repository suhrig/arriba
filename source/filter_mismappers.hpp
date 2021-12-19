#ifndef FILTER_MISMAPPER_H
#define FILTER_MISMAPPER_H 1

#include <string>
#include <unordered_map>
#include <vector>
#include "common.hpp"
#include "annotation.hpp"
#include "assembly.hpp"

using namespace std;

typedef unsigned int kmer_as_int_t; // represent kmer as integer
typedef unordered_map< kmer_as_int_t, vector<int> > kmer_index_t; // store coordinates of kmers
typedef vector<kmer_index_t> kmer_indices_t; // one index per contig

kmer_as_int_t kmer_to_int(const string& kmer, const string::size_type position, const char kmer_length);
void make_kmer_index(const fusions_t& fusions, const assembly_t& assembly, int padding, const char kmer_length, kmer_indices_t& kmer_indices);

unsigned int filter_mismappers(fusions_t& fusions, const kmer_indices_t& kmer_indices, const char kmer_length, const assembly_t& assembly, const exon_annotation_index_t& exon_annotation_index, const float max_mismapper_fraction, const int max_mate_gap);

#endif /* FILTER_MISMAPPERS_H */
