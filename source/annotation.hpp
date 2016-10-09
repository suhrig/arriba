#ifndef _ANNOTATION_H
#define _ANNOTATION_H 1

#include <unordered_map>
#include <map>
#include <vector>
#include <set>
#include <string>
#include "common.hpp"

using namespace std;

string removeChr(string contig);
string addChr(string contig);

void read_annotation_gtf(const string& filename, const contigs_t& contigs, gene_annotation_t& gene_annotation, exon_annotation_t& exon_annotation, unordered_map<string,gene_t>& gene_names);

template <class T> void make_annotation_index(annotation_t<T>& annotation, annotation_index_t<T*>& annotation_index, const contigs_t& contigs);

template <class T> void annotation_multiset_to_set(annotation_multiset_t<T> annotation_multiset, annotation_set_t<T>& annotation_set);
template <class T> annotation_set_t<T> annotation_multiset_to_set(annotation_multiset_t<T> annotation_multiset);

bool is_breakpoint_spliced(const gene_t gene, const direction_t direction, const contig_t contig, const position_t breakpoint, const exon_annotation_index_t& exon_annotation_index);

template <class T> void combine_annotations(const annotation_set_t<T>& genes1, const annotation_set_t<T>& genes2, annotation_set_t<T>& combined, bool make_union = true);

template <class T> void get_annotation_by_coordinate(const contig_t contig, const position_t start, const position_t end, annotation_set_t<T>& annotation_set, const annotation_index_t<T>& annotation_index);

void get_annotation_by_alignment(const alignment_t& alignment, gene_set_t& gene_set, const exon_annotation_index_t& exon_annotation_index);

void get_boundaries_of_biggest_gene(gene_set_t& genes, position_t& start, position_t& end);

void fetch_gene_sequences_from_fasta(const string& assembly_file_path, fusions_t& fusions, const vector<string>& contigs_by_id);

void dna_to_reverse_complement(string& dna, string& reverse_complement);

string dna_to_reverse_complement(string& dna);

// include template functions
#include "annotation.t.hpp"

#endif /*_ANNOTATION_H*/
