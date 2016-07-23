#ifndef _ANNOTATION_H
#define _ANNOTATION_H 1

#include <unordered_map>
#include <map>
#include <vector>
#include <set>
#include <string>
#include "common.hpp"

using namespace std;

struct annotation_record_t {
	unsigned int id;
	string name;
	string sequence;
	contig_t contig;
	position_t start;
	position_t end;
	strand_t strand;
	int exonic_length; // sum of the length of all exons in a gene
	bool is_dummy;
	// function to sort annotation records by coordinate
	inline bool operator < (const annotation_record_t& x) const {
		if (contig != x.contig) return contig < x.contig;
		return end < x.end;
	}
};
typedef vector<annotation_record_t> annotation_t;
typedef map<position_t,gene_multiset_t> contig_annotation_index_t;
typedef vector<contig_annotation_index_t> annotation_index_t; //TODO make this a tuple<contig,position>

string removeChr(string contig);
string addChr(string contig);

void read_annotation_gtf(const string& filename, annotation_t& gene_annotation, annotation_t& exon_annotation, contigs_t& contigs);

void make_annotation_index(const annotation_t& annotation, annotation_index_t& annotation_index, const contigs_t& contigs);

void gene_multiset_to_set(gene_multiset_t gene_multiset, gene_set_t& gene_set);
gene_set_t gene_multiset_to_set(gene_multiset_t gene_multiset);

void combine_annotations(const gene_set_t& genes1, const gene_set_t& genes2, gene_set_t& combined, bool make_union = true);

void get_annotation_by_coordinate(const contig_t contig, const position_t start, const position_t end, gene_set_t& gene_set, annotation_index_t& annotation_index);

void get_boundaries_of_biggest_gene(gene_set_t& genes, const annotation_t& gene_annotation, position_t& start, position_t& end);

void fetch_gene_sequences_from_fasta(const string& assembly_file_path, fusions_t& fusions, annotation_t& gene_annotation, const vector<string>& contigs_by_id);

void dna_to_reverse_complement(string& dna, string& reverse_complement);

string dna_to_reverse_complement(string& dna);

bool is_breakpoint_spliced(const gene_t gene, const direction_t direction, const contig_t contig, const position_t breakpoint, const annotation_index_t& exon_annotation_index);


#endif /*_ANNOTATION_H*/
