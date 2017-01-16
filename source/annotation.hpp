#ifndef _ANNOTATION_H
#define _ANNOTATION_H 1

#include <unordered_map>
#include <map>
#include <vector>
#include <set>
#include <string>
#include "common.hpp"

using namespace std;

// coordinates that are at most this far away from a splice-site are considered to be at the splice-site
const unsigned int MAX_SPLICE_SITE_DISTANCE = 2;

struct gtf_features_t {
	string gene_name;
	string gene_id;
	string transcript_id;
	string gene_status;
	string status_known;
	string gene_type;
	string type_protein_coding;
	string feature_exon;
	string feature_utr;
	string feature_gene;
};
const string DEFAULT_GTF_FEATURES = "gene_name=gene_name gene_id=gene_id transcript_id=transcript_id gene_status=gene_status status_KNOWN=KNOWN gene_type=gene_type type_protein_coding=protein_coding feature_exon=exon feature_UTR=UTR feature_gene=gene";

bool parse_gtf_features(string gtf_features_string, gtf_features_t& gtf_features);

string removeChr(string contig);
string addChr(string contig);

void read_annotation_gtf(const string& filename, const contigs_t& contigs, const string& gtf_features_string, gene_annotation_t& gene_annotation, exon_annotation_t& exon_annotation, unordered_map<string,gene_t>& gene_names);

template <class T> void make_annotation_index(annotation_t<T>& annotation, annotation_index_t<T*>& annotation_index, const contigs_t& contigs);

template <class T> void annotation_multiset_to_set(annotation_multiset_t<T> annotation_multiset, annotation_set_t<T>& annotation_set);
template <class T> annotation_set_t<T> annotation_multiset_to_set(annotation_multiset_t<T> annotation_multiset);

bool is_breakpoint_spliced(const gene_t gene, const direction_t direction, const contig_t contig, const position_t breakpoint, const exon_annotation_index_t& exon_annotation_index);

template <class T> void combine_annotations(const annotation_set_t<T>& genes1, const annotation_set_t<T>& genes2, annotation_set_t<T>& combined, bool make_union = true);

template <class T> void get_annotation_by_coordinate(const contig_t contig, const position_t start, const position_t end, annotation_set_t<T>& annotation_set, const annotation_index_t<T>& annotation_index);

void annotate_alignments(mates_t& mates, const exon_annotation_index_t& exon_annotation_index);

void get_boundaries_of_biggest_gene(gene_set_t& genes, position_t& start, position_t& end);

void fetch_gene_sequences_from_fasta(const string& assembly_file_path, fusions_t& fusions, const vector<string>& contigs_by_id);

void dna_to_reverse_complement(string& dna, string& reverse_complement);

string dna_to_reverse_complement(string& dna);

int get_spliced_distance(const contig_t contig, const position_t position1, const position_t position2, const direction_t direction1, const direction_t direction2, const gene_t gene, const exon_annotation_index_t& exon_annotation_index);

// get the complementary strand
inline strand_t complement_strand(const strand_t strand) {
	return (strand == FORWARD) ? REVERSE : FORWARD;
}

// only return the complementary strand, if condition is met
inline strand_t complement_strand_if(const strand_t strand, const bool condition) {
	if (condition) {
		return complement_strand(strand);
	} else {
		return strand;
	}
}

// include template functions
#include "annotation.t.hpp"

#endif /*_ANNOTATION_H*/
