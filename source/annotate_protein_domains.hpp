#ifndef _ANNOTATE_PROTEIN_DOMAINS_H
#define _ANNOTATE_PROTEIN_DOMAINS_H 1

#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_map>
#include "common.hpp"

using namespace std;

struct protein_domain_annotation_record_t: public annotation_record_t {
	gene_t gene;
	string name;
};
typedef annotation_t<protein_domain_annotation_record_t> protein_domain_annotation_t;
typedef protein_domain_annotation_record_t* protein_domain_t;
typedef contig_annotation_index_t<protein_domain_t> protein_domain_contig_annotation_index_t;
typedef annotation_index_t<protein_domain_t> protein_domain_annotation_index_t;

void load_protein_domains(const string& filename, const contigs_t& contigs, const gene_annotation_t& gene_annotation, const unordered_map<string,gene_t>& gene_names, protein_domain_annotation_t& protein_domain_annotation, protein_domain_annotation_index_t& protein_domain_annotation_index);

string annotate_retained_protein_domains(const contig_t contig, const position_t breakpoint, const strand_t predicted_strand, const bool predicted_strand_ambiguous, const gene_t gene, const direction_t direction, const protein_domain_annotation_index_t& protein_domain_annotation_index);

char dna_to_protein(const string& triplet);

int get_reading_frame(const vector<position_t>& transcribed_bases, const int from, const int to, const transcript_t transcript, const gene_t gene, const assembly_t& assembly, exon_t& exon_with_start_codon);

string get_fusion_peptide_sequence(const string& transcript_sequence, const vector<position_t>& positions, const gene_t gene_5, const gene_t gene_3, const transcript_t transcript_5, const transcript_t transcript_3, const strand_t predicted_strand_3, const exon_annotation_index_t& exon_annotation_index, const assembly_t& assembly);

string is_in_frame(const string& fusion_peptide_sequence);

#endif /* _ANNOTATE_PROTEIN_DOMAINS_H */
