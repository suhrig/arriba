#ifndef _COMMON_H
#define _COMMON_H 1

#include <list>
#include <map>
#include <string>
#include <set>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include "sam.h"

using namespace std;

typedef bool strand_t;
const strand_t FORWARD = true;
const strand_t REVERSE = false;

typedef const string* filter_t;
extern unordered_map<string,filter_t> FILTERS;

typedef short int contig_t;
typedef unordered_map<string,contig_t> contigs_t;
typedef int position_t;

typedef unordered_map<contig_t,string> assembly_t;

struct annotation_record_t {
	contig_t contig;
	position_t start;
	position_t end;
	strand_t strand;
	// function to sort annotation records by coordinate
	inline bool operator < (const annotation_record_t& x) const {
		if (contig != x.contig) return contig < x.contig;
		if (end != x.end) return end < x.end;
		return start < x.start;
	}
	void copy(const annotation_record_t& x) { *this = x; }
	inline unsigned int length() const { return this->end - this->start; }
};
template <class T> class annotation_set_t: public vector<T> {
	public:
		typename annotation_set_t<T>::iterator insert(const T& value) {
			typename annotation_set_t<T>::iterator existing_element = lower_bound(this->begin(), this->end(), value);
			if (existing_element == this->end() || *existing_element != value)
				return this->insert(upper_bound(this->begin(), this->end(), value), value);
			else
				return existing_element;
		};
		void insert(typename annotation_set_t<T>::const_iterator first, typename annotation_set_t<T>::const_iterator last) {
			this->reserve(this->size() + distance(first, last));
			for (auto annotation_record = first; annotation_record != last; ++annotation_record)
				this->insert(*annotation_record);
		};
		using vector<T>::insert;
};
template <class T> class annotation_t: public list<T> {};
template <class T> class contig_annotation_index_t: public map< position_t, annotation_set_t<T> > {};
template <class T> class annotation_index_t: public vector< contig_annotation_index_t<T> > {};

struct gene_annotation_record_t: public annotation_record_t {
	unsigned int id;
	string name;
	int exonic_length; // sum of the length of all exons in a gene
	bool is_dummy;
	bool is_protein_coding;
};
typedef gene_annotation_record_t* gene_t;
typedef annotation_set_t<gene_t> gene_set_t;
typedef annotation_t<gene_annotation_record_t> gene_annotation_t;
typedef contig_annotation_index_t<gene_t> gene_contig_annotation_index_t;
typedef annotation_index_t<gene_t> gene_annotation_index_t;

struct transcript_annotation_record_t {
	unsigned int id;
	position_t start;
	position_t end;
};
typedef annotation_t<transcript_annotation_record_t> transcript_annotation_t;
typedef transcript_annotation_record_t* transcript_t;

struct exon_annotation_record_t: public annotation_record_t {
	gene_t gene;
	transcript_t transcript;
	exon_annotation_record_t* previous_exon, * next_exon;
	position_t coding_region_start, coding_region_end;
};
typedef exon_annotation_record_t* exon_t;
typedef annotation_set_t<exon_t> exon_set_t;
typedef annotation_t<exon_annotation_record_t> exon_annotation_t;
typedef contig_annotation_index_t<exon_t> exon_contig_annotation_index_t;
typedef annotation_index_t<exon_t> exon_annotation_index_t;

class cigar_t: public vector<uint32_t> {
	public:
		uint32_t operation(unsigned int index) const { return this->at(index) & 15; }; // select lower 4 bits to get the operation of the CIGAR element
		uint32_t op_length(unsigned int index) const { return this->at(index) >> 4; }; // remove lower 4 bits to get the length of the CIGAR element
};

struct alignment_t {
	bool supplementary;
	bool first_in_pair;
	bool exonic;
	strand_t strand; // strand which the read aligns to
	strand_t predicted_strand; // strand which is predicted to be transcribed
	bool predicted_strand_ambiguous; // true, if transcribed strand cannot be predicted reliably
	contig_t contig;
	position_t start;
	position_t end;
	cigar_t cigar;
	string sequence;
	gene_set_t genes;
	alignment_t(): supplementary(false), first_in_pair(false), exonic(false), predicted_strand_ambiguous(true) {};
	unsigned int preclipping() const { return (cigar.operation(0) == BAM_CSOFT_CLIP || cigar.operation(0) == BAM_CHARD_CLIP) ? cigar.op_length(0) : 0; };
	unsigned int postclipping() const { return (cigar.operation(cigar.size()-1) == BAM_CSOFT_CLIP || cigar.operation(cigar.size()-1) == BAM_CHARD_CLIP) ? cigar.op_length(cigar.size()-1) : 0; };
};
const unsigned int MATE1 = 0;
const unsigned int MATE2 = 1;
const unsigned int SPLIT_READ = 1;
const unsigned int SUPPLEMENTARY = 2;
class mates_t: public vector<alignment_t> {
	public:
		filter_t filter; // name of the filter which discarded the reads (NULL means not discarded)
		bool single_end;
		mates_t(): filter(NULL) {};
};
typedef unordered_map<string,mates_t> chimeric_alignments_t;

typedef unsigned char confidence_t;
const confidence_t CONFIDENCE_LOW = 0;
const confidence_t CONFIDENCE_MEDIUM = 1;
const confidence_t CONFIDENCE_HIGH = 2;

typedef bool direction_t;
const direction_t UPSTREAM = true;
const direction_t DOWNSTREAM = false;

typedef bool transcript_start_t;
const transcript_start_t TRANSCRIPT_START_GENE1 = true;
const transcript_start_t TRANSCRIPT_START_GENE2 = false;

struct fusion_t {
	direction_t direction1:1, direction2:1;
	strand_t predicted_strand1:1, predicted_strand2:1;
	bool predicted_strands_ambiguous:1;
	transcript_start_t transcript_start:1;
	bool transcript_start_ambiguous:1;
	bool exonic1:1, exonic2:1;
	bool spliced1:1, spliced2:1;
	confidence_t confidence:2;
	contig_t contig1, contig2;
	short unsigned int split_reads1, split_reads2, discordant_mates;
	float evalue; // expected number of fusions with the given properties by random chance
	position_t breakpoint1, breakpoint2;
	position_t anchor_start1, anchor_start2;
	position_t closest_genomic_breakpoint1, closest_genomic_breakpoint2;
	gene_t gene1, gene2;
	vector<chimeric_alignments_t::iterator> split_read1_list, split_read2_list, discordant_mate_list;
	filter_t filter; // name of the filter which discarded the fusion (NULL means not discarded)
	fusion_t(): exonic1(false), exonic2(false), split_reads1(0), split_reads2(0), discordant_mates(0), anchor_start1(0), anchor_start2(0), closest_genomic_breakpoint1(-1), closest_genomic_breakpoint2(-1), filter(NULL) {};
	unsigned int supporting_reads() const { return split_reads1 + split_reads2 + discordant_mates; };
	bool breakpoint_overlaps_both_genes(const unsigned int which_breakpoint = 0) const {
		if (which_breakpoint == 1) return breakpoint1 >= gene2->start && breakpoint1 <= gene2->end;
		if (which_breakpoint == 2) return breakpoint2 >= gene1->start && breakpoint2 <= gene1->end;
		return breakpoint_overlaps_both_genes(1) || breakpoint_overlaps_both_genes(2);
	};
	bool is_read_through() const { return contig1 == contig2 && breakpoint2 - breakpoint1 < 400000 && direction1 == DOWNSTREAM && direction2 == UPSTREAM; };
	bool is_intragenic() const {
		return gene1 == gene2 ||
		       breakpoint1 >= gene2->start - 10000 && breakpoint1 <= gene2->end + 10000 &&
		       breakpoint2 >= gene1->start - 10000 && breakpoint2 <= gene1->end + 10000;
	};
	bool both_breakpoints_spliced() {
		return spliced1 && spliced2 &&
		       (gene1->strand == gene2->strand && direction1 != direction2 ||
		        gene1->strand != gene2->strand && direction1 == direction2);
	};
};
typedef unordered_map< tuple<unsigned int /*gene1 id*/, unsigned int /*gene2 id*/, contig_t /*contig1*/, contig_t /*contig2*/, position_t /*breakpoint1*/, position_t /*breakpoint2*/, direction_t /*direction1*/, direction_t /*direction2*/>,fusion_t > fusions_t;

typedef char strandedness_t;
const strandedness_t STRANDEDNESS_NO = 0;
const strandedness_t STRANDEDNESS_YES = 1;
const strandedness_t STRANDEDNESS_REVERSE = 2;
const strandedness_t STRANDEDNESS_AUTO = 3;

// implement hash() function for tuples so they can be used as keys in unordered_maps
namespace std {

#define TUPLE_TYPES std::tuple<TT...>
#define TUPLE_ELEMENT_TYPE std::tuple_element< element,std::tuple<TT...> >::type

        template <typename ... TT> struct hash< TUPLE_TYPES > {

		size_t operator()(const TUPLE_TYPES& tuple, std::integral_constant<int, std::tuple_size< TUPLE_TYPES >::value>) const {
			return 0;
		}
    
		template<int element = 0> size_t operator()(const TUPLE_TYPES& tuple, std::integral_constant<int,element> = std::integral_constant<int,0>()) const {
			return std::hash< typename TUPLE_ELEMENT_TYPE >()(std::get<element>(tuple)) ^ operator()(tuple, std::integral_constant<int,element+1>()) <<4;
		}
	};

#undef TUPLE_TYPES
#undef TUPLE_ELEMENT_TYPE

}

#endif /* _COMMON_H */
