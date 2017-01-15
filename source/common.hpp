#ifndef _COMMON_H
#define _COMMON_H 1

#include <deque>
#include <map>
#include <string>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>
#include "sam.h"

using namespace std;

typedef bool strand_t;
const strand_t FORWARD = true;
const strand_t REVERSE = false;

typedef const string* filter_t;
extern unordered_map<string,filter_t> FILTERS;

typedef unsigned short int contig_t;
typedef unordered_map<string,contig_t> contigs_t;
typedef int position_t;

struct annotation_record_t {
	contig_t contig;
	position_t start;
	position_t end;
	strand_t strand;
	// function to sort annotation records by coordinate
	inline bool operator < (const annotation_record_t& x) const {
		if (contig != x.contig) return contig < x.contig;
		return end < x.end;
	}
	void copy(const annotation_record_t& x) {
		*this = x;
	}
};
template <class T> class annotation_set_t: public set<T> {};
template <class T> class annotation_multiset_t: public multiset<T> {};
template <class T> class annotation_t: public deque<T> {}; // must be a deque for stable pointers
template <class T> class contig_annotation_index_t: public map< position_t, annotation_multiset_t<T> > {};
template <class T> class annotation_index_t: public vector< contig_annotation_index_t<T> > {};

struct gene_annotation_record_t: public annotation_record_t {
	string name;
	string sequence;
	int exonic_length; // sum of the length of all exons in a gene
	bool is_dummy;
	bool is_known;
};
typedef gene_annotation_record_t* gene_t;
typedef annotation_set_t<gene_t> gene_set_t;
typedef annotation_multiset_t<gene_t> gene_multiset_t;
typedef annotation_t<gene_annotation_record_t> gene_annotation_t;
typedef contig_annotation_index_t<gene_t> gene_contig_annotation_index_t;
typedef annotation_index_t<gene_t> gene_annotation_index_t;

typedef unsigned int transcript_t;
struct exon_annotation_record_t: public annotation_record_t {
	gene_t gene;
	transcript_t transcript;
	bool is_transcript_start, is_transcript_end;
};
typedef exon_annotation_record_t* exon_t;
typedef annotation_set_t<exon_t> exon_set_t;
typedef annotation_multiset_t<exon_t> exon_multiset_t;
typedef annotation_t<exon_annotation_record_t> exon_annotation_t;
typedef contig_annotation_index_t<exon_t> exon_contig_annotation_index_t;
typedef annotation_index_t<exon_t> exon_annotation_index_t;

class cigar_t: public vector<uint32_t> {
	public:
		uint32_t operation(unsigned int index) const { return this->at(index) & 15; }; // select lower 4 bits to get the operation of the CIGAR element
		uint32_t op_length(unsigned int index) const { return this->at(index) >> 4; }; // remove lower 4 bits to get the length of the CIGAR element
		cigar_t operator=(const bam1_t* bam_record) { for (int i = 0; i < bam_record->core.n_cigar; ++i) { this->push_back(bam1_cigar(bam_record)[i]); }; return *this; };
};

struct alignment_t {
	bool supplementary;
	bool first_in_pair;
	bool exonic;
	strand_t strand;
	strand_t predicted_strand;
	bool predicted_strand_ambiguous;
	contig_t contig;
	position_t start;
	position_t end;
	cigar_t cigar;
	string sequence;
	gene_set_t genes;
	alignment_t(): supplementary(false), first_in_pair(false), exonic(false), predicted_strand_ambiguous(true) {};
	unsigned int preclipping() { return (cigar.operation(0) & (BAM_CSOFT_CLIP | BAM_CHARD_CLIP)) ? cigar.op_length(0) : 0; };
	unsigned int postclipping() { return (cigar.operation(cigar.size()-1) & (BAM_CSOFT_CLIP | BAM_CHARD_CLIP)) ? cigar.op_length(cigar.size()-1) : 0; };
};
const unsigned int MATE1 = 0;
const unsigned int MATE2 = 1;
const unsigned int SPLIT_READ = 1;
const unsigned int SUPPLEMENTARY = 2;
class mates_t: public vector<alignment_t> {
	public:
		string name;
		filter_t filter; // name of the filter which discarded the reads (NULL means not discarded)
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
	gene_t gene1, gene2;
	contig_t contig1, contig2;
	position_t breakpoint1, breakpoint2;
	direction_t direction1, direction2;
	strand_t predicted_strand1, predicted_strand2;
	bool predicted_strands_ambiguous:1;
	transcript_start_t transcript_start;
	bool transcript_start_ambiguous:1;
	bool exonic1:1, exonic2:1;
	bool overlap_duplicate1:1, overlap_duplicate2:1;
	bool spliced1:1, spliced2:1;
	confidence_t confidence:2;
	unsigned int split_reads1, split_reads2;
	unsigned int discordant_mates;
	position_t closest_genomic_breakpoint1, closest_genomic_breakpoint2;
	position_t anchor_start1, anchor_start2;
	vector<mates_t*> split_read1_list, split_read2_list, discordant_mate_list;
	float evalue; // expected number of fusions with the given properties by random chance
	filter_t filter; // name of the filter which discarded the fusion (NULL means not discarded)
	fusion_t(): split_reads1(0), split_reads2(0), discordant_mates(0), anchor_start1(0), anchor_start2(0), closest_genomic_breakpoint1(-1), closest_genomic_breakpoint2(-1), filter(NULL) {};
	unsigned int supporting_reads() const { return split_reads1 + split_reads2 + discordant_mates; };
	bool breakpoint_overlaps_both_genes(const unsigned int which_breakpoint = 0) const {
		if (which_breakpoint == 1) return breakpoint1 >= gene2->start && breakpoint1 <= gene2->end;
		if (which_breakpoint == 2) return breakpoint2 >= gene1->start && breakpoint2 <= gene1->end;
		return breakpoint_overlaps_both_genes(1) || breakpoint_overlaps_both_genes(2);
	};
	bool is_read_through() const { return contig1 == contig2 && breakpoint2 - breakpoint1 < 400000 && direction1 == DOWNSTREAM && direction2 == UPSTREAM; };
};
typedef unordered_map< tuple<gene_t /*gene1*/, gene_t /*gene2*/, contig_t /*contig1*/, contig_t /*contig2*/, position_t /*breakpoint1*/, position_t /*breakpoint2*/, direction_t /*direction1*/, direction_t /*direction2*/>,fusion_t > fusions_t;

typedef char strandedness_t;
const strandedness_t STRANDEDNESS_NO = 0;
const strandedness_t STRANDEDNESS_YES = 1;
const strandedness_t STRANDEDNESS_REVERSE = 2;

namespace std{
    namespace
    {

	// Code from boost
	// Reciprocal of the golden ratio helps spread entropy
	//     and handles duplicates.
	// See Mike Seymour in magic-numbers-in-boosthash-combine:
	//     http://stackoverflow.com/questions/4948780

	template <class T>
	inline void hash_combine(std::size_t& seed, T const& v)
	{
	    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
	}

	// Recursive template code derived from Matthieu M.
	template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
	struct HashValueImpl
	{
	  static void apply(size_t& seed, Tuple const& tuple)
	  {
	    HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
	    hash_combine(seed, std::get<Index>(tuple));
	  }
	};

	template <class Tuple>
	struct HashValueImpl<Tuple,0>
	{
	  static void apply(size_t& seed, Tuple const& tuple)
	  {
	    hash_combine(seed, std::get<0>(tuple));
	  }
	};
    }

    template <typename ... TT>
    struct hash<std::tuple<TT...>> 
    {
	size_t
	operator()(std::tuple<TT...> const& tt) const
	{					      
	    size_t seed = 0;			     
	    HashValueImpl<std::tuple<TT...> >::apply(seed, tt);    
	    return seed;				 
	}					      

    };
}

#endif /* _COMMON_H */
