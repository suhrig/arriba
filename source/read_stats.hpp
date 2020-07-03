#ifndef _READ_STATS_H
#define _READ_STATS_H 1

#include <vector>
#include "common.hpp"
#include "annotation.hpp"

using namespace std;

bool estimate_fragment_length(const chimeric_alignments_t& chimeric_alignments, float& mate_gap_mean, float& mate_gap_stddev, float& read_length_mean, const gene_annotation_index_t& gene_annotation_index, const exon_annotation_index_t& exon_annotation_index);

strandedness_t detect_strandedness(const chimeric_alignments_t& chimeric_alignments, const gene_annotation_index_t& gene_annotation_index, const exon_annotation_index_t& exon_annotation_index);

const int COVERAGE_RESOLUTION = 20; // at what resolution in bp to calculate the coverage
// for each contig store for every window of <COVERAGE_RESOLUTION> bp whether a read starts/ends here
// this information is needed by the 'no_coverage' filter
class coverage_t {
	public:
		vector< vector<bool> > fragment_starts; // for each window, store if a fragment starts here
		vector< vector<bool> > fragment_ends; // for each window, store if a fragment ends here
		vector< vector<unsigned short int> > coverage; // for each window, store the coverage
		coverage_t(const contigs_t& contigs, const assembly_t& assembly);
		void add_fragment(bam1_t* mate1, bam1_t* mate2, bool is_chimeric);
		bool fragment_starts_here(const contig_t contig, const position_t start, const position_t end) const;
		bool fragment_ends_here(const contig_t contig, const position_t start, const position_t end) const;
		int get_coverage(const contig_t contig, const position_t position, const direction_t direction) const;
};

#endif /* _READ_STATS_H */
