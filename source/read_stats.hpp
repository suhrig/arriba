#ifndef _READ_STATS_H
#define _READ_STATS_H 1

#include <string>
#include <vector>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"

using namespace std;

bool estimate_mate_gap_distribution(const chimeric_alignments_t& chimeric_alignments, float& mate_gap_mean, float& mate_gap_stddev, const gene_annotation_index_t& gene_annotation_index, const exon_annotation_index_t& exon_annotation_index);

strandedness_t detect_strandedness(const chimeric_alignments_t& chimeric_alignments, const gene_annotation_index_t& gene_annotation_index);

const unsigned int COVERAGE_RESOLUTION = 20; // at what resolution in bp to calculate the coverage
// for each contig store for every window of <COVERAGE_RESOLUTION> bp whether a read starts/ends here
// this information is needed by the 'no_coverage' filter
// at a later point, this structure may also store coverage information
class coverage_t {
	private:
		vector< vector<bool> > fragment_starts;
		vector< vector<bool> > fragment_ends;
	public:

		// initialize data structure to compute coverage for windows of size <COVERAGE_RESOLUTION>
		coverage_t(const contigs_t& contigs, const assembly_t& assembly) {
			fragment_starts.resize(contigs.size());
			fragment_ends.resize(contigs.size());
			for (assembly_t::const_iterator contig = assembly.begin(); contig != assembly.end(); ++contig) {
				if (!contig->second.empty()) {
					fragment_starts[contig->first].resize((contig->second.size()+COVERAGE_RESOLUTION-1)/COVERAGE_RESOLUTION);
					fragment_ends[contig->first].resize((contig->second.size()+COVERAGE_RESOLUTION-1)/COVERAGE_RESOLUTION);
				}
			}
		}

		// add alignment to coverage
		void add_alignment(const bam1_t* bam_record) {
			if (!(bam_record->core.flag & BAM_FREVERSE) || !(bam_record->core.flag & BAM_FPAIRED))
				fragment_starts[bam_record->core.tid][bam_record->core.pos/COVERAGE_RESOLUTION] = true;
			if ((bam_record->core.flag & BAM_FREVERSE) || !(bam_record->core.flag & BAM_FPAIRED))
				fragment_ends[bam_record->core.tid][(bam_endpos(bam_record)-1)/COVERAGE_RESOLUTION] = true;
		};

		// returns true, if a fragment begins at the given position
		bool fragment_starts_here(const contig_t contig, const position_t start, const position_t end) const {
			if (contig >= fragment_starts.size())
				return false;
			for (unsigned int window = start/COVERAGE_RESOLUTION + 1; window <= end/COVERAGE_RESOLUTION; ++window)
				if (fragment_starts[contig][window])
					return true;
			return false;
		};

		// returns true, if a fragment ends at the given position
		bool fragment_ends_here(const contig_t contig, const position_t start, const position_t end) const {
			if (contig >= fragment_ends.size())
				return false;
			for (unsigned int window = start/COVERAGE_RESOLUTION; window < end/COVERAGE_RESOLUTION; ++window)
				if (fragment_ends[contig][window])
					return true;
			return false;
		};
};

#endif /* _READ_STATS_H */
