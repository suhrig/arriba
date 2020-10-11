#include <vector>
#include "common.hpp"
#include "read_stats.hpp"
#include "filter_no_coverage.hpp"

using namespace std;

// when a tumor is truly infected by a virus, there is fairly homogeneous coverage of the respective viral contig
// in contrast, viral contigs that attract alignment artifacts have very focal coverage
// => remove viral contigs if there is high coverage, but the coverage is focal
unsigned int filter_low_coverage_viral_contigs(chimeric_alignments_t& chimeric_alignments, const coverage_t& coverage, const vector<bool>& viral_contigs, const float min_covered_fraction) {

	// compute average coverage for each viral contig
	vector<float> average_coverage(viral_contigs.size());
	for (contig_t contig = 0; contig < viral_contigs.size(); ++contig) {
		for (auto window = coverage.coverage[contig].begin(); window != coverage.coverage[contig].end(); ++window)
			average_coverage[contig] += *window;
		average_coverage[contig] /= coverage.coverage[contig].size();
	}

	// for each viral contig, determine fraction that is covered at least as highly as 0.05 * average coverage
	vector<float> fraction_with_sufficient_coverage(viral_contigs.size());
	for (contig_t contig = 0; contig < viral_contigs.size(); ++contig) {
		for (auto window = coverage.coverage[contig].begin(); window != coverage.coverage[contig].end(); ++window)
			if (*window > 0.05 * average_coverage[contig])
				fraction_with_sufficient_coverage[contig]++;
		fraction_with_sufficient_coverage[contig] /= coverage.coverage[contig].size();
	}

	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != FILTER_none)
			continue; // alignment has already been filtered

		// remove alignments mapping to viral contigs with focal coverage
		for (mates_t::iterator mate = chimeric_alignment->second.begin(); mate != chimeric_alignment->second.end(); ++mate) {
			if (viral_contigs[mate->contig] && fraction_with_sufficient_coverage[mate->contig] < min_covered_fraction) {
				chimeric_alignment->second.filter = FILTER_low_coverage_viral_contigs;
				goto next_read;
			}
		}

		remaining++;

		next_read: continue;
	}

	return remaining;
}

