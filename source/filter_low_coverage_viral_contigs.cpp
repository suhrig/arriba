#include <string>
#include "common.hpp"
#include "read_stats.hpp"
#include "filter_no_coverage.hpp"

using namespace std;

// when a tumor is truly infected by a virus, there is fairly homogeneous coverage of the respective viral contig
// in contrast, viral contigs that attract alignment artifacts have very focal coverage
// => remove viral contigs if there is high coverage, but the coverage is focal
unsigned int filter_low_coverage_viral_contigs(chimeric_alignments_t& chimeric_alignments, const coverage_t& coverage, const contigs_t& contigs, const string& viral_contigs, const float min_covered_fraction) {

	// compute average coverage for each viral contig
	vector<float> average_coverage(contigs.size());
	for (auto contig = contigs.begin(); contig != contigs.end(); ++contig) {
		for (auto window = coverage.coverage[contig->second].begin(); window != coverage.coverage[contig->second].end(); ++window)
			average_coverage[contig->second] += *window;
		average_coverage[contig->second] /= coverage.coverage[contig->second].size();
	}

	// for each viral contig, determine fraction that is covered at least as highly as 0.05 * average coverage
	vector<float> fraction_with_sufficient_coverage(contigs.size());
	for (auto contig = contigs.begin(); contig != contigs.end(); ++contig) {
		for (auto window = coverage.coverage[contig->second].begin(); window != coverage.coverage[contig->second].end(); ++window)
			if (*window > 0.05 * average_coverage[contig->second])
				fraction_with_sufficient_coverage[contig->second]++;
		fraction_with_sufficient_coverage[contig->second] /= coverage.coverage[contig->second].size();
	}

	// convert viral_contigs to vector of booleans for faster lookup
	vector<bool> viral_contigs_bool(contigs.size());
	for (contigs_t::const_iterator contig = contigs.begin(); contig != contigs.end(); ++contig)
		viral_contigs_bool[contig->second] = is_interesting_contig(contig->first, viral_contigs);

	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != FILTER_none)
			continue; // alignment has already been filtered

		// remove alignments mapping to viral contigs with focal coverage
		for (mates_t::iterator mate = chimeric_alignment->second.begin(); mate != chimeric_alignment->second.end(); ++mate) {
			if (viral_contigs_bool[mate->contig] && fraction_with_sufficient_coverage[mate->contig] < min_covered_fraction) {
				chimeric_alignment->second.filter = FILTER_low_coverage_viral_contigs;
				goto next_read;
			}
		}

		remaining++;

		next_read: continue;
	}

	return remaining;
}

