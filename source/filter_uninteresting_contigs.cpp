#include "common.hpp"
#include "filter_uninteresting_contigs.hpp"

using namespace std;

unsigned int filter_uninteresting_contigs(chimeric_alignments_t& chimeric_alignments, const contigs_t& contigs, const contigs_t& interesting_contigs) {

	// convert interesting_contigs to vector of booleans for faster lookup
	vector<bool> interesting_contigs_bool(contigs.size());
	for (contigs_t::const_iterator contig = contigs.begin(); contig != contigs.end(); ++contig)
		interesting_contigs_bool[contig->second] = (interesting_contigs.find(contig->first) != interesting_contigs.end());

	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != FILTER_none)
			continue; // the read has already been filtered

		// all mates must be on an interesting contig
		for (mates_t::iterator mate = chimeric_alignment->second.begin(); ; ++mate) {
			if (mate == chimeric_alignment->second.end()) {
				++remaining;
				break;
			} else if (!interesting_contigs_bool[mate->contig]) {
				chimeric_alignment->second.filter = FILTER_uninteresting_contigs;
				break;
			}
		}
	}
	return remaining;
}

