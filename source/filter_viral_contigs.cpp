#include <string>
#include <vector>
#include "common.hpp"
#include "filter_viral_contigs.hpp"

using namespace std;

unsigned int filter_viral_contigs(chimeric_alignments_t& chimeric_alignments, const contigs_t& contigs, const string& viral_contigs) {

	// convert viral_contigs to vector of booleans for faster lookup
	vector<bool> viral_contigs_bool(contigs.size());
	for (contigs_t::const_iterator contig = contigs.begin(); contig != contigs.end(); ++contig)
		viral_contigs_bool[contig->second] = is_interesting_contig(contig->first, viral_contigs);

	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != FILTER_none)
			continue; // the read has already been filtered

		// at least one mate must map to host genome
		for (mates_t::iterator mate = chimeric_alignment->second.begin(); mate != chimeric_alignment->second.end(); ++mate) {
			if (!viral_contigs_bool[mate->contig]) {
				remaining++;
				goto next_read;
			}
		}

		chimeric_alignment->second.filter = FILTER_viral_contigs;

		next_read: continue;
	}
	return remaining;
}

