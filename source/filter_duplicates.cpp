#include <unordered_map>
#include <tuple>
#include "common.hpp"
#include "filter_duplicates.hpp"

using namespace std;

unsigned int filter_duplicates(chimeric_alignments_t& chimeric_alignments, const bool external_duplicate_marking) {
	unsigned int remaining = 0;
	unordered_map< tuple<contig_t,contig_t,position_t,position_t> , unsigned int> duplicate_count;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {
		if (chimeric_alignment->second.filter != FILTER_none)
			continue; // read has already been filtered

		if (external_duplicate_marking) {

			// rely on duplicate marking by a preceding program (e.g., when UMIs are used)
			if (chimeric_alignment->second.duplicate)
				chimeric_alignment->second.filter = FILTER_duplicates;
			else
				remaining++;

		} else { // perform our own duplicate marking

			// get start coordinates of reads
			position_t position1 = static_cast<position_t>(
				(chimeric_alignment->second[MATE1].strand == FORWARD) ?
				chimeric_alignment->second[MATE1].start - chimeric_alignment->second[MATE1].preclipping() :
				chimeric_alignment->second[MATE1].end   + chimeric_alignment->second[MATE1].postclipping()
			);
			unsigned int mate2 = (chimeric_alignment->second.size() == 2) ? MATE2 : SUPPLEMENTARY;
			position_t position2 = static_cast<position_t>(
				(chimeric_alignment->second[mate2].strand == FORWARD) ?
				chimeric_alignment->second[mate2].start - chimeric_alignment->second[mate2].preclipping() :
				chimeric_alignment->second[mate2].end   + chimeric_alignment->second[mate2].postclipping()
			);
			contig_t contig1 = chimeric_alignment->second[MATE1].contig;
			contig_t contig2 = chimeric_alignment->second[mate2].contig;

			// always put the mate with the lower coordinate in first position
			// or else we might not recognize the duplicate
			if (position1 > position2) {
				swap(position1, position2);
				swap(contig1, contig2);
			}

			if (duplicate_count[make_tuple(contig1, contig2, position1, position2)]++ > 0)
				chimeric_alignment->second.filter = FILTER_duplicates;
			else
				++remaining;
		}
	}

	return remaining;
}

