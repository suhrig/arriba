#include "common.hpp"
#include "filter_both_intronic.hpp"

using namespace std;

bool list_contains_exonic_reads(const vector<chimeric_alignments_t::iterator>& read_list) {
	for (auto chimeric_alignments = read_list.begin(); chimeric_alignments != read_list.end(); ++chimeric_alignments)
		if ((**chimeric_alignments).second.filter == FILTER_none)
			for (mates_t::iterator mate = (**chimeric_alignments).second.begin(); mate != (**chimeric_alignments).second.end(); ++mate)
				if (mate->exonic)
					return true;
	return false;
}

unsigned int filter_both_intronic(fusions_t& fusions) {
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter != FILTER_none)
			continue; // read has already been filtered

		if (!list_contains_exonic_reads(fusion->second.split_read1_list) &&
		    !list_contains_exonic_reads(fusion->second.split_read2_list) &&
		    !list_contains_exonic_reads(fusion->second.discordant_mate_list)) {
			fusion->second.filter = FILTER_intronic;
		} else {
			++remaining;
		}
	}

	return remaining;
}

