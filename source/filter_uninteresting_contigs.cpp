#include <vector>
#include "common.hpp"
#include "filter_uninteresting_contigs.hpp"

using namespace std;

unsigned int filter_uninteresting_contigs(chimeric_alignments_t& chimeric_alignments, vector<bool>& interesting_contigs) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {
		if (!i->second.filters.empty())
			continue; // the read has already been filtered

		for (mates_t::iterator j = i->second.begin(); ; ++j) {
			if (j == i->second.end()) {
				++remaining;
				break;
			} else if (!interesting_contigs[j->contig]) {
				i->second.filters.insert(FILTERS.at("uninteresting_contigs"));
				break;
			}
		}
	}
	return remaining;
}

