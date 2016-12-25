#include "common.hpp"
#include "filter_homopolymer.hpp"

using namespace std;



unsigned int filter_homopolymer(chimeric_alignments_t& chimeric_alignments, const unsigned int homopolymer_length) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {
		if (i->second.filter != NULL)
			continue; // read has already been filtered

		if (i->second.size() == 3) { // these are alignments of a split read

			// get sequences near breakpoint
			string sequence = "";
			if (i->second[SPLIT_READ].strand == FORWARD) {
				if (i->second[SPLIT_READ].preclipping() >= homopolymer_length)
					sequence += i->second[SPLIT_READ].sequence.substr(i->second[SPLIT_READ].preclipping() - homopolymer_length, homopolymer_length) + " ";
				if (i->second[SPLIT_READ].sequence.length() - i->second[SPLIT_READ].preclipping() >= homopolymer_length)
					sequence += i->second[SPLIT_READ].sequence.substr(i->second[SPLIT_READ].preclipping(), homopolymer_length) + " ";
			} else { // i->second[SPLIT_READ].strand == REVERSE
				if (i->second[SPLIT_READ].postclipping() >= homopolymer_length)
					sequence += i->second[SPLIT_READ].sequence.substr(i->second[SPLIT_READ].sequence.length() - i->second[SPLIT_READ].postclipping(), homopolymer_length) + " ";
				if (i->second[SPLIT_READ].sequence.length() - i->second[SPLIT_READ].postclipping() >= homopolymer_length)
					sequence += i->second[SPLIT_READ].sequence.substr(i->second[SPLIT_READ].sequence.length() - i->second[SPLIT_READ].postclipping() - homopolymer_length, homopolymer_length) + " ";
			}

			// check for homopolymers
			unsigned int run = 1;
			for (unsigned int c = 1; c < sequence.length(); c++) {
				if (sequence[c-1] == sequence[c]) {
					run++;
					if (run == homopolymer_length) {
						i->second.filter = FILTERS.at("homopolymer");
						goto next_read;
					}
				} else {
					run = 1;
				}
			}
	
		}

		++remaining;

		next_read: NULL; // NULL is a dummy statement for the goto label
	}

	return remaining;
}

