#include <cmath>
#include <map>
#include <tuple>
#include "common.hpp"
#include "recover_both_spliced.hpp"

using namespace std;

direction_t opposite_direction(const direction_t direction) {
	return ((direction == DOWNSTREAM) ? UPSTREAM : DOWNSTREAM);
}

unsigned int recover_both_spliced(fusions_t& fusions, const unsigned int max_fusions_to_recover) {

	// look for any supporting reads between two genes
	map< tuple<gene_t,gene_t,direction_t,direction_t>, vector<fusion_t*> > fusions_by_gene_pair;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter == NULL ||
		    (fusion->second.both_breakpoints_spliced() && fusion->second.filter != FILTERS.at("merge_adjacent")) ||
		    (fusion->second.both_breakpoints_spliced() && fusion->second.filter == FILTERS.at("pcr_fusions")) || // when there is risk of PCR-mediated fusions, only consider spliced events
		    fusion->second.filter == FILTERS.at("intronic") ||
		    fusion->second.filter == FILTERS.at("relative_support") ||
		    fusion->second.filter == FILTERS.at("min_support")) {
			fusions_by_gene_pair[make_tuple(fusion->second.gene1, fusion->second.gene2, (direction_t) fusion->second.direction1, (direction_t) fusion->second.direction2)].push_back(&fusion->second);
		}

	unsigned int remaining = 0;
	const char MODE_COUNTING = 0;
	const char MODE_RECOVER = 1;
	map<unsigned int /*supporting reads*/, unsigned int /*number of recovered fusions*/> recovered_fusions_by_supporting_reads;
	unsigned int min_supporting_reads = 0;
	for (char mode = MODE_COUNTING; mode <= MODE_RECOVER; ++mode) {
		for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

			if (fusion->second.filter == NULL) { // fusion has not been filtered, no need to recover
				if (mode == MODE_RECOVER)
					remaining++;
				continue;
			}

			if (fusion->second.gene1 == fusion->second.gene2)
				continue; // don't recover intragenic events (this would produce too many hits)

			if (fusion->second.filter != NULL &&
			    fusion->second.filter != FILTERS.at("relative_support") &&
			    fusion->second.filter != FILTERS.at("min_support") &&
			    fusion->second.filter != FILTERS.at("pcr_fusions") && fusion->second.discordant_mates <= fusion->second.split_reads1 + fusion->second.split_reads2)
				continue; // we won't recover fusions which were not discarded due to low support

			if (!fusion->second.both_breakpoints_spliced())
				continue; // only recover fusions with two spliced breakpoints

			// count all supporting reads of all fusions between the pair of genes
			unsigned int sum_of_supporting_reads = 0;

			// look for other reads with the same orientation
			auto fusions_of_given_gene_pair = fusions_by_gene_pair.find(make_tuple(fusion->second.gene1, fusion->second.gene2, (direction_t) fusion->second.direction1, (direction_t) fusion->second.direction2));
			if (fusions_of_given_gene_pair != fusions_by_gene_pair.end())
				for (auto another_fusion = fusions_of_given_gene_pair->second.begin(); another_fusion != fusions_of_given_gene_pair->second.end(); ++another_fusion)
					if (fusion->second.filter == FILTERS.at("pcr_fusions")) {
						if ((**another_fusion).both_breakpoints_spliced() && (**another_fusion).discordant_mates <= (**another_fusion).split_reads1 + (**another_fusion).split_reads2)
							sum_of_supporting_reads++; // if there is risk of PCR-mediated fusions, ignore the number of supporting reads and count the event as 1 read
					} else { // the event is probably not PCR-mediated => actually count the number of supporting reads
						sum_of_supporting_reads += max((unsigned int) 1, (**another_fusion).supporting_reads());
					}

			// consider reciprocal translocations
			auto reciprocal_fusions_of_given_gene_pair = fusions_by_gene_pair.find(make_tuple(fusion->second.gene1, fusion->second.gene2, opposite_direction(fusion->second.direction1), opposite_direction(fusion->second.direction2)));
			if (reciprocal_fusions_of_given_gene_pair != fusions_by_gene_pair.end())
				for (auto another_fusion = reciprocal_fusions_of_given_gene_pair->second.begin(); another_fusion != reciprocal_fusions_of_given_gene_pair->second.end(); ++another_fusion)
					if (!(**another_fusion).is_read_through())
						// if the fusion is supported by 2 events of which all breakpoints are spliced, we don't question it
						// if the other fusion is not spliced, then we check if the reciprocal fusions support a common genomic breakpoint
						if ((**another_fusion).both_breakpoints_spliced() ||
						    (((fusion->second.direction1 == DOWNSTREAM) != /*xor*/ (fusion->second.breakpoint1 > (**another_fusion).breakpoint1)) &&
						     ((fusion->second.direction2 == DOWNSTREAM) != /*xor*/ (fusion->second.breakpoint2 > (**another_fusion).breakpoint2))))
							if (fusion->second.filter == FILTERS.at("pcr_fusions")) {
								if ((**another_fusion).both_breakpoints_spliced() && (**another_fusion).discordant_mates <= (**another_fusion).split_reads1 + (**another_fusion).split_reads2)
									sum_of_supporting_reads++; // if there is risk of PCR-mediated fusions, ignore the number of supporting reads and count the event as 1 read
							} else { // the event is probably not PCR-mediated => actually count the number of supporting reads
								sum_of_supporting_reads += max((unsigned int) 1, (**another_fusion).supporting_reads());
							}

			if (sum_of_supporting_reads >= 2) { // require at least two reads or else the false positive rate sky-rockets
				if (mode == MODE_RECOVER) { // we are in recover mode => actually recover the fusion by clearing the filters
					if (fusion->second.supporting_reads() >= min_supporting_reads && !fusion->second.is_read_through() && !fusion->second.breakpoint_overlaps_both_genes()) {
						fusion->second.filter = NULL;
						remaining++;
					}
				} else { // mode == MODE_COUNTING
					// We are in dry-run mode, where we only count how many fusions would be recovered.
					// If too many would be recovered, we increase the minimum required supporting reads,
					// because spliced transcripts might simply be an artifact of the RNA extraction protocol.
					unsigned int stratum = pow(2, ceil(log2(max((unsigned int) 1, fusion->second.supporting_reads()))));
					recovered_fusions_by_supporting_reads[stratum]++;
				}
			}
		}

		if (mode == MODE_COUNTING) {
			// The map "recovered_fusions_by_supporting_reads" now holds for each number of supporting reads
			// the number of fusions that would be recovered (e.g., X fusions with Y supporting reads would be recovered).
			// In some samples we observe an extraordinary number of recovered fusions with less than 3 supporting reads.
			// We don't want to recover those, since most of them are probably false positives.
			// Instead, we increase the minimum number of required supporting reads, e.g., if we observe a gazillion
			// of recovered fusions with 2 supporting reads or less, we only recover ones with more than 2, unless
			// we also see a gazillion of those, in which case we only recover ones with more than 3, and so on.
			// => Find the stratum of fusions, where a moderate number of fusions will be recovered.
			for (auto stratum = recovered_fusions_by_supporting_reads.rbegin(); stratum != recovered_fusions_by_supporting_reads.rend(); ++stratum) {
				if (stratum->second >= max_fusions_to_recover) {
					min_supporting_reads = stratum->first + 1;
					break;
				}
			}
		}
	}
	return remaining;
}
