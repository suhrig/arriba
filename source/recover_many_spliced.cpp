#include <map>
#include <tuple>
#include "common.hpp"
#include "recover_many_spliced.hpp"

using namespace std;

bool fusion_has_mismappers(fusion_t& fusion) {
	for (auto mates = fusion.split_read1_list.begin(); mates != fusion.split_read1_list.end(); ++mates)
		if ((**mates).filter == FILTERS.at("mismappers"))
			return true;
	for (auto mates = fusion.split_read2_list.begin(); mates != fusion.split_read2_list.end(); ++mates)
		if ((**mates).filter == FILTERS.at("mismappers"))
			return true;
	return false;
}

unsigned int recover_many_spliced(fusions_t& fusions, const unsigned int min_spliced_events) {

	// look for any spliced reads between two genes
	map< tuple<gene_t,gene_t>, unsigned int > spliced_fusions_by_gene_pair;
	map< tuple<gene_t,gene_t>, unsigned int > spliced_fusions_without_mismappers_by_gene_pair;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (!fusion->second.is_read_through() &&
		    (fusion->second.spliced1 || fusion->second.spliced2) &&
		    fusion->second.gene1 != fusion->second.gene2 &&
		    !fusion->second.breakpoint_overlaps_both_genes() &&
		    (fusion->second.filter == NULL ||
		     fusion->second.filter == FILTERS.at("promiscuous_genes") ||
		     fusion->second.filter == FILTERS.at("min_support") ||
		     fusion->second.filter == FILTERS.at("end_to_end") ||
		     fusion->second.filter == FILTERS.at("select_best") ||
		     fusion->second.filter == FILTERS.at("short_anchor") ||
		     fusion->second.filter == FILTERS.at("mismappers") ||
		     fusion->second.filter == FILTERS.at("non_expressed"))) {
			spliced_fusions_by_gene_pair[make_tuple(fusion->second.gene1, fusion->second.gene2)]++;
			// homologous genes have many spliced fusions, but ALL of them are discarded by the 'mismappers' filter
			// we only recover fusions with many spliced events, if at least one was not discarded by the 'mismappers' filter
			if (!fusion_has_mismappers(fusion->second))
				spliced_fusions_without_mismappers_by_gene_pair[make_tuple(fusion->second.gene1, fusion->second.gene2)]++;
		}

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter == NULL) { // fusion has not been filtered, no need to recover
			remaining++;
			continue;
		}

		if (fusion->second.is_read_through() ||
		    fusion->second.gene1 == fusion->second.gene2 ||
		    fusion->second.breakpoint_overlaps_both_genes())
			continue;

		if (fusion->second.filter != FILTERS.at("short_anchor") &&
		    fusion->second.filter != FILTERS.at("end_to_end") &&
		    fusion->second.filter != FILTERS.at("mismappers") &&
		    fusion->second.filter != FILTERS.at("non_expressed"))
			continue; // only undo some speculative filters

		if (spliced_fusions_by_gene_pair[make_tuple(fusion->second.gene1, fusion->second.gene2)] >= min_spliced_events &&
		    spliced_fusions_without_mismappers_by_gene_pair[make_tuple(fusion->second.gene1, fusion->second.gene2)] > 0) {
			fusion->second.filter = NULL;
			remaining++;
		}
	}
	return remaining;
}
