#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "common.hpp"
#include "annotation.hpp"
#include "filter_blacklisted_ranges.hpp"
#include "read_compressed_file.hpp"
#include "read_stats.hpp"
#include "recover_known_fusions.hpp"

using namespace std;

unsigned int recover_known_fusions(fusions_t& fusions, const string& known_fusions_file_path, const contigs_t& contigs, const unordered_map<string,gene_t>& genes, const coverage_t& coverage, const int max_mate_gap) {

	// the known fusions file has the same format as the blacklist file => we can use the same code
	// known fusions are indexed by coordinate using a hash map for efficient lookup
	unordered_map< tuple<contig_t,position_t>, vector< pair<blacklist_item_t,blacklist_item_t> > > known_fusions_by_coordinate;

	// load known fusions from file
	autodecompress_file_t known_fusions_file(known_fusions_file_path);
	string line;
	while (known_fusions_file.getline(line)) {
		if (!line.empty() && line[0] != '#') {
			tsv_stream_t tsv(line);
			string range1, range2;
			tsv >> range1 >> range2;
			blacklist_item_t item1, item2;
			if (!parse_blacklist_item(range1, item1, contigs, genes, false) ||
			    !parse_blacklist_item(range2, item2, contigs, genes, false))
				continue;
			genome_bins_t genome_bins;
			get_genome_bins_from_range(item1.contig, item1.start, item1.end, genome_bins);
			get_genome_bins_from_range(item2.contig, item2.start, item2.end, genome_bins);
			for (auto genome_bin = genome_bins.begin(); genome_bin != genome_bins.end(); ++genome_bin)
				known_fusions_by_coordinate[*genome_bin].push_back(make_pair(item1, item2));
		}
	}

	// look for known fusions with low support which were filtered
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter == FILTER_none)
			continue; // fusion has not been filtered, no need to recover

		if (fusion->second.gene1 == fusion->second.gene2)
			continue; // don't recover intragenic events

		if (fusion->second.filter != FILTER_relative_support && fusion->second.filter != FILTER_min_support)
			continue; // we won't recover fusions which were not discarded due to low support

		// check if fusion is in list of known fusions
		genome_bins_t genome_bins;
		get_genome_bins_from_range(fusion->second.contig1, fusion->second.breakpoint1, fusion->second.breakpoint1, genome_bins);
		get_genome_bins_from_range(fusion->second.contig2, fusion->second.breakpoint2, fusion->second.breakpoint2, genome_bins);
		get_genome_bins_from_range(fusion->second.contig1, fusion->second.gene1->start, fusion->second.gene1->end, genome_bins);
		get_genome_bins_from_range(fusion->second.contig2, fusion->second.gene2->start, fusion->second.gene2->end, genome_bins);
		for (auto genome_bin = genome_bins.begin(); genome_bin != genome_bins.end(); ++genome_bin) {
			auto known_fusions = known_fusions_by_coordinate.find(*genome_bin);
			if (known_fusions != known_fusions_by_coordinate.end()) {
				for (auto known_fusion = known_fusions->second.begin(); known_fusion != known_fusions->second.end(); ++known_fusion) {
					if (matches_blacklist_item(known_fusion->first,  fusion->second, 1, max_mate_gap) &&
					    matches_blacklist_item(known_fusion->second, fusion->second, 2, max_mate_gap) ||
					    matches_blacklist_item(known_fusion->first,  fusion->second, 2, max_mate_gap) &&
					    matches_blacklist_item(known_fusion->second, fusion->second, 1, max_mate_gap)) {

						if (known_fusion->first.type == BLACKLIST_POSITION && known_fusion->second.type == BLACKLIST_POSITION || // when the whitelist specifies two exact breakpoints, the event is always rescued
						    fusion->second.supporting_reads() >= 2 || // otherwise, we require at least two reads, or else there will be too many false positives
						    fusion->second.both_breakpoints_spliced() && // unless the breakpoints are at splice-sites
						    coverage.get_coverage(fusion->second.contig1, fusion->second.breakpoint1, (fusion->second.direction1 == UPSTREAM) ? DOWNSTREAM : UPSTREAM) +
						    coverage.get_coverage(fusion->second.contig2, fusion->second.breakpoint2, (fusion->second.direction2 == UPSTREAM) ? DOWNSTREAM : UPSTREAM) < 200 &&
						    (fusion->second.contig1 != fusion->second.contig2 || abs(fusion->second.breakpoint2 - fusion->second.breakpoint1) > 1000000))
								fusion->second.filter = FILTER_none;

					}
				}
			}
		}
	}

	// count remaining fusions
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter == FILTER_none)
			++remaining;
	return remaining;
}
