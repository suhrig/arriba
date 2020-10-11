#include <string>
#include <unordered_map>
#include "common.hpp"
#include "filter_blacklisted_ranges.hpp"
#include "read_compressed_file.hpp"
#include "annotate_tags.hpp"

using namespace std;

void load_tags(const string& tags_file_path, const contigs_t& contigs, const unordered_map<string,gene_t>& genes, tags_t& tags) {
	autodecompress_file_t tags_file(tags_file_path);
	string line;
	while (tags_file.getline(line)) {
		if (!line.empty() && line[0] != '#') {

			// parse line
			tsv_stream_t tsv(line);
			string range1, range2, tag;
			tsv >> range1 >> range2 >> tag;
			blacklist_item_t item1, item2;
			if (!parse_blacklist_item(range1, item1, contigs, genes, false) ||
			    !parse_blacklist_item(range2, item2, contigs, genes, false))
				continue;

			// replace special characters with underscores
			for (string::size_type pos = 0; pos < tag.size(); ++pos)
				if (tag[pos] < '!' || tag[pos] > '~' || tag[pos] == ',')
					tag[pos] = '_';

			// index tags by coordinate
			genome_bins_t genome_bins;
			get_genome_bins_from_range(item1.contig, item1.start, item1.end, genome_bins);
			get_genome_bins_from_range(item2.contig, item2.start, item2.end, genome_bins);
			for (auto genome_bin = genome_bins.begin(); genome_bin != genome_bins.end(); ++genome_bin)
				tags[*genome_bin].push_back(make_tuple(item1, item2, tag));
		}
	}
}

string annotate_tags(const fusion_t& fusion, const tags_t& tags, const int max_mate_gap) {

	// determine bins around fusion breakpoints; the tags candidates are located in these bins
	genome_bins_t genome_bins;
	get_genome_bins_from_range(fusion.contig1, fusion.breakpoint1, fusion.breakpoint1, genome_bins);
	get_genome_bins_from_range(fusion.contig2, fusion.breakpoint2, fusion.breakpoint2, genome_bins);
	get_genome_bins_from_range(fusion.contig1, fusion.gene1->start, fusion.gene1->end, genome_bins);
	get_genome_bins_from_range(fusion.contig2, fusion.gene2->start, fusion.gene2->end, genome_bins);

	// iterate through tag candidates in bins and find those that truly match
	set<string> matching_tags;
	for (auto genome_bin = genome_bins.begin(); genome_bin != genome_bins.end(); ++genome_bin) {
		auto tag_candidates = tags.find(*genome_bin);
		if (tag_candidates != tags.end()) {
			for (auto tag = tag_candidates->second.begin(); tag != tag_candidates->second.end(); ++tag) {
				// 5' gene of predicted fusion must match gene in 1st column of tags list
				// 3' gene of predicted fusion must match gene in 2nd column of tags list
				const unsigned char gene_5 = (fusion.transcript_start == TRANSCRIPT_START_GENE1) ? 1 : 2;
				const unsigned char gene_3 = (fusion.transcript_start != TRANSCRIPT_START_GENE1) ? 1 : 2;
				if (matches_blacklist_item(get<0>(*tag), fusion, gene_5, max_mate_gap) &&
				    matches_blacklist_item(get<1>(*tag), fusion, gene_3, max_mate_gap)) {
					matching_tags.insert(get<2>(*tag));
				}
			}
		}
	}

	// concatenate tags to comma-separated string
	string result;
	for (auto tag = matching_tags.begin(); tag != matching_tags.end(); ++tag) {
		if (!result.empty())
			result += ",";
		result += *tag;
	}
	if (result.empty())
		return ".";
	else
		return result;
}

