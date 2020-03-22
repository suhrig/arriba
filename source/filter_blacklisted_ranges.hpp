#ifndef _FILTER_BLACKLISTED_RANGES_H
#define _FILTER_BLACKLISTED_RANGES_H 1

#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include "common.hpp"

using namespace std;

// parse string representation of a blacklist rule
enum blacklist_item_type_t { BLACKLIST_RANGE, BLACKLIST_POSITION, BLACKLIST_GENE, BLACKLIST_ANY, BLACKLIST_SPLIT_READ_DONOR, BLACKLIST_SPLIT_READ_ACCEPTOR, BLACKLIST_SPLIT_READ_ANY, BLACKLIST_DISCORDANT_MATES, BLACKLIST_READ_THROUGH, BLACKLIST_LOW_SUPPORT, BLACKLIST_FILTER_SPLICED, BLACKLIST_NOT_BOTH_SPLICED };
struct blacklist_item_t {
	blacklist_item_type_t type;
	bool strand_defined;
	strand_t strand;
	contig_t contig;
	position_t start;
	position_t end;
	gene_t gene;
};
bool parse_blacklist_item(const string& text, blacklist_item_t& blacklist_item, const contigs_t& contigs, const unordered_map<string,gene_t>& genes, const bool allow_keyword);

// check if the breakpoint of a fusion match an entry in the blacklist
bool matches_blacklist_item(const blacklist_item_t& blacklist_item, const fusion_t& fusion, const unsigned char which_breakpoint, const int max_mate_gap, const float evalue_cutoff = 0);

// divide the genome into fixed size bins
typedef vector< tuple<contig_t,position_t> > genome_bins_t;
void get_genome_bins_from_range(const contig_t contig, const position_t start, const position_t end, genome_bins_t& genome_bins);

unsigned int filter_blacklisted_ranges(fusions_t& fusions, const string& blacklist_file_path, const contigs_t& contigs, const unordered_map<string,gene_t>& genes, const float evalue_cutoff, const int max_mate_gap);

#endif /* _FILTER_BLACKLISTED_RANGES_H */
