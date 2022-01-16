#ifndef ANNOTATE_TAGS_H
#define ANNOTATE_TAGS_H 1

#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include "common.hpp"
#include "filter_blacklisted_ranges.hpp"

using namespace std;

typedef unordered_map< genome_bin_t, vector< tuple<blacklist_item_t,blacklist_item_t,string> > > tags_t;

void load_tags(const string& tags_file_path, const contigs_t& contigs, const unordered_map<string,gene_t>& genes, tags_t& tags);

string annotate_tags(const fusion_t& fusion, const tags_t& tags, const int max_mate_gap);

#endif /* ANNOTATE_TAGS_H */
