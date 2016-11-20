#ifndef _FILTER_BLACKLISTED_RANGES_H
#define _FILTER_BLACKLISTED_RANGES_H 1

#include <string>
#include <unordered_map>
#include "common.hpp"

using namespace std;

unsigned int filter_blacklisted_ranges(fusions_t& fusions, const string& blacklist_file_path, const contigs_t& contigs, const unordered_map<string,gene_t>& genes);

#endif /* _FILTER_BLACKLISTED_RANGES_H */
