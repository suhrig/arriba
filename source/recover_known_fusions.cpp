#include <iostream>
#include <fstream>
#include <sstream>
#include <tuple>
#include <string>
#include <set>
#include <unordered_map>
#include "common.hpp"
#include "annotation.hpp"
#include "read_compressed_file.hpp"
#include "read_stats.hpp"
#include "recover_known_fusions.hpp"

using namespace std;

unsigned int recover_known_fusions(fusions_t& fusions, const string& known_fusions_file_path, const unordered_map<string,gene_t>& genes, const coverage_t& coverage) {

	// load known fusions from file
	stringstream known_fusions_file;
	autodecompress_file(known_fusions_file_path, known_fusions_file);
	set< tuple<gene_t,gene_t> > known_fusions;
	string line;
	while (getline(known_fusions_file, line)) {
		if (!line.empty() && line[0] != '#') {
			istringstream iss(line);
			string gene1, gene2;
			iss >> gene1 >> gene2;
			if (genes.find(gene1) != genes.end()) {
				if (genes.find(gene2) != genes.end()) {
					known_fusions.insert(make_tuple(genes.at(gene1), genes.at(gene2)));
					known_fusions.insert(make_tuple(genes.at(gene2), genes.at(gene1)));
				} else {
					cerr << "WARNING: unknown gene in known fusions list: " << gene2 << endl;
				}
			} else {
				cerr << "WARNING: unknown gene in known fusions list: " << gene1 << endl;
			}
		}
	}

	// look for known fusions with low support which were filtered
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter == FILTER_none || // fusion has not been filtered, no need to recover
		    fusion->second.gene1 == fusion->second.gene2) // don't recover intragenic events
			continue;

		if (fusion->second.filter != FILTER_none && // fusion has been filtered
		    fusion->second.filter != FILTER_relative_support && fusion->second.filter != FILTER_min_support) // reason is not low support
			continue; // we won't recover fusions which were not discarded due to low support

		if (fusion->second.supporting_reads() >= 2 || // we still require at least two reads, otherwise there will be too many false positives
		    fusion->second.both_breakpoints_spliced() && // we accept less than two reads, when the breakpoints are at splice-sites
		    coverage.get_coverage(fusion->second.contig1, fusion->second.breakpoint1, (fusion->second.direction1 == UPSTREAM) ? DOWNSTREAM : UPSTREAM) < 100 &&
		    coverage.get_coverage(fusion->second.contig2, fusion->second.breakpoint2, (fusion->second.direction2 == UPSTREAM) ? DOWNSTREAM : UPSTREAM) < 100 &&
		    (fusion->second.contig1 != fusion->second.contig2 || abs(fusion->second.breakpoint2 - fusion->second.breakpoint1) > 1000000))
			if (known_fusions.find(make_tuple(fusion->second.gene1, fusion->second.gene2)) != known_fusions.end()) // fusion is known
				fusion->second.filter = FILTER_none;
	}

	// count remaining fusions
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.filter == FILTER_none)
			++remaining;
	return remaining;
}
