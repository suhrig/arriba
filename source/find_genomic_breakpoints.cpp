#include <algorithm>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include "common.hpp"
#include "annotation.hpp"
#include "fusions.hpp"
#include "find_genomic_breakpoints.hpp"

using namespace std;

unsigned int find_genomic_breakpoints(fusions_t& fusions, annotation_t& gene_annotation) {

	typedef tuple<gene_t /*gene1*/, gene_t /*gene2*/, direction_t /*direction1*/, direction_t /*direction2*/> gene_pair_t;

	// for each pair of genes, find the breakpoints that have passed all filters
	// for these genes, we will look for additional breakpoints which were removed due to filtering
	unordered_map<gene_pair_t,fusions_t::iterator> best_breakpoints;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i)
		if (i->second.filters.empty())
			best_breakpoints[make_tuple(i->second.gene1, i->second.gene2, i->second.direction1, i->second.direction2)] = i;

	// for each pair of genes, look for breakpoints away from exon boundaries (non-spliced breakpoints)
	unordered_map< gene_pair_t, vector< fusions_t::iterator > > non_spliced_breakpoints;

	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		gene_pair_t gene_pair = make_tuple(i->second.gene1, i->second.gene2, i->second.direction1, i->second.direction2);

		if (i->second.filters.empty() || // the fusion has not been filtered and does not need to be recovered
		    best_breakpoints.find(gene_pair) == best_breakpoints.end()) // there are no trustworthy fusions between the given pair of genes => don't look for genomic breakpoints
			continue;

		// look for split-reads which are not spliced => they could point to the genomic breakpoint
		// non-spliced breakpoints must fulfill these conditions:
		// - it must be a split-read (discordant mates are not spliced by definition)
		// - the gene must be annotated (otherwise we have no information about splice sites)
		// - it must not be within 3bp of an exon boundary
		unsigned int split_reads = i->second.split_reads1 + i->second.split_reads2;
		if (split_reads == 0) { // look for discarded split reads
			for (auto chimeric_alignment = i->second.chimeric_alignments.begin(); chimeric_alignment != i->second.chimeric_alignments.end(); ++chimeric_alignment) {
				if ((**chimeric_alignment).size() == 3) { // is a split read
					if ((**chimeric_alignment).filters.find(FILTERS.at("same_gene")) == (**chimeric_alignment).filters.end() && // ignore intragenic split-reads (common in overlapping genes)
					    (**chimeric_alignment).filters.find(FILTERS.at("pcr_fusions")) == (**chimeric_alignment).filters.end() && // don't dig out potential PCR fusions (they would clutter the results list)
					    (**chimeric_alignment).filters.find(FILTERS.at("blacklist")) == (**chimeric_alignment).filters.end() && // don't dig out blacklisted genes (they would clutter the results list)
					    (**chimeric_alignment).filters.find(FILTERS.at("merge_adjacent")) == (**chimeric_alignment).filters.end()) { // don't dig out merged fusions (the one merged fusion is enough)
						split_reads++;
						break; // one split read is enough to call it a genomic breakpoint
					}
				}
			}
		}
		if (split_reads > 0 && // is a split read
		    !gene_annotation[i->second.gene1].is_dummy && // gene1 is annotated
		    !gene_annotation[i->second.gene2].is_dummy && // gene2 is annotated
		    !i->second.spliced1 && !i->second.spliced2) { // breakpoints are not near an exon boundary
			non_spliced_breakpoints[gene_pair].push_back(i);
		}

	}

	// recover non-spliced breakpoints
	for (auto i = non_spliced_breakpoints.begin(); i != non_spliced_breakpoints.end(); ++i) {
		for (auto j = i->second.begin(); j != i->second.end(); ++j) {
			if (abs((**j).second.breakpoint1 - best_breakpoints[i->first]->second.breakpoint1) > 10 &&
			    abs((**j).second.breakpoint2 - best_breakpoints[i->first]->second.breakpoint2) > 10) {
				(**j).second.filters.clear(); // recover discarded breakpoint
			}
		}
	}

	// count remaining fusions
	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i)
		if (i->second.filters.empty())
			remaining++;
	return remaining;
}

