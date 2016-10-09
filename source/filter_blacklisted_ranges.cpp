#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "common.hpp"
#include "annotation.hpp"
#include "read_compressed_file.hpp"
#include "filter_blacklisted_ranges.hpp"

using namespace std;

bool parse_range(string range, const contigs_t& contigs, contig_t& contig, position_t& start, position_t& end) {
	istringstream iss;

	// extract contig from range
	string contig_name;
	replace(range.begin(), range.end(), ':', ' ');
	iss.str(range);
	iss >> contig_name;
	if (contig_name.empty()) {
		cerr << "WARNING: unknown gene or malformed range: " << range << endl;
		return false;
	}
	if (contigs.find(contig_name) == contigs.end()) {
		cerr << "WARNING: unknown gene or malformed range: " << contig_name << endl;
		return false;
	} else {
		contig = contigs.at(contig_name);
	}

	// extract start (and end) of range
	if (range.find("-") != string::npos) { // range has start and end (chr:start-end)
		replace(range.begin(), range.end(), '-', ' ');
		iss.str(range);
		iss >> contig_name; // discard contig
		if ((iss >> start).fail() || (iss >> end).fail()) {
			cerr << "WARNING: unknown gene or malformed range: " << range << endl;
			return false;
		}
		start--; // convert to zero-based coordinate
		end--;

	} else { // range is a single base (chr:position)
		if ((iss >> start).fail()) {
			cerr << "WARNING: unknown gene or malformed range: " << range << endl;
			return false;
		}
		start--; // convert to zero-based coordinate
		end = start;
	}

	return true;
}

bool blacklist_fusion(const contig_t contig1, const position_t start1, const position_t end1, const string& range1, const string& range2, const contigs_t& contigs, const unordered_map<string,gene_t>& genes, const contig_t fusion_contig1, const contig_t fusion_contig2, const position_t breakpoint1, const position_t breakpoint2, const direction_t direction1, const direction_t direction2, const gene_t gene1, const gene_t gene2, const unsigned int donor_split_reads, const unsigned int acceptor_split_reads, fusion_t* fusion) {

	if (genes.find(range1) != genes.end() && genes.at(range1) == gene1 || // match location by gene name
	    fusion_contig1 == contig1 && breakpoint1 >= start1 && breakpoint1 <= end1) { // match location by coordinate

		if (range2 == "any") { // remove the fusion if one breakpoint is within a region that is completely blacklisted
			fusion->filters.insert(FILTERS.at("blacklist"));
			return true;

		} else if (range2 == "split_read_donor") { // remove fusions which are only supported by blacklisted donor split reads
			if (fusion->discordant_mates + acceptor_split_reads == 0) {
				fusion->filters.insert(FILTERS.at("blacklist"));
				return true;
			}

		} else if (range2 == "split_read_acceptor") { // remove fusions which are only supported by blacklisted acceptor split reads
			if (fusion->discordant_mates + donor_split_reads == 0) {
				fusion->filters.insert(FILTERS.at("blacklist"));
				return true;
			}

		} else if (range2 == "split_read_any") { // remove fusions which are only supported by blacklisted split reads
			if (fusion->discordant_mates == 0) {
				fusion->filters.insert(FILTERS.at("blacklist"));
				return true;
			}

		} else if (range2 == "discordant_mates") { // remove fusions which are only supported by blacklisted discordant mates
			if (donor_split_reads + acceptor_split_reads == 0) {
				fusion->filters.insert(FILTERS.at("blacklist"));
				return true;
			}

		} else if (genes.find(range2) != genes.end()) { // remove fusions by name of gene
			if (genes.at(range2) == gene2) {
				fusion->filters.insert(FILTERS.at("blacklist"));
				return true;
			}

		} else { // range2 really contains a range (not a keyword)
			contig_t contig2;
			position_t start2, end2;
			if (parse_range(range2, contigs, contig2, start2, end2)) {
				if (fusion_contig2 == contig2 && breakpoint2 >= start2 && breakpoint2 <= end2) {
					fusion->filters.insert(FILTERS.at("blacklist"));
					return true;
				}
			}
		}
	}

	// if the fusion is supported only by discordant mates, then we discard
	// it, provided that the discordant mates point towards the blacklisted breakpoints
	if (donor_split_reads + acceptor_split_reads == 0 && // fusion is supported by discordant mates only
	    start1 == end1 && // blacklisted breakpoint is a single base position
	    fusion_contig1 == contig1 && (breakpoint1 <= start1 && breakpoint1 + 200 >= start1 && direction1 == DOWNSTREAM || breakpoint1 >= start1 && breakpoint1 <= start1 + 200 && direction1 == UPSTREAM) && // read1 points towards blacklisted breakpoint
	    genes.find(range1) == genes.end() && genes.find(range2) == genes.end() && // range1 and range2 do not contain a gene name
	    range2.find(":") != string::npos) { // range2 really contains a range and not a keyword
		contig_t contig2;
		position_t start2, end2;
		if (parse_range(range2, contigs, contig2, start2, end2)) {
			if (start2 == end2 && // blacklisted breakpont is a single base position
			    fusion_contig2 == contig2 && (breakpoint2 <= start2 && breakpoint2 + 200 >= start2 && direction2 == DOWNSTREAM || breakpoint2 >= start2 && breakpoint2 <= start2 + 200 && direction2 == UPSTREAM)) { // read2 points towards blacklisted breakpoint
				fusion->filters.insert(FILTERS.at("blacklist"));
				return true;
			}
		}
	}

	return false; // fusion was not filtered
}

unsigned int filter_blacklisted_ranges(fusions_t& fusions, const string blacklist_file_path, const contigs_t& contigs, const unordered_map<string,gene_t>& genes) {

	// sort fusions by coordinate of gene1
	map< position_t, vector<fusion_t*> > fusions_by_position;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
		if (!i->second.filters.empty())
			continue; // fusion has already been filtered
		fusions_by_position[i->second.breakpoint1].push_back(&i->second);
		fusions_by_position[i->second.breakpoint2].push_back(&i->second);
	}

	// hash fusions by gene
	unordered_map< gene_t, vector<fusion_t*> > fusions_by_gene;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
		if (!i->second.filters.empty())
			continue; // fusion has already been filtered
		fusions_by_gene[i->second.gene1].push_back(&i->second);
		fusions_by_gene[i->second.gene2].push_back(&i->second);
	}

	// load blacklist from file
	stringstream blacklist_file;
	autodecompress_file(blacklist_file_path, blacklist_file);
	string line;
	while (getline(blacklist_file, line)) {
		if (!line.empty() && line[0] != '#') {

			// parse line
			istringstream iss(line);
			string range1, range2;
			iss >> range1 >> range2;
			contig_t contig1;
			position_t start1, end1;

			// find genes falling into range1
			vector<fusion_t*> fusions_to_check;
			if (genes.find(range1) != genes.end()) { // range1 is a gene
				fusions_to_check = fusions_by_gene[genes.at(range1)];
			} else if (parse_range(range1, contigs, contig1, start1, end1)) { // range1 is a range
				// if the blacklisted range is a single base position, we loosen up the boundaries a little to
				// also catch discordant mates
				position_t loose_start1 = (start1 == end1) ? start1 - 200 : start1;
				position_t loose_end1 = (start1 == end1) ? end1 + 200 : end1;
				for (auto i = fusions_by_position.lower_bound(loose_start1); i != fusions_by_position.end() && i->first <= loose_end1; ++i)
					fusions_to_check.insert(fusions_to_check.end(), i->second.begin(), i->second.end());
			}

			// check for all fusions falling into range1 whether they also fall into range2
			for (auto fusion = fusions_to_check.begin(); fusion != fusions_to_check.end(); ++fusion) {

				if (!(**fusion).filters.empty())
					continue; // fusion has already been filtered

				// check if breakpoint1 is in range1 and breakpoint2 is in range2
				if (blacklist_fusion(contig1, start1, end1, range1, range2, contigs, genes, (**fusion).contig1, (**fusion).contig2, (**fusion).breakpoint1, (**fusion).breakpoint2, (**fusion).direction1, (**fusion).direction2, (**fusion).gene1, (**fusion).gene2, (**fusion).split_reads1, (**fusion).split_reads2, *fusion))
					continue;

				// check if breakpoint2 is in range1 and breakpoint1 is in range2
				if (blacklist_fusion(contig1, start1, end1, range1, range2, contigs, genes, (**fusion).contig2, (**fusion).contig1, (**fusion).breakpoint2, (**fusion).breakpoint1, (**fusion).direction2, (**fusion).direction1, (**fusion).gene2, (**fusion).gene1, (**fusion).split_reads2, (**fusion).split_reads1, *fusion))
					continue;

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

