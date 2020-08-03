#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include "common.hpp"
#include "annotation.hpp"
#include "read_compressed_file.hpp"
#include "filter_blacklisted_ranges.hpp"

using namespace std;

// convert string representation of a range into coordinates
bool parse_range(string range, const contigs_t& contigs, blacklist_item_t& blacklist_item) {

	// extract contig from range
	tsv_stream_t tsv(range, ':');
	string contig_name, start_and_end_position;
	tsv >> contig_name >> start_and_end_position;
	if (tsv.fail() || contig_name.empty() || start_and_end_position.empty()) {
		cerr << "WARNING: unknown gene or malformed range: " << range << endl;
		return false;
	}
	if (contig_name[0] == '+') {
		blacklist_item.strand_defined = true;
		blacklist_item.strand = FORWARD;
		contig_name = contig_name.substr(1);
	} else if (contig_name[0] == '-') {
		blacklist_item.strand_defined = true;
		blacklist_item.strand = REVERSE;
		contig_name = contig_name.substr(1);
	} else {
		blacklist_item.strand_defined = false;
	}
	contig_name = removeChr(contig_name);
	if (contigs.find(contig_name) == contigs.end()) {
		cerr << "WARNING: unknown gene or malformed range: " << range << endl;
		return false;
	} else {
		blacklist_item.contig = contigs.at(contig_name);
	}

	// extract start (and end) of range
	tsv_stream_t tsv2(start_and_end_position, '-');
	if (start_and_end_position.find('-') != string::npos) { // range has start and end (chr:start-end)
		if ((tsv2 >> blacklist_item.start >> blacklist_item.end).fail()) {
			cerr << "WARNING: unknown gene or malformed range: " << range << endl;
			return false;
		}
		blacklist_item.start--; // convert to zero-based coordinate
		blacklist_item.end--;

	} else { // range is a single base (chr:position)
		if ((tsv2 >> blacklist_item.start).fail()) {
			cerr << "WARNING: unknown gene or malformed range: " << range << endl;
			return false;
		}
		blacklist_item.start--; // convert to zero-based coordinate
		blacklist_item.end = blacklist_item.start;
	}

	return true;
}

// parse string representation of a blacklist rule
bool parse_blacklist_item(const string& text, blacklist_item_t& blacklist_item, const contigs_t& contigs, const unordered_map<string,gene_t>& genes, const bool allow_keyword) {

	if (allow_keyword) {
		if (text == "any") { blacklist_item.type = BLACKLIST_ANY; return true; }
		else if (text == "split_read_donor") { blacklist_item.type = BLACKLIST_SPLIT_READ_DONOR; return true; }
		else if (text == "split_read_acceptor") { blacklist_item.type = BLACKLIST_SPLIT_READ_ACCEPTOR; return true; }
		else if (text == "split_read_any") { blacklist_item.type = BLACKLIST_SPLIT_READ_ANY; return true; }
		else if (text == "discordant_mates") { blacklist_item.type = BLACKLIST_DISCORDANT_MATES; return true; }
		else if (text == "read_through") { blacklist_item.type = BLACKLIST_READ_THROUGH; return true; }
		else if (text == "low_support") { blacklist_item.type = BLACKLIST_LOW_SUPPORT; return true; }
		else if (text == "filter_spliced") { blacklist_item.type = BLACKLIST_FILTER_SPLICED; return true; }
		else if (text == "not_both_spliced") { blacklist_item.type = BLACKLIST_NOT_BOTH_SPLICED; return true; }
	}

	auto gene = genes.find(text); // check if text is a gene name
	if (gene != genes.end()) {
		blacklist_item.type = BLACKLIST_GENE;
		blacklist_item.gene = gene->second;
		blacklist_item.contig = gene->second->contig;
		blacklist_item.start = gene->second->start;
		blacklist_item.end = gene->second->end;
	} else if (parse_range(text, contigs, blacklist_item)) { // text is a range
		if (blacklist_item.start == blacklist_item.end) {
			blacklist_item.type = BLACKLIST_POSITION;
		} else {
			blacklist_item.type = BLACKLIST_RANGE;
		}
	} else {
		return false; // failed to parse text
	}
	return true; // parsed text successfully

}

// returns the fraction of range1 that overlaps range2
float overlapping_fraction(const position_t start1, const position_t end1, const position_t start2, const position_t end2) {
	if (start1 >= start2 && end1 <= end2) {
		return 1;
	} else if (start1 < start2 && end1 > end2) {
		return 1.0 * (end2 - start2) / (end1 - start1 + 1);
	} else if (start1 >= start2 && start1 <= end2) {
		return 1.0 * (end2 - start1) / (end1 - start1 + 1);
	} else if (end1 >= start2 && end1 <= end2) {
		return 1.0 * (end1 - start2) / (end1 - start1 + 1);
	} else {
		return 0;
	}
}

// check if the breakpoint of a fusion match an entry in the blacklist
bool matches_blacklist_item(const blacklist_item_t& blacklist_item, const fusion_t& fusion, const unsigned char which_breakpoint, const int max_mate_gap, const float evalue_cutoff) {

	switch (blacklist_item.type) {
		case BLACKLIST_ANY: // remove the fusion if one breakpoint is within a region that is completely blacklisted
			return true;
		case BLACKLIST_SPLIT_READ_DONOR: // remove fusions which are only supported by donor split reads
			return (which_breakpoint == 1 && fusion.discordant_mates + fusion.split_reads1 == 0 ||
			        which_breakpoint == 2 && fusion.discordant_mates + fusion.split_reads2 == 0);
		case BLACKLIST_SPLIT_READ_ACCEPTOR: // remove fusions which are only supported by acceptor split reads
			return (which_breakpoint == 1 && fusion.discordant_mates + fusion.split_reads2 == 0 ||
			        which_breakpoint == 2 && fusion.discordant_mates + fusion.split_reads1 == 0);
		case BLACKLIST_SPLIT_READ_ANY: // remove fusions which are only supported by split reads
			return (fusion.discordant_mates == 0);
		case BLACKLIST_DISCORDANT_MATES: // remove fusions which are only supported by discordant mates
			return (fusion.split_reads1 + fusion.split_reads2 == 0);
		case BLACKLIST_READ_THROUGH: // remove read-through fusions
			return fusion.is_read_through();
		case BLACKLIST_LOW_SUPPORT: // remove recurrent speculative fusions that were recovered for one or the other reason
			return (fusion.evalue > evalue_cutoff);
		case BLACKLIST_FILTER_SPLICED: // remove recurrent speculative fusions that were recovered by the 'spliced' filter
			return (fusion.evalue > evalue_cutoff && fusion.spliced1 && fusion.spliced2);
		case BLACKLIST_NOT_BOTH_SPLICED: // remove fusions which do not have both breakpoints at splice-sites
			return (!fusion.spliced1 || !fusion.spliced2);
		case BLACKLIST_GENE: // remove blacklisted gene
			return (which_breakpoint == 1 && fusion.gene1 == blacklist_item.gene ||
			        which_breakpoint == 2 && fusion.gene2 == blacklist_item.gene);
		case BLACKLIST_POSITION: { // remove blacklisted breakpoint

			// contig must match
			contig_t contig = (which_breakpoint == 1) ? fusion.contig1 : fusion.contig2;
			if (contig != blacklist_item.contig)
				return false;

			// strand must match, if defined
			if (blacklist_item.strand_defined) {
				if (!fusion.predicted_strands_ambiguous) { // assume match, if strands could not be predicted
					strand_t strand = (which_breakpoint == 1) ? fusion.predicted_strand1 : fusion.predicted_strand2;
					if (strand != blacklist_item.strand)
						return false;
				}
			}

			// exact breakpoint must match
			position_t breakpoint = (which_breakpoint == 1) ? fusion.breakpoint1 : fusion.breakpoint2;
			if (breakpoint == blacklist_item.start)
				return true;

			// if the fusion has no split reads, then we discard it, if the discordant mates are near a blacklisted breakpoint
			// and point towards it
			if (fusion.split_reads1 + fusion.split_reads2 == 0) {
				direction_t direction = (which_breakpoint == 1) ? fusion.direction1 : fusion.direction2;
				if (direction == DOWNSTREAM && breakpoint <= blacklist_item.start && breakpoint >= blacklist_item.start - max_mate_gap ||
				    direction == UPSTREAM   && breakpoint >= blacklist_item.start && breakpoint <= blacklist_item.start + max_mate_gap)
					return true;
			}

			return false; // blacklist item does not match
		}

		case BLACKLIST_RANGE: { // remove blacklisted range

			// contig must match
			contig_t contig = (which_breakpoint == 1) ? fusion.contig1 : fusion.contig2;
			if (contig != blacklist_item.contig)
				return false;

			// strand must match, if defined
			if (blacklist_item.strand_defined) {
				if (!fusion.predicted_strands_ambiguous) { // assume match, if strands could not be predicted
					strand_t strand = (which_breakpoint == 1) ? fusion.predicted_strand1 : fusion.predicted_strand2;
					if (strand != blacklist_item.strand)
						return false;
				}
			}

			// check if the gene that the breakpoint is associated with overlaps the blacklisted range
			gene_t gene = (which_breakpoint == 1) ? fusion.gene1 : fusion.gene2;
			if (overlapping_fraction(gene->start, gene->end, blacklist_item.start, blacklist_item.end) > 0.5)
				return true;

			return false; // blacklist item does not match
		}
	}

	return false; // blacklist item does not match
}

// divide the genome into fixed size bins
void get_genome_bins_from_range(const contig_t contig, const position_t start, const position_t end, genome_bins_t& genome_bins) {
	const int bucket_size = 100000; // bp
	for (position_t position = start/bucket_size; position <= (end+bucket_size-1)/bucket_size /*integer ceil*/; ++position)
		genome_bins.push_back(make_tuple(contig, position*bucket_size));
}

unsigned int filter_blacklisted_ranges(fusions_t& fusions, const string& blacklist_file_path, const contigs_t& contigs, const unordered_map<string,gene_t>& genes, const float evalue_cutoff, const int max_mate_gap) {

	// index fusions by coordinate
	unordered_map< tuple<contig_t,position_t>, set<fusion_t*> > fusions_by_coordinate;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != FILTER_none && fusion->second.closest_genomic_breakpoint1 < 0)
			continue; // fusion has already been filtered and won't be recovered by the 'genomic_support' filter

		// assign fusions within a window of ...bps to the same bin
		// for fast lookup of fusions by coordinate
		genome_bins_t genome_bins;
		get_genome_bins_from_range(fusion->second.contig1, fusion->second.breakpoint1, fusion->second.breakpoint1, genome_bins);
		get_genome_bins_from_range(fusion->second.contig2, fusion->second.breakpoint2, fusion->second.breakpoint2, genome_bins);
		get_genome_bins_from_range(fusion->second.contig1, fusion->second.gene1->start, fusion->second.gene1->end, genome_bins);
		get_genome_bins_from_range(fusion->second.contig2, fusion->second.gene2->start, fusion->second.gene2->end, genome_bins);
		for (auto genome_bin = genome_bins.begin(); genome_bin != genome_bins.end(); ++genome_bin)
			fusions_by_coordinate[*genome_bin].insert(&(fusion->second));
	}

	// load blacklist from file
	autodecompress_file_t blacklist_file(blacklist_file_path);
	string line;
	while (blacklist_file.getline(line)) {

		// skip comment lines
		if (line.empty() || line[0] == '#')
			continue;

		// parse line
		tsv_stream_t tsv(line);
		string range1, range2;
		tsv >> range1 >> range2;
		blacklist_item_t item1, item2;
		if (!parse_blacklist_item(range1, item1, contigs, genes, false) ||
		    !parse_blacklist_item(range2, item2, contigs, genes, true))
			continue;

		// find all fusions with breakpoints in the vicinity of the blacklist items
		genome_bins_t genome_bins;
		if (item1.type == BLACKLIST_POSITION || item1.type == BLACKLIST_RANGE || item1.type == BLACKLIST_GENE)
			get_genome_bins_from_range(item1.contig, item1.start-max_mate_gap, item1.end+max_mate_gap, genome_bins);
		if (item2.type == BLACKLIST_POSITION || item2.type == BLACKLIST_RANGE || item2.type == BLACKLIST_GENE)
			get_genome_bins_from_range(item2.contig, item2.start-max_mate_gap, item2.end+max_mate_gap, genome_bins);
		for (auto genome_bin = genome_bins.begin(); genome_bin != genome_bins.end(); ++genome_bin) {
			auto fusions_near_coordinate = fusions_by_coordinate.find(*genome_bin);
			if (fusions_near_coordinate != fusions_by_coordinate.end()) {
				for (auto fusion = fusions_near_coordinate->second.begin(); fusion != fusions_near_coordinate->second.end();) {
					if (matches_blacklist_item(item1, **fusion, 1, max_mate_gap, evalue_cutoff) &&
					    matches_blacklist_item(item2, **fusion, 2, max_mate_gap, evalue_cutoff) ||
					    matches_blacklist_item(item1, **fusion, 2, max_mate_gap, evalue_cutoff) &&
					    matches_blacklist_item(item2, **fusion, 1, max_mate_gap, evalue_cutoff)) {
						(**fusion).filter = FILTER_blacklist;
						fusions_near_coordinate->second.erase(fusion++); // remove fusion from index, so we don't check it again
					} else {
						++fusion;
					}
				}
			}
		}
	}

	// count remaining fusions
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (!fusion->second.filter != FILTER_none)
			remaining++;
	return remaining;
}

