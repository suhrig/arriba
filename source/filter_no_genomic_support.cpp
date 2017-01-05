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
#include "filter_no_genomic_support.hpp"

using namespace std;

void parse_breakpoint(string breakpoint, const contigs_t& contigs, contig_t& contig, position_t& position) {
	istringstream iss;

	// extract contig from breakpoint
	string contig_name;
	replace(breakpoint.begin(), breakpoint.end(), ':', ' ');
	iss.str(breakpoint);
	iss >> contig_name;
	if (contigs.find(contig_name) == contigs.end()) {
		cerr << "ERROR: unknown contig: " << contig_name << endl;
		exit(1);
	} else {
		contig = contigs.at(contig_name);
	}

	// extract position from breakpoint
	if ((iss >> position).fail()) {
		cerr << "ERROR: malformed breakpoint: " << breakpoint << endl;
		exit(1);
	}
	position--; // convert to zero-based coordinate
}

void parse_direction(const string& direction_string, direction_t& direction) {
	if (direction_string == "upstream") {
		direction = UPSTREAM;
	} else if (direction_string == "downstream") {
		direction = DOWNSTREAM;
	} else {
		cerr << "ERROR: invalid value for direction: " << direction_string << endl;
		exit(1);
	}
}

bool is_genomic_breakpoint_close_enough(const direction_t direction, const position_t genomic_breakpoint, const position_t fusion_breakpoint, const gene_t gene, const int max_distance) {
	// calculate most distal genomic position to still consider it as supporting
	position_t most_distal_genomic_position;
	if (direction == UPSTREAM) {
		if (gene->is_dummy)
			most_distal_genomic_position = fusion_breakpoint - max_distance;
		else
			most_distal_genomic_position = gene->start - max_distance;
		return genomic_breakpoint >= most_distal_genomic_position && genomic_breakpoint <= fusion_breakpoint + 5;
	} else {
		if (gene->is_dummy)
			most_distal_genomic_position = fusion_breakpoint + max_distance;
		else
			most_distal_genomic_position = gene->end + max_distance;
		return genomic_breakpoint <= most_distal_genomic_position && genomic_breakpoint >= fusion_breakpoint - 5;
	}
}

unsigned int mark_genomic_support(fusions_t& fusions, const string& genomic_breakpoints_file_path, const contigs_t& contigs, const int max_distance) {

	// make index structure for genomic breakpoints
	unordered_map< tuple<contig_t, contig_t, direction_t, direction_t>, map< position_t/*breakpoint1*/, vector<position_t/*breakpoint2*/> > > genomic_breakpoints;

	// load genomic breakpoints from file into index
	stringstream genomic_breakpoints_file;
	autodecompress_file(genomic_breakpoints_file_path, genomic_breakpoints_file);
	string line;
	while (getline(genomic_breakpoints_file, line)) {
		if (!line.empty() && line[0] != '#') {

			// parse line
			istringstream iss(line);
			string breakpoint1, breakpoint2;
			string string_direction1, string_direction2;
			iss >> breakpoint1 >> breakpoint2 >> string_direction1 >> string_direction2;
			contig_t contig1, contig2;
			position_t position1, position2;
			direction_t direction1, direction2;
			parse_breakpoint(breakpoint1, contigs, contig1, position1);
			parse_breakpoint(breakpoint2, contigs, contig2, position2);
			parse_direction(string_direction1, direction1);
			parse_direction(string_direction2, direction2);

			// make sure we index by the smaller coordinate
			if (contig2 < contig1 || contig2 == contig1 && position2 < position1) {
				swap(contig1, contig2);
				swap(position1, position2);
				swap(direction1, direction2);
			}

			// add genomic breakpoint to index
			genomic_breakpoints[make_tuple(contig1, contig2, direction1, direction2)][position1].push_back(position2);
		}
	}

	// for each fusion, check if it is supported by a genomic breakpoint
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		auto genomic_breakpoints_on_same_contigs = genomic_breakpoints.find(make_tuple(fusion->second.contig1, fusion->second.contig2, fusion->second.direction1, fusion->second.direction2));
		if (genomic_breakpoints_on_same_contigs != genomic_breakpoints.end()) {

			auto closeby_genomic_breakpoints = genomic_breakpoints_on_same_contigs->second.lower_bound(fusion->second.breakpoint1 + ((fusion->second.direction1 == UPSTREAM) ? +5 : -5)); // +/-5 allows for some alignment flexibility

			if (fusion->second.direction1 == UPSTREAM) {
				if (closeby_genomic_breakpoints == genomic_breakpoints_on_same_contigs->second.begin())
					continue; // there is no genomic breakpoint upstream of the fusion breakpoint
				--closeby_genomic_breakpoints; // move one genomic breakpoint upstream, since lower_bound() always returns the next downstream genomic breakpoint
			} else { // fusion->second.direction1 == DOWNSTREAM
				if (closeby_genomic_breakpoints == genomic_breakpoints_on_same_contigs->second.end())
					continue; // there is no genomic breakpoint downstream of the fusion breakpoint
			}

			// find the closest genomic breakpoints
			while (is_genomic_breakpoint_close_enough(fusion->second.direction1, closeby_genomic_breakpoints->first, fusion->second.breakpoint1, fusion->second.gene1, max_distance)) {

				for (auto closeby_genomic_breakpoint2 = closeby_genomic_breakpoints->second.begin(); closeby_genomic_breakpoint2 != closeby_genomic_breakpoints->second.end(); ++closeby_genomic_breakpoint2) {
					if (is_genomic_breakpoint_close_enough(fusion->second.direction2, *closeby_genomic_breakpoint2, fusion->second.breakpoint2, fusion->second.gene2, max_distance) &&
					    (fusion->second.contig1 != fusion->second.contig2 || // we need to make extra checks for deletions and inversions:
					     fusion->second.direction1 == UPSTREAM && fusion->second.direction2 == DOWNSTREAM || // (but not duplications)
					     fusion->second.direction1 == DOWNSTREAM && fusion->second.direction2 == UPSTREAM && closeby_genomic_breakpoints->first < fusion->second.breakpoint2 && *closeby_genomic_breakpoint2 > fusion->second.breakpoint1 || // for deletions, both genomic breakpoints must be between the transcriptomic breakpoints
					     fusion->second.direction1 == UPSTREAM && fusion->second.direction2 == UPSTREAM && *closeby_genomic_breakpoint2 > fusion->second.breakpoint1 || // for inversions, one genomic breakpoint must be between the transcriptomic breakpoints
					     fusion->second.direction1 == DOWNSTREAM && fusion->second.direction2 == DOWNSTREAM && closeby_genomic_breakpoints->first < fusion->second.breakpoint2)) { // for inversions, one genomic breakpoint must be between the transcriptomic breakpoints
					                                                                                                                                      // (this avoids false associations in the case of small deletions/inversions)
						// we consider a pair of genomic breakpoints to be closer than a given one,
						// if the sum of the distances between genomic and transcriptomic breakpoints is lower
						if (fusion->second.closest_genomic_breakpoint1 < 0 || fusion->second.closest_genomic_breakpoint2 < 0 ||
						    abs(fusion->second.breakpoint1 - fusion->second.closest_genomic_breakpoint1) + abs(fusion->second.breakpoint2 - fusion->second.closest_genomic_breakpoint2) > abs(closeby_genomic_breakpoints->first - fusion->second.breakpoint1) + abs(fusion->second.breakpoint2 - *closeby_genomic_breakpoint2)) {
							fusion->second.closest_genomic_breakpoint1 = closeby_genomic_breakpoints->first;
							fusion->second.closest_genomic_breakpoint2 = *closeby_genomic_breakpoint2;
						}
					}
				}

				if (closeby_genomic_breakpoints != genomic_breakpoints_on_same_contigs->second.begin())
					--closeby_genomic_breakpoints;
				else
					break;
			}

		}
	}

	// count number of fusions with genomic support
	unsigned int marked = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (fusion->second.closest_genomic_breakpoint1 >= 0)
			marked++;
	return marked;
}

void assign_confidence(fusions_t& fusions) {

	// order fusions by gene pair
	// => the confidence in an event is increased, when there are other events between the same pair of genes
	map< gene_t, vector<fusion_t*> > fusions_by_gene;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		fusions_by_gene[fusion->second.gene1].push_back(&fusion->second);
		fusions_by_gene[fusion->second.gene2].push_back(&fusion->second);
	}

	// assign a confidence to each fusion
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != NULL) {
			fusion->second.confidence = CONFIDENCE_LOW;
		} else {
			if (fusion->second.evalue > 0.3) {
				fusion->second.confidence = CONFIDENCE_LOW;
				if (fusion->second.spliced1 && fusion->second.spliced2 && !fusion->second.is_read_through()) {
					// look for multiple spliced translocations
					unsigned int number_of_spliced_breakpoints = 0;
					auto fusions_of_gene = fusions_by_gene.find(fusion->second.gene1);
					for (auto fusion_of_gene = fusions_of_gene->second.begin(); fusion_of_gene != fusions_of_gene->second.end(); ++fusion_of_gene) {
						if ((**fusion_of_gene).gene1 == fusion->second.gene1 && (**fusion_of_gene).gene2 == fusion->second.gene2 && (**fusion_of_gene).spliced1 && (**fusion_of_gene).spliced2)
							++number_of_spliced_breakpoints;
					}
					fusions_of_gene = fusions_by_gene.find(fusion->second.gene2);
					for (auto fusion_of_gene = fusions_of_gene->second.begin(); fusion_of_gene != fusions_of_gene->second.end(); ++fusion_of_gene) {
						if ((**fusion_of_gene).gene1 == fusion->second.gene1 && (**fusion_of_gene).gene2 == fusion->second.gene2 && (**fusion_of_gene).spliced1 && (**fusion_of_gene).spliced2)
							++number_of_spliced_breakpoints;
					}
					if (number_of_spliced_breakpoints >= 2)
						fusion->second.confidence = CONFIDENCE_MEDIUM;
				}
			} else if (fusion->second.is_read_through()) {
				fusion->second.confidence = CONFIDENCE_LOW;
				if ((fusion->second.split_reads1 > 0 && fusion->second.split_reads2 > 0 || fusion->second.split_reads1 > 0 && fusion->second.discordant_mates > 0 || fusion->second.split_reads2 > 0 && fusion->second.discordant_mates > 0) && fusion->second.supporting_reads() > 9) {
					fusion->second.confidence = CONFIDENCE_MEDIUM;
				} else {
					// look for multiple deletions involving the same gene
					unsigned int number_of_deletions = 0;
					auto fusions_of_gene = fusions_by_gene.find(fusion->second.gene1);
					for (auto fusion_of_gene = fusions_of_gene->second.begin(); fusion_of_gene != fusions_of_gene->second.end(); ++fusion_of_gene) {
						if ((**fusion_of_gene).filter == NULL &&
						    (**fusion_of_gene).split_reads1 + (**fusion_of_gene).split_reads2 > 0 &&
						    (**fusion_of_gene).direction1 == DOWNSTREAM && (**fusion_of_gene).direction2 == UPSTREAM &&
						    ((**fusion_of_gene).gene1 == fusion->second.gene1 && (**fusion_of_gene).gene2 != fusion->second.gene2 || // don't count different isoforms
						     (**fusion_of_gene).gene1 != fusion->second.gene1 && (**fusion_of_gene).gene2 == fusion->second.gene2) &&
						    ((**fusion_of_gene).breakpoint1 != fusion->second.breakpoint1 || (**fusion_of_gene).breakpoint2 != fusion->second.breakpoint2) &&
						    (**fusion_of_gene).breakpoint2 > fusion->second.breakpoint1 && (**fusion_of_gene).breakpoint1 < fusion->second.breakpoint2) {
							++number_of_deletions;
						}
					}
					fusions_of_gene = fusions_by_gene.find(fusion->second.gene2);
					for (auto fusion_of_gene = fusions_of_gene->second.begin(); fusion_of_gene != fusions_of_gene->second.end(); ++fusion_of_gene) {
						if ((**fusion_of_gene).filter == NULL &&
						    (**fusion_of_gene).split_reads1 + (**fusion_of_gene).split_reads2 > 0 &&
						    (**fusion_of_gene).direction1 == DOWNSTREAM && (**fusion_of_gene).direction2 == UPSTREAM &&
						    ((**fusion_of_gene).gene1 == fusion->second.gene1 && (**fusion_of_gene).gene2 != fusion->second.gene2 || // don't count different isoforms
						     (**fusion_of_gene).gene1 != fusion->second.gene1 && (**fusion_of_gene).gene2 == fusion->second.gene2) &&
						    ((**fusion_of_gene).breakpoint1 != fusion->second.breakpoint1 || (**fusion_of_gene).breakpoint2 != fusion->second.breakpoint2) &&
						    (**fusion_of_gene).breakpoint2 > fusion->second.breakpoint1 && (**fusion_of_gene).breakpoint1 < fusion->second.breakpoint2) {
							++number_of_deletions;
							}
					}
					if (number_of_deletions >= 1)
						fusion->second.confidence = CONFIDENCE_MEDIUM;
				}

			} else if (fusion->second.breakpoint_overlaps_both_genes() || fusion->second.gene1 == fusion->second.gene2) { // intragenic event
				if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0) {
					fusion->second.confidence = CONFIDENCE_LOW;
				} else if (!fusion->second.exonic1 && !fusion->second.exonic2) {
					if (fusion->second.split_reads1 > 0 && fusion->second.split_reads2 > 0) {
						fusion->second.confidence = CONFIDENCE_HIGH;
					} else {
						fusion->second.confidence = CONFIDENCE_MEDIUM;
					}
				} else if (!fusion->second.exonic1 || !fusion->second.exonic2) {
					if (fusion->second.split_reads1 > 3 && fusion->second.split_reads2 > 3) {
						fusion->second.confidence = CONFIDENCE_HIGH;
					} else {
						fusion->second.confidence = CONFIDENCE_MEDIUM;
					}
				} else {
					fusion->second.confidence = CONFIDENCE_LOW;
				}
			} else if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0 || fusion->second.split_reads1 + fusion->second.discordant_mates == 0 || fusion->second.split_reads2 + fusion->second.discordant_mates == 0) {
				fusion->second.confidence = CONFIDENCE_MEDIUM;
			} else {
				fusion->second.confidence = CONFIDENCE_HIGH;
			}

			if (fusion->second.evalue > 0.2) { // decrease the confidence, when the e-value is not overwhelming
				if (fusion->second.confidence == CONFIDENCE_HIGH)
					fusion->second.confidence = CONFIDENCE_MEDIUM;
				else if (fusion->second.confidence == CONFIDENCE_MEDIUM)
					fusion->second.confidence = CONFIDENCE_LOW;
			}

			if (fusion->second.closest_genomic_breakpoint1 >= 0) { // has genomic support
				if (fusion->second.confidence == CONFIDENCE_LOW)
					fusion->second.confidence = CONFIDENCE_MEDIUM;
				else if (fusion->second.confidence == CONFIDENCE_MEDIUM)
					fusion->second.confidence = CONFIDENCE_HIGH;
			}
		}
	}
}

// filter speculative fusions without support from WGS
unsigned int filter_no_genomic_support(fusions_t& fusions) {

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.closest_genomic_breakpoint1 < 0 && // no genomic support
		     fusion->second.confidence == CONFIDENCE_LOW)
			fusion->second.filter = FILTERS.at("no_genomic_support");
		else
			remaining++;
	}

	return remaining;
}

