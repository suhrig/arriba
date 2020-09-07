#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
#include "common.hpp"
#include "annotation.hpp"
#include "read_compressed_file.hpp"
#include "read_stats.hpp"
#include "filter_genomic_support.hpp"

using namespace std;

bool parse_breakpoint(string breakpoint, const contigs_t& contigs, contig_t& contig, position_t& position) {
	unsigned int separator = breakpoint.find_last_of(':'); // split by last colon, because there could be colons in the contig name
	if (separator != string::npos)
		breakpoint[separator] = '\t';
	tsv_stream_t tsv(breakpoint);

	// extract contig from breakpoint
	string contig_name;
	tsv >> contig_name;
	contig_name = removeChr(contig_name);
	if (contigs.find(contig_name) == contigs.end())
		return false;
	contig = contigs.at(contig_name);

	// extract position from breakpoint
	if ((tsv >> position).fail())
		return false;
	position--; // convert to zero-based coordinate

	return true;
}

bool parse_direction(const string& direction_string, direction_t& direction) {
	if (direction_string == "upstream" || direction_string == "-") {
		direction = UPSTREAM;
	} else if (direction_string == "downstream" || direction_string == "+") {
		direction = DOWNSTREAM;
	} else {
		return false;
	}
	return true;
}

bool parse_vcf_info(const string& vcf_info, const string& field, string& value) {
	size_t start;
	if (vcf_info.substr(0, field.size() + 1) == field + "=") { // field is the first one in VCF INFO column
		start = field.size() + 1; // +1 for "="
	} else {
		start = vcf_info.find(";" + field + "="); // if it's not the first, there must be a semi-colon preceding it
		if (start == string::npos)
			return false;
		start += field.size() + 2; // +2 for ";" and "="
	}
	value = vcf_info.substr(start, vcf_info.find(';', start) - start);
	return true;
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
	autodecompress_file_t genomic_breakpoints_file(genomic_breakpoints_file_path);
	string line;
	while (genomic_breakpoints_file.getline(line)) {
		if (!line.empty() && line[0] != '#') {

			// try to parse line as Arriba's four-column format (contig1:position1\tcontig2:position2\tdirection1\tdirection2)
			tsv_stream_t tsv(line);
			string breakpoint1, breakpoint2;
			string string_direction1, string_direction2;
			tsv >> breakpoint1 >> breakpoint2 >> string_direction1 >> string_direction2;
			contig_t contig1, contig2;
			position_t position1, position2;
			direction_t direction1, direction2;
			string vcf_sv_type = "";
			if (!(parse_breakpoint(breakpoint1, contigs, contig1, position1) &&
			      parse_breakpoint(breakpoint2, contigs, contig2, position2) &&
			      parse_direction(string_direction1, direction1) &&
			      parse_direction(string_direction2, direction2))) {

				// parsing as Arriba's four-column format failed => try VCF
				tsv_stream_t tsv2(line);
				string vcf_chrom, vcf_pos, vcf_alt, vcf_info, vcf_filter, ignore;
				tsv2 >> vcf_chrom >> vcf_pos >> ignore >> ignore >> vcf_alt >> ignore >> vcf_filter >> vcf_info;
				if (!parse_vcf_info(vcf_info, "SVTYPE", vcf_sv_type))
					goto failed_to_parse_line;
				if (vcf_sv_type == "BND") {
					size_t opening_bracket = vcf_alt.find('[');
					size_t closing_bracket = vcf_alt.find(']');
					char bracket = (opening_bracket < closing_bracket) ? '[' : ']';
					size_t bracket_pos1 = min(opening_bracket, closing_bracket);
					size_t bracket_pos2 = vcf_alt.find(bracket, bracket_pos1 + 1);
					if (bracket_pos1 == string::npos || bracket_pos2 == string::npos)
						if (!vcf_alt.empty() && (vcf_alt[0] == '.' || vcf_alt[vcf_alt.size()-1] == '.')) // is it a single breakend?
							continue; // silently ignore single breakend
						else
							goto failed_to_parse_line;
					direction1 = (bracket_pos1 == 0) ? UPSTREAM : DOWNSTREAM;
					direction2 = (bracket == '[') ? UPSTREAM : DOWNSTREAM;
					breakpoint2 = vcf_alt.substr(bracket_pos1 + 1, bracket_pos2 - bracket_pos1 - 1);
				} else {
					string vcf_info_end;
					if (!parse_vcf_info(vcf_info, "END", vcf_info_end))
						goto failed_to_parse_line;
					breakpoint2 = vcf_chrom + ":" + vcf_info_end;
					if (vcf_sv_type == "INV") {
						direction1 = DOWNSTREAM;
						direction2 = DOWNSTREAM;
					} else if (vcf_sv_type == "DEL") {
						direction1 = DOWNSTREAM;
						direction2 = UPSTREAM;
					} else if (vcf_sv_type == "DUP") {
						direction1 = UPSTREAM;
						direction2 = DOWNSTREAM;
					} else
						goto failed_to_parse_line;
				}
				if (!parse_breakpoint(vcf_chrom + ":" + vcf_pos, contigs, contig1, position1) ||
				    !parse_breakpoint(breakpoint2, contigs, contig2, position2))
					goto failed_to_parse_line;

				if (vcf_filter != "PASS")
					continue;
			}

			// make sure we index by the smaller coordinate
			if (contig2 < contig1 || contig2 == contig1 && position2 < position1) {
				swap(contig1, contig2);
				swap(position1, position2);
				swap(direction1, direction2);
			}

			// add genomic breakpoint to index
			genomic_breakpoints[make_tuple(contig1, contig2, direction1, direction2)][position1].push_back(position2);
			// the VCF SVTYPE "INV" encodes two separate breakpoints
			if (vcf_sv_type == "INV")
				genomic_breakpoints[make_tuple(contig1, contig2, UPSTREAM, UPSTREAM)][position1].push_back(position2);
		}

		continue;
		failed_to_parse_line:
			cerr << "WARNING: failed to parse line: " << line << endl;
	}

	// for each fusion, check if it is supported by a genomic breakpoint
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		auto genomic_breakpoints_on_same_contigs = genomic_breakpoints.find(make_tuple(fusion->second.contig1, fusion->second.contig2, (direction_t) fusion->second.direction1, (direction_t) fusion->second.direction2));
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

void assign_confidence(fusions_t& fusions, const coverage_t& coverage) {

	// order fusions by gene pair
	// => the confidence in an event is increased, when there are other events between the same pair of genes
	map< gene_t, vector<fusion_t*> > fusions_by_gene;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		fusions_by_gene[fusion->second.gene1].push_back(&fusion->second);
		fusions_by_gene[fusion->second.gene2].push_back(&fusion->second);
	}

	// assign a confidence to each fusion
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		// calculate supporting reads as fraction of coverage
		// we use the sizes of the read lists rather than supporting_reads() because the latter ignores duplicates, but the coverage does not
		int coverage1 = coverage.get_coverage(fusion->second.contig1, fusion->second.breakpoint1, (fusion->second.direction1 == UPSTREAM) ? DOWNSTREAM : UPSTREAM);
		int coverage2 = coverage.get_coverage(fusion->second.contig2, fusion->second.breakpoint2, (fusion->second.direction2 == UPSTREAM) ? DOWNSTREAM : UPSTREAM);
		float coverage_fraction = ((float) (fusion->second.split_read1_list.size() + fusion->second.split_read2_list.size() + fusion->second.discordant_mate_list.size())) / max(1, max(coverage1, coverage2));

		if (fusion->second.filter != FILTER_none) {
			// discarded events get low confidence, no matter what
			fusion->second.confidence = CONFIDENCE_LOW;

		} else {

			// non-discarded events get high confidence by default, which my be reduced by penalties
			fusion->second.confidence = CONFIDENCE_HIGH;

			if (fusion->second.evalue > 0.3 || fusion->second.supporting_reads() < 2) {
				fusion->second.confidence = CONFIDENCE_LOW; // events with poor support get low confidence, regardless of the location

			} else if (fusion->second.is_read_through()) {
				// read-through events get low confidence by default,
				// but this may improve, if they have good support or
				// if there are many read-through fusions all hinting at the same genomic event
				fusion->second.confidence = CONFIDENCE_LOW;

				// increase confidence of read-through fusion, if is has many supporting reads
				if ((fusion->second.split_reads1 > 0 && fusion->second.split_reads2 > 0 ||
				     fusion->second.split_reads1 > 0 && fusion->second.discordant_mates > 0 ||
				     fusion->second.split_reads2 > 0 && fusion->second.discordant_mates > 0) &&
				     fusion->second.supporting_reads() >= 10) {
					if (fusion->second.split_reads1 + fusion->second.split_reads2 >= 10 && coverage_fraction > 0.07)
						fusion->second.confidence = CONFIDENCE_HIGH;
					else
						fusion->second.confidence = CONFIDENCE_MEDIUM;

				} else {
					// look for multiple deletions involving the same gene
					unsigned int number_of_deletions = 0;
					auto fusions_of_gene = fusions_by_gene.find(fusion->second.gene1);
					for (auto fusion_of_gene = fusions_of_gene->second.begin(); fusion_of_gene != fusions_of_gene->second.end(); ++fusion_of_gene) {
						if ((**fusion_of_gene).filter == FILTER_none &&
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
						if ((**fusion_of_gene).filter == FILTER_none &&
						    (**fusion_of_gene).split_reads1 + (**fusion_of_gene).split_reads2 > 0 &&
						    (**fusion_of_gene).direction1 == DOWNSTREAM && (**fusion_of_gene).direction2 == UPSTREAM &&
						    ((**fusion_of_gene).gene1 == fusion->second.gene1 && (**fusion_of_gene).gene2 != fusion->second.gene2 || // don't count different isoforms
						     (**fusion_of_gene).gene1 != fusion->second.gene1 && (**fusion_of_gene).gene2 == fusion->second.gene2) &&
						    ((**fusion_of_gene).breakpoint1 != fusion->second.breakpoint1 || (**fusion_of_gene).breakpoint2 != fusion->second.breakpoint2) &&
						    (**fusion_of_gene).breakpoint2 > fusion->second.breakpoint1 && (**fusion_of_gene).breakpoint1 < fusion->second.breakpoint2) {
							++number_of_deletions;
							}
					}
					if (number_of_deletions >= 1) // there is at least one additional fusion which may arise from the same deletion
						fusion->second.confidence = CONFIDENCE_MEDIUM;
				}

			} else if (fusion->second.breakpoint_overlaps_both_genes() || fusion->second.gene1 == fusion->second.gene2) { // intragenic event

				// by default intragenic events get low confidence, because there are so many artifacts
				fusion->second.confidence = CONFIDENCE_LOW;

				if (fusion->second.split_reads1 + fusion->second.split_reads2 > 0) {
					if (!fusion->second.exonic1 && !fusion->second.exonic2) {
						// most intragenic artifacts have both breakpoints in exons
						// => increase confidence a bit if the breakpoints are in introns
						if (fusion->second.split_reads1 > 0 && fusion->second.split_reads2 > 0)
							fusion->second.confidence = CONFIDENCE_HIGH;
						else
							fusion->second.confidence = CONFIDENCE_MEDIUM;
					} else if (!fusion->second.exonic1 || !fusion->second.exonic2) {
						// one breakpoint is in an intron, the other isn't
						// => also increase confidence a bit, but require some more split reads
						if (fusion->second.split_reads1 > 3 && fusion->second.split_reads2 > 3)
							fusion->second.confidence = CONFIDENCE_HIGH;
						else
							fusion->second.confidence = CONFIDENCE_MEDIUM;
					}
				}

			}

			// lift confidence score of rescued internal tandem duplications
			if (fusion->second.confidence == CONFIDENCE_LOW &&
			    fusion->second.gene1 == fusion->second.gene2 &&
			    fusion->second.exonic1 && fusion->second.exonic2 &&
			    !fusion->second.spliced1 && !fusion->second.spliced2 &&
			    fusion->second.breakpoint2 - fusion->second.breakpoint1 < 100 &&
			    fusion->second.split_reads1 > 0 && fusion->second.split_reads2 > 0 &&
			    fusion->second.split_reads1 + fusion->second.split_reads2 >= 10 &&
			    coverage_fraction > 0.15 &&
			    fusion->second.direction1 == UPSTREAM && fusion->second.direction2 == DOWNSTREAM)
				fusion->second.confidence = CONFIDENCE_MEDIUM;

			// increase confidence, when there are multiple spliced events between the same pair of genes
			if (fusion->second.confidence < CONFIDENCE_HIGH &&
			    fusion->second.spliced1 && fusion->second.spliced2 && !fusion->second.is_read_through() && fusion->second.gene1 != fusion->second.gene2) {
				unsigned int number_of_spliced_breakpoints = 0;
				auto fusions_of_gene = fusions_by_gene.find(fusion->second.gene1);
				for (auto fusion_of_gene = fusions_of_gene->second.begin(); fusion_of_gene != fusions_of_gene->second.end(); ++fusion_of_gene) {
					if ((**fusion_of_gene).gene1 == fusion->second.gene1 && (**fusion_of_gene).gene2 == fusion->second.gene2 &&
					    (**fusion_of_gene).spliced1 && (**fusion_of_gene).spliced2 &&
					    (abs((**fusion_of_gene).breakpoint1 - fusion->second.breakpoint1) > 2 || abs((**fusion_of_gene).breakpoint2 - fusion->second.breakpoint2) > 2))
						++number_of_spliced_breakpoints;
				}
				fusions_of_gene = fusions_by_gene.find(fusion->second.gene2);
				for (auto fusion_of_gene = fusions_of_gene->second.begin(); fusion_of_gene != fusions_of_gene->second.end(); ++fusion_of_gene) {
					if ((**fusion_of_gene).gene1 == fusion->second.gene1 && (**fusion_of_gene).gene2 == fusion->second.gene2 &&
					    (**fusion_of_gene).spliced1 && (**fusion_of_gene).spliced2 &&
					    (abs((**fusion_of_gene).breakpoint1 - fusion->second.breakpoint1) > 2 || abs((**fusion_of_gene).breakpoint2 - fusion->second.breakpoint2) > 2))
						++number_of_spliced_breakpoints;
				}
				if (number_of_spliced_breakpoints > 0)
					fusion->second.confidence++;
			}

			// true events are likely to have at least one spliced breakpoint => penalize if not
			// unless the event is intragenic, because then non-spliced breakpoints are actually more credible, because it's not a circRNA
			if (fusion->second.gene1 != fusion->second.gene2)
				if (fusion->second.confidence > CONFIDENCE_LOW)
					if (!fusion->second.spliced1 && !fusion->second.spliced2)
						fusion->second.confidence--;

			// events with excellent support get high confidence, regardless of whether their breakpoints are spliced
			if (fusion->second.split_reads1 > 20 && fusion->second.split_reads2 > 20 && fusion->second.supporting_reads() > 60)
				fusion->second.confidence = CONFIDENCE_HIGH;

			// decrease the confidence, when something does not look right with the number of supporting reads
			if (fusion->second.confidence > CONFIDENCE_LOW) {

				// there should be supporting reads from both ends of the fusion
				if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0 ||
				    fusion->second.split_reads1 + fusion->second.discordant_mates == 0 ||
				    fusion->second.split_reads2 + fusion->second.discordant_mates == 0) {
					fusion->second.confidence--;

				// the ratio of split reads and discordant mates should be somewhat balanced
				} else if ((fusion->second.split_reads1 + fusion->second.split_reads2) * 20 < fusion->second.discordant_mates) {
					fusion->second.confidence--;

				// the number of supporting reads is not overwhelming compared to the coverage
				} else if (fusion->second.evalue > 0.2 || coverage_fraction < 0.01) {
					fusion->second.confidence = CONFIDENCE_MEDIUM;
				}

			}

			// increase confidence when there is a supporting SV
			if (fusion->second.confidence < CONFIDENCE_HIGH &&
			    fusion->second.closest_genomic_breakpoint1 >= 0 && // has genomic support and
			    (fusion->second.evalue < 0.3 && fusion->second.supporting_reads() >= 2 || // has good e-value or
			     fusion->second.spliced1 && fusion->second.spliced2 && fusion->second.gene1 != fusion->second.gene2 || // was recovered due to splicing or
			     abs(fusion->second.breakpoint1 - fusion->second.closest_genomic_breakpoint1) + abs(fusion->second.breakpoint2 - fusion->second.closest_genomic_breakpoint2) < 20000 || // genomic breakpoints are very close to transcriptomic breakpoints
			     fusion->second.contig1 != fusion->second.contig2 || (abs(fusion->second.breakpoint2 - fusion->second.breakpoint1) > 1000000 && fusion->second.gene1 != fusion->second.gene2))) // distant translocation
				fusion->second.confidence++;

		}
	}
}

// filter speculative fusions without support from WGS
unsigned int filter_no_genomic_support(fusions_t& fusions) {

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter == FILTER_none) {
			if (fusion->second.closest_genomic_breakpoint1 < 0 && // no genomic support
			     fusion->second.confidence == CONFIDENCE_LOW)
				fusion->second.filter = FILTER_no_genomic_support;
			else
				remaining++;
		}
	}

	return remaining;
}

unsigned int recover_genomic_support(fusions_t& fusions) {

	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter == FILTER_none) {
			remaining++;
			continue; // no need to recover fusions that were not filtered
		}

		if (fusion->second.closest_genomic_breakpoint1 >= 0 && // fusion has genomic support
		    (fusion->second.filter == FILTER_end_to_end ||
		     fusion->second.filter == FILTER_intronic ||
		     fusion->second.filter == FILTER_mismappers ||
		     fusion->second.filter == FILTER_no_coverage ||
		     fusion->second.filter == FILTER_in_vitro ||
		     fusion->second.filter == FILTER_relative_support)) {
			fusion->second.filter = FILTER_none;
			remaining++;
		}
	}

	return remaining;
}

