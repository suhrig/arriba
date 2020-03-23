#include <algorithm>
#include <cctype>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "assembly.hpp"
#include "output_fusions.hpp"
#include "read_stats.hpp"

using namespace std;

typedef map< position_t, map<string/*base*/,unsigned int/*frequency*/> > pileup_t;

void pileup_chimeric_alignments(vector<chimeric_alignments_t::iterator>& chimeric_alignments, const unsigned int mate, const bool reverse_complement, const direction_t direction, const position_t breakpoint, pileup_t& pileup) {

	unordered_map< tuple<position_t,position_t>/*intron boundaries*/, unsigned int/*frequency*/> introns;

	for (auto chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if ((**chimeric_alignment).second.filter == FILTER_duplicates)
			continue; // skip duplicates

		alignment_t& read = (**chimeric_alignment).second[mate]; // introduce alias for cleaner code

		if ((**chimeric_alignment).second.size() == 2) // discordant mate
			if (!(direction == DOWNSTREAM && read.strand == FORWARD && read.end   <= breakpoint+2 && read.end   >= breakpoint-200 ||
			      direction == UPSTREAM   && read.strand == REVERSE && read.start >= breakpoint-2 && read.start <= breakpoint+200)) // only consider discordant mates close to the breakpoints (we don't care about the ones in other exons)
				continue;

		string read_sequence = (mate == SUPPLEMENTARY) ? (**chimeric_alignment).second[SPLIT_READ].sequence : read.sequence;
		if (reverse_complement)
			read_sequence = dna_to_reverse_complement(read_sequence);

		position_t read_offset = 0;
		position_t reference_offset = read.start;
		int subtract_from_next_element = 0;
		position_t intron_start;
		for (unsigned int cigar_element = 0; cigar_element < read.cigar.size(); cigar_element++) {
			switch (read.cigar.operation(cigar_element)) {
				case BAM_CINS:
					pileup[reference_offset][read_sequence.substr(read_offset, read.cigar.op_length(cigar_element)+1)]++;
					read_offset += read.cigar.op_length(cigar_element) + 1; // +1, because we take one base from the next element
					++reference_offset; // +1, because we take one base from the next element
					subtract_from_next_element = 1; // because we took one base from the next element
					break;
				case BAM_CREF_SKIP:
					intron_start = reference_offset;
					reference_offset += read.cigar.op_length(cigar_element) - subtract_from_next_element;
					introns[make_tuple(intron_start, reference_offset-1)]++;
					subtract_from_next_element = 0;
					break;
				case BAM_CDEL:
					for (position_t base = 0; base < (int) read.cigar.op_length(cigar_element) - subtract_from_next_element; ++base, ++reference_offset)
						pileup[reference_offset]["-"]++; // indicate deletion by dash
					subtract_from_next_element = 0;
					break;
				case BAM_CSOFT_CLIP:
					if ((**chimeric_alignment).second.size() == 3 && mate == SPLIT_READ &&
					    (cigar_element == 0 && read.strand == FORWARD || cigar_element == read.cigar.size()-1 && read.strand == REVERSE)) {
						if (cigar_element == 0 && read.strand == FORWARD)
							reference_offset -= read.cigar.op_length(cigar_element);
						// fall through to next branch (we want the clipped segment to be part of the pileup to look for non-template bases
					} else {
						read_offset += read.cigar.op_length(cigar_element) - subtract_from_next_element;
						break;
					}
				case BAM_CMATCH:
				case BAM_CEQUAL:
				case BAM_CDIFF:
					for (position_t base = 0; base < (int) read.cigar.op_length(cigar_element) - subtract_from_next_element; ++base, ++read_offset, ++reference_offset)
						pileup[reference_offset][read_sequence.substr(read_offset,1)]++;
					subtract_from_next_element = 0;
					break;
			}
		}
	}

	// mark boundaries of introns
	for (auto intron = introns.begin(); intron != introns.end(); ++intron) {
		position_t intron_start = get<0>(intron->first);
		position_t intron_end = get<1>(intron->first);
		pileup[intron_start][">"] += intron->second; // intron start is represented as ">"
		pileup[intron_end]["<"] += intron->second; // intron end is represented as "<"
		for (auto i = intron_start+1; i < intron_end; ++i)
			pileup[i]["_"]+= intron->second; // intron is represented as "_"
	}
}

void get_sequence_from_pileup(const pileup_t& pileup, const position_t breakpoint, const direction_t direction, const gene_t gene, const assembly_t& assembly, string& sequence, vector<position_t>& positions, string& clipped_sequence) {

	// for each position, find the most frequent allele in the pileup
	bool intron_open = false; // keep track of whether the current position is in an intron
	bool intron_closed = true; // keep track of whether the current position is in an intron
	for (pileup_t::const_iterator position = pileup.begin(); position != pileup.end(); ++position) {

		if (position != pileup.begin() && prev(position)->first < position->first - 1 && !intron_open) {
			sequence += "..."; // indicate uncovered stretches with an ellipsis
			positions.resize(positions.size() + 3, -1);
		}

		// find out base in reference to mark SNPs/SNVs
		string reference_base = "N";
		assembly_t::const_iterator contig_sequence = assembly.find(gene->contig);
		if (contig_sequence != assembly.end())
			reference_base = contig_sequence->second[position->first];

		// find most frequent allele at current position and compute coverage
		auto most_frequent_base = position->second.end();
		unsigned int coverage = 0;
		for (auto base = position->second.begin(); base != position->second.end(); ++base) {
			bool base_is_intron = base->first == "_" || base->first == ">" || base->first == "<";
			if (most_frequent_base == position->second.end() ||
			    base->second > most_frequent_base->second ||
			    (base->second == most_frequent_base->second &&
			     ((base->first == reference_base && most_frequent_base->first != "_" && most_frequent_base->first != ">" && most_frequent_base->first != "<") ||
			      (base->first == "<" && most_frequent_base->first != "_" && most_frequent_base->first != ">") ||
			      (base->first == "_" || base->first == ">"))))
				most_frequent_base = base;
			if (!base_is_intron)
				coverage += base->second;
		}

		// we trust the base, if it has a frequency of >= 75%
		// or if it is an intron with a frequency higher than the coverage
		// or if the base matches the reference
		string most_frequent_base2 = ((most_frequent_base->first == "_" || most_frequent_base->first == ">" || most_frequent_base->first == "<") && most_frequent_base->second >= coverage ||
		                              most_frequent_base->second >= 0.75 * coverage ||
		                              most_frequent_base->first == reference_base) ? most_frequent_base->first : "?";


		if (most_frequent_base2 == "_") { // in the middle of an intron

			if (!intron_open) { // we are in an intron without ever having seen the start of it
				sequence += "...___"; // indicate missing stretch up until exon boundary
				positions.resize(positions.size() + 6, -1);
				intron_open = true;
				intron_closed = false;
			}

		} else if (most_frequent_base2 == ">") { // start of intron

			if (!intron_open) {
				sequence += "___";
				positions.resize(positions.size() + 3, -1);
				intron_open = true;
				intron_closed = false;
			}

		} else if (most_frequent_base2 == "<") { // end of intron

			if (!intron_open) { // we are in an intron without ever having seen the start of it
				sequence += "...___"; // indicate missing stretch up until exon boundary
				positions.resize(positions.size() + 6, -1);
			}
			intron_open = true;
			intron_closed = true;

		} else { // not an intron

			if (!intron_closed) { // we are not in an intron anymore without ever having seen the end of it
				sequence += "..."; // indicate missing stretch up until exon boundary
				positions.resize(positions.size() + 3, -1);
			}
			intron_open = false;
			intron_closed = true;

			// mark SNPs/SNVs via lowercase letters
			if (most_frequent_base2.size() > 1 /*insertion*/ || most_frequent_base2 != reference_base && reference_base != "N")
				std::transform(most_frequent_base2.begin(), most_frequent_base2.end(), most_frequent_base2.begin(), (int (*)(int))std::tolower);

			// mark insertions via brackets
			if (most_frequent_base2.size() > 1) {
				most_frequent_base2 = "[" + most_frequent_base2.substr(0, most_frequent_base2.size()-1) + "]" + most_frequent_base2[most_frequent_base2.size()-1];
				positions.resize(positions.size() + most_frequent_base2.size() - 1, -1);
				if (toupper(most_frequent_base2[most_frequent_base2.size()-1]) == reference_base[0])
					most_frequent_base2[most_frequent_base2.size()-1] = toupper(most_frequent_base2[most_frequent_base2.size()-1]);
			}

			if (direction == UPSTREAM && position->first < breakpoint || direction == DOWNSTREAM && position->first > breakpoint) {
				clipped_sequence += most_frequent_base2;
			} else {
				sequence += most_frequent_base2;
				positions.push_back(position->first);
			}

		}

	}
}

void get_fusion_transcript_sequence(fusion_t& fusion, const assembly_t& assembly, string& sequence, vector<position_t>& positions) {

	if (fusion.predicted_strands_ambiguous || fusion.transcript_start_ambiguous) {
		sequence = "."; // sequence is unknown, because the strands cannot be determined
		positions.push_back(-1);
		return;
	}

	// get the sequences next to the breakpoints
	pileup_t pileup1, pileup2;
	pileup_chimeric_alignments(fusion.split_read1_list, SPLIT_READ, false, fusion.direction1, fusion.breakpoint1, pileup1);
	pileup_chimeric_alignments(fusion.split_read1_list, MATE1, false, fusion.direction1, fusion.breakpoint1, pileup1);
	pileup_chimeric_alignments(fusion.split_read1_list, SUPPLEMENTARY, fusion.direction1 == fusion.direction2, fusion.direction2, fusion.breakpoint2, pileup2);
	pileup_chimeric_alignments(fusion.split_read2_list, SPLIT_READ, false, fusion.direction2, fusion.breakpoint2, pileup2);
	pileup_chimeric_alignments(fusion.split_read2_list, MATE1, false, fusion.direction2, fusion.breakpoint2, pileup2);
	pileup_chimeric_alignments(fusion.split_read2_list, SUPPLEMENTARY, fusion.direction1 == fusion.direction2, fusion.direction1, fusion.breakpoint1, pileup1);
	pileup_chimeric_alignments(fusion.discordant_mate_list, MATE1, false, fusion.direction1, fusion.breakpoint1, pileup1);
	pileup_chimeric_alignments(fusion.discordant_mate_list, MATE2, false, fusion.direction1, fusion.breakpoint1, pileup1);
	pileup_chimeric_alignments(fusion.discordant_mate_list, MATE1, false, fusion.direction2, fusion.breakpoint2, pileup2);
	pileup_chimeric_alignments(fusion.discordant_mate_list, MATE2, false, fusion.direction2, fusion.breakpoint2, pileup2);

	// look for non-template bases inserted between the fused genes
	unsigned int non_template_bases = 0;
	if (!fusion.spliced1 && !fusion.spliced2) {

		map<unsigned int/*number of non-template bases*/, unsigned int/*number of reads with given number of non-template bases*/> non_template_bases_count;
		for (auto read = fusion.split_read1_list.begin(); read != fusion.split_read2_list.end(); ++read) {

			// continue with split_read2_list if we have processed the split_read1_list
			if (read == fusion.split_read1_list.end()) {
				read = fusion.split_read2_list.begin();
				if (read == fusion.split_read2_list.end())
					break;
			}

			// there are non-template bases, if the sum of the clipped bases of split read and supplementary alignment are greater than the read length
			unsigned int clipped_split_read = ((**read).second[SPLIT_READ].strand == FORWARD) ? (**read).second[SPLIT_READ].preclipping() : (**read).second[SPLIT_READ].postclipping();
			unsigned int clipped_supplementary = ((**read).second[SUPPLEMENTARY].strand == FORWARD) ? (**read).second[SUPPLEMENTARY].postclipping() : (**read).second[SUPPLEMENTARY].preclipping();
			if (clipped_split_read + clipped_supplementary >= (**read).second[SPLIT_READ].sequence.size()) {
				unsigned int unmapped_bases = clipped_split_read + clipped_supplementary - (**read).second[SPLIT_READ].sequence.size();
				if (++non_template_bases_count[unmapped_bases] > non_template_bases_count[non_template_bases])
					non_template_bases = unmapped_bases;
			}
		}

	}

	// determine most frequent bases in pileup
	string sequence1, sequence2, clipped_sequence1, clipped_sequence2;
	vector<position_t> positions1, positions2;
	get_sequence_from_pileup(pileup1, fusion.breakpoint1, fusion.direction1, fusion.gene1, assembly, sequence1, positions1, clipped_sequence1);
	get_sequence_from_pileup(pileup2, fusion.breakpoint2, fusion.direction2, fusion.gene2, assembly, sequence2, positions2, clipped_sequence2);

	// if we have no split reads, the exact breakpoints are unknown => use ellipsis to indicate uncertainty
	if (fusion.split_read1_list.size() + fusion.split_read2_list.size() == 0) {
		if (fusion.direction1 == DOWNSTREAM) {
			sequence1 += "...";
			positions1.resize(positions1.size()+3, -1);
		} else { // fusion.direction1 == UPSTREAM
			sequence1 = "..." + sequence1;
			positions1.insert(positions1.begin(), 3, -1);
		}
		if (fusion.direction2 == DOWNSTREAM) {
			sequence2 += "...";
			positions2.resize(positions2.size()+3, -1);
		} else { // fusion.direction2 == UPSTREAM
			sequence2 = "..." + sequence2;
			positions2.insert(positions2.begin(), 3, -1);
		}
	}

	// add non-template bases (if there are any)
	if (non_template_bases > 0) {
		if (clipped_sequence1.size() >= non_template_bases) {
			std::transform(clipped_sequence1.begin(), clipped_sequence1.end(), clipped_sequence1.begin(), (int (*)(int))std::tolower);
			if (fusion.direction1 == UPSTREAM) {
				sequence1 = clipped_sequence1.substr(clipped_sequence1.size() - non_template_bases) + sequence1;
				positions1.insert(positions1.begin(), clipped_sequence1.size() - non_template_bases, -1);
			} else {
				sequence1 += clipped_sequence1.substr(0, non_template_bases);
				positions1.resize(positions1.size()+non_template_bases, -1);
			}
		} else if (clipped_sequence2.size() >= non_template_bases) {
			std::transform(clipped_sequence2.begin(), clipped_sequence2.end(), clipped_sequence2.begin(), (int (*)(int))std::tolower);
			if (fusion.direction2 == UPSTREAM) {
				sequence2 = clipped_sequence2.substr(clipped_sequence2.size() - non_template_bases) + sequence2;
				positions2.insert(positions2.begin(), clipped_sequence2.size() - non_template_bases, -1);
			} else {
				sequence2 += clipped_sequence2.substr(0, non_template_bases);
				positions2.resize(positions2.size()+non_template_bases, -1);
			}
		}
	}

	// look for mismatched bases (i.e., lowercase letters) next to the breakpoints, these are usually non-template bases
	bool sequence1_has_non_template_bases = false;
	bool sequence2_has_non_template_bases = false;
	if (fusion.direction1 == UPSTREAM) {
		int base = 0;
		while (base < (int) sequence1.size() && (sequence1[base] == 'a' || sequence1[base] == 't' || sequence1[base] == 'c' || sequence1[base] == 'g'))
			++base;
		if (base > 0 && base < (int) sequence1.size()) {
			sequence1 = sequence1.substr(0, base) + "|" + sequence1.substr(base);
			fill(positions1.begin(), positions1.begin()+base, -1); // mark non-reference bases
			positions1.insert(positions1.begin()+base, -1); // add position for control character
			sequence1_has_non_template_bases = true;
		}
	} else if (fusion.direction1 == DOWNSTREAM) {
		int base = sequence1.size()-1;
		while (base >= 0 && (sequence1[base] == 'a' || sequence1[base] == 't' || sequence1[base] == 'c' || sequence1[base] == 'g'))
			--base;
		if (base+1 < (int) sequence1.size() && base >= 0) {
			sequence1 = sequence1.substr(0, base+1) + "|" + sequence1.substr(base+1);
			fill(positions1.begin()+base+1, positions1.end(), -1); // mark non-reference bases
			positions1.insert(positions1.begin()+base+1, -1); // add position for control character
			sequence1_has_non_template_bases = true;
		}
	}
	if (fusion.direction2 == UPSTREAM) {
		int base = 0;
		while (base < (int) sequence2.size() && (sequence2[base] == 'a' || sequence2[base] == 't' || sequence2[base] == 'c' || sequence2[base] == 'g'))
			++base;
		if (base > 0 && base < (int) sequence2.size()) {
			sequence2 = sequence2.substr(0, base) + "|" + sequence2.substr(base);
			fill(positions2.begin(), positions2.begin()+base, -1); // mark non-reference bases
			positions2.insert(positions2.begin()+base, -1); // add position for control character
			sequence2_has_non_template_bases = true;
		}
	} else if (fusion.direction2 == DOWNSTREAM) {
		int base = sequence2.size()-1;
		while (base >= 0 && (sequence2[base] == 'a' || sequence2[base] == 't' || sequence2[base] == 'c' || sequence2[base] == 'g'))
			--base;
		if (base+1 < (int) sequence2.size() && base >= 0) {
			sequence2 = sequence2.substr(0, base+1) + "|" + sequence2.substr(base+1);
			fill(positions2.begin()+base+1, positions2.end(), -1); // mark non-reference bases
			positions2.insert(positions2.begin()+base+1, -1); // add position for control character
			sequence2_has_non_template_bases = true;
		}
	}

	if (fusion.transcript_start == TRANSCRIPT_START_GENE1) {
		sequence = ((fusion.predicted_strand1 == FORWARD) ? sequence1 : dna_to_reverse_complement(sequence1));
		if (fusion.predicted_strand1 != FORWARD)
			reverse(positions1.begin(), positions1.end());
		positions = positions1;
		if (!sequence1_has_non_template_bases || !sequence2_has_non_template_bases) { // otherwise the concatenated sequence might have three pipes
			sequence += "|";
			positions.push_back(-1);
		}
		sequence += ((fusion.direction2 == UPSTREAM) ? sequence2 : dna_to_reverse_complement(sequence2));
		if (fusion.direction2 != UPSTREAM)
			reverse(positions2.begin(), positions2.end());
		positions.insert(positions.end(), positions2.begin(), positions2.end());
	} else { // fusion.transcript_start == TRANSCRIPT_START_GENE2)
		sequence = ((fusion.predicted_strand2 == FORWARD) ? sequence2 : dna_to_reverse_complement(sequence2));
		if (fusion.predicted_strand2 != FORWARD)
			reverse(positions2.begin(), positions2.end());
		positions = positions2;
		if (!sequence2_has_non_template_bases || !sequence1_has_non_template_bases) { // otherwise the concatenated sequence might have three pipes
			sequence += "|";
			positions.push_back(-1);
		}
		sequence += ((fusion.direction1 == UPSTREAM) ? sequence1 : dna_to_reverse_complement(sequence1));
		if (fusion.direction1 != UPSTREAM)
			reverse(positions1.begin(), positions1.end());
		positions.insert(positions.end(), positions1.begin(), positions1.end());
	}

	// Arriba uses question marks to denote ambiguous bases;
	// sequencers typically use 'N's
	// => convert 'N's to '?'
	for (size_t i = 0; i < sequence.size(); ++i)
		if (sequence[i] == 'n' || sequence[i] == 'N')
			sequence[i] = '?';
}

bool sort_fusions_by_support(const fusion_t* x, const fusion_t* y) {
	if (x->confidence != y->confidence)
		return x->confidence > y->confidence;
	else if (x->supporting_reads() != y->supporting_reads())
		return x->supporting_reads() > y->supporting_reads();
	else if (x->evalue != y->evalue)
		return x->evalue < y->evalue;
	else if (x->gene1->id != y->gene1->id)      // the following rules don't really sort, they only ensure deterministic sorting and
		return x->gene1->id < y->gene1->id; // that fusions between the same pair of genes are grouped together if e-value and
	else if (x->gene2->id != y->gene2->id)      // supporting reads are equal
		return x->gene2->id < y->gene2->id;
	else if (x->breakpoint1 != y->breakpoint1)
		return x->breakpoint1 < y->breakpoint1;
	else
		return x->breakpoint2 < y->breakpoint2;
}
// make helper struct which groups events between the same pair of genes together,
// such that events with only few supporting reads are listed near the best event (the one with the most supporting reads)
struct sort_fusions_by_rank_of_best_t {
	unordered_map< tuple<gene_t/*gene1*/,gene_t/*gene2*/>, fusion_t* >* best;
	bool operator()(const fusion_t* x, const fusion_t* y) {
		fusion_t* best_x = best->at(make_tuple(x->gene1, x->gene2));
		fusion_t* best_y = best->at(make_tuple(y->gene1, y->gene2));
		if (best_x != best_y)
			return sort_fusions_by_support(best_x, best_y);
		else
			return sort_fusions_by_support(x, y);
	}
};

string gene_to_name(const gene_t gene, const contig_t contig, const position_t breakpoint, gene_annotation_index_t& gene_annotation_index) {
	// if the gene is not a dummy gene (intergenic region), simply return the name of the gene
	if (!gene->is_dummy) {
		return gene->name;
	} else { // is a dummy gene => find flanking genes and report distances to them

		string result;

		// lookup position in gene annotation index
		gene_contig_annotation_index_t::iterator index_hit2 = gene_annotation_index[contig].lower_bound(breakpoint);
		gene_contig_annotation_index_t::reverse_iterator index_hit1(index_hit2);

		// go upstream until we find a non-dummy gene
		while (index_hit1 != gene_annotation_index[contig].rend() && (index_hit1->second.empty() || (**(index_hit1->second.begin())).is_dummy))
			++index_hit1;

		// append upstream flanking genes with distances to gene name
		if (index_hit1 != gene_annotation_index[contig].rend()) {
			for (gene_set_t::iterator gene = index_hit1->second.begin(); gene != index_hit1->second.end(); gene = upper_bound(index_hit1->second.begin(), index_hit1->second.end(), *gene)) {
				if (!(**gene).is_dummy) {
					if (!result.empty())
						result += ",";
					result += (**gene).name + "(" + to_string(static_cast<long long int>(breakpoint - (**gene).end)) + ")";
				}
			}
		}

		// go downstream until we find a non-dummy gene
		while (index_hit2 != gene_annotation_index[contig].end() && (index_hit2->second.empty() || (**(index_hit2->second.begin())).is_dummy))
			++index_hit2;

		// append downstream flanking genes with distances to gene name
		if (index_hit2 != gene_annotation_index[contig].end()) {
			for (gene_set_t::iterator gene = index_hit2->second.begin(); gene != index_hit2->second.end(); gene = upper_bound(index_hit2->second.begin(), index_hit2->second.end(), *gene)) {
				if (!(**gene).is_dummy) {
					if (!result.empty())
						result += ",";
					result += (**gene).name + "(" + to_string(static_cast<long long int>((**gene).start - breakpoint)) + ")";
				}
			}
		}

		if (result.empty())
			result = ".";

		return result;
	}
}

string get_fusion_type(const fusion_t& fusion) {
	if (fusion.contig1 != fusion.contig2) {
		if (fusion.gene1->is_dummy || fusion.gene2->is_dummy ||
		    fusion.direction1 == fusion.direction2 && fusion.gene1->strand != fusion.gene2->strand ||
		    fusion.direction1 != fusion.direction2 && fusion.gene1->strand == fusion.gene2->strand) {
			return "translocation";
		} else {
			if ((fusion.direction1 == UPSTREAM && fusion.gene1->strand == FORWARD || fusion.direction1 == DOWNSTREAM && fusion.gene1->strand == REVERSE) &&
			    (fusion.direction2 == UPSTREAM && fusion.gene2->strand == FORWARD || fusion.direction2 == DOWNSTREAM && fusion.gene2->strand == REVERSE)) {
				return "translocation/3'-3'"; // tail-to-tail fusion
			} else {
				return "translocation/5'-5'"; // head-to-head fusion
			}
		}
	} else { // fusion.contig1 == fusion.contig2
		if (fusion.direction1 == DOWNSTREAM && fusion.direction2 == UPSTREAM) {
			if (fusion.gene1->is_dummy || fusion.gene2->is_dummy ||
			    (fusion.gene1->strand == fusion.gene2->strand)) {
				if (fusion.is_read_through()) {
					return "deletion/read-through";
				} else {
					return "deletion";
				}
			} else if (fusion.gene1->strand == FORWARD || fusion.gene2->strand == REVERSE) {
				if (fusion.is_read_through()) {
					return "deletion/read-through/5'-5'";
				} else {
					return "deletion/5'-5'";
				}
			} else {
				if (fusion.is_read_through()) {
					return "deletion/read-through/3'-3'";
				} else {
					return "deletion/3'-3'"; // tail-to-tail fusion
				}
			}
		} else if (fusion.direction1 == fusion.direction2) {
			if (fusion.gene1->is_dummy || fusion.gene2->is_dummy ||
			    (fusion.gene1->strand != fusion.gene2->strand)) {
				return "inversion";
			} else {
				if (fusion.direction1 == UPSTREAM && fusion.gene1->strand == REVERSE ||
				    fusion.direction1 == UPSTREAM && fusion.gene1->strand == REVERSE) {
					return "inversion/5'-5'";
				} else {
					return "inversion/3'-3'";
				}
			}
		} else { // fusion.direction1 == UPSTREAM && fusion.direction2 == DOWNSTREAM
			if (fusion.gene1->is_dummy || fusion.gene2->is_dummy ||
			    (fusion.gene1->strand == fusion.gene2->strand)) {
				if (fusion.gene1 == fusion.gene2 && fusion.spliced1 && fusion.spliced2) {
					return "duplication/non-canonical_splicing";
				} else {
					return "duplication";
				}
			} else {
				if (fusion.gene1->strand == REVERSE) {
					return "duplication/5'-5'";
				} else {
					return "duplication/3'-3'";
				}
			}
		}
	}
}

string get_fusion_strand(const strand_t strand, const gene_t gene, const bool predicted_strands_ambiguous) {
	string result;

	// determine strand of gene as per annotation
	if (gene->is_dummy)
		result = ".";
	else
		result = (gene->strand == FORWARD) ? "+" : "-";

	// separator between gene strand and fusion strand
	result += "/";

	// output strand of fusion, if it could be predicted
	if (predicted_strands_ambiguous)
		result += ".";
	else
		result += (strand == FORWARD) ? "+" : "-";

	return result;
}

string get_fusion_site(const gene_t gene, const bool spliced, const bool exonic, const contig_t contig, const position_t breakpoint, const exon_annotation_index_t& exon_annotation_index) {
	string site;
	if (gene->is_dummy) {
		site = "intergenic";
	} else if (spliced) {
		site = "splice-site";
	} else if (exonic) {
		// re-annotate exonic breakpoints
		if (breakpoint < gene->start) {
			if (gene->strand == FORWARD)
				site = "5'UTR";
			else
				site = "3'UTR";
		} else if (breakpoint > gene->end) {
			if (gene->strand == FORWARD)
				site = "3'UTR";
			else
				site = "5'UTR";
		} else {
			exon_set_t exons;
			get_annotation_by_coordinate(contig, breakpoint, breakpoint, exons, exon_annotation_index);
			bool has_overlapping_exon = false;
			bool is_utr = true;
			unsigned int is_3_end = 0;
			unsigned int is_5_end = 0;
			for (exon_set_t::iterator exon = exons.begin(); exon != exons.end(); ++exon) {
				if ((**exon).gene == gene) {
					has_overlapping_exon = true;
					if ((**exon).coding_region_start <= breakpoint && (**exon).coding_region_end >= breakpoint)
						is_utr = false;
					if (is_utr && gene->is_protein_coding) {
						// detect if we are in a 5' or 3' UTR by going upstream/downstream and
						// checking whether we first hit a protein-coding exon or the transcript end
						if ((**exon).coding_region_start != -1 && (**exon).coding_region_start > breakpoint) {
							if (gene->strand == FORWARD)
								++is_5_end;
							else
								++is_3_end;
						} else if ((**exon).coding_region_end != -1 && (**exon).coding_region_end < breakpoint) {
							if (gene->strand == REVERSE)
								++is_5_end;
							else
								++is_3_end;
						} else {
							exon_annotation_record_t* next_exon = (**exon).next_exon;
							while (next_exon != NULL && next_exon->coding_region_start == -1)
								next_exon = next_exon->next_exon;
							exon_annotation_record_t* previous_exon = (**exon).previous_exon;
							while (previous_exon != NULL && previous_exon->coding_region_start == -1)
								previous_exon = previous_exon->previous_exon;
							if (previous_exon != NULL || next_exon != NULL) { // is true, if the transcript contains a coding region
								if ((next_exon == NULL) != /*xor*/ (gene->strand == REVERSE))
									++is_3_end;
								else
									++is_5_end;
							}
						}
					}
				}
			}
			if (!has_overlapping_exon) {
				site = "intron";
			} else if (gene->is_protein_coding) {
				if (is_utr) {
					if (is_3_end > is_5_end) {
						site = "3'UTR";
					} else if (is_3_end < is_5_end) {
						site = "5'UTR";
					} else if (is_3_end + is_5_end == 0) {
						site = "exon";
					} else {
						site = "UTR";
					}
				} else {
					site = "CDS";
				}
			} else {
				site = "exon";
			}
		}
	} else {
		site = "intron";
	}
	return site;
}

transcript_t get_transcript(const string& transcript_sequence, const vector<position_t>& transcribed_bases, const gene_t gene, const strand_t strand, const bool strand_ambiguous, const unsigned char which_end, const exon_annotation_index_t& exon_annotation_index) {

	if (strand_ambiguous || strand != gene->strand)
		return NULL; // in case of anti-sense transcription, it makes no sense to determine the transcript

	// determine start and end of 5' or 3' moiety of transcript sequence (whichever is requested)
	size_t from, to, breakpoint;
	if (which_end == 5) {
		from = 0;
		to = transcript_sequence.find('|');
		if (to == string::npos)
			return NULL;
		while (to > 0 && transcribed_bases[to] == -1) // skip control characters in transcript sequence, such as "..."
			to--;
		if (transcribed_bases[to] == -1)
			return NULL; // we get here when the sequence consists of only control characters (should be impossible, but who knows)
		breakpoint = to;
	} else { // which_end == 3
		from = transcript_sequence.find_last_of('|');
		while (from < transcript_sequence.size() && transcribed_bases[from] == -1) // skip control characters in transcript sequence, such as "..."
			from++;
		if (from >= transcript_sequence.size())
			return NULL; // we get here when the sequence is empty or consists of only control characters (should be impossible, but who knows)
		breakpoint = from;
		to = transcript_sequence.size() - 1;
	}
	if (transcribed_bases[from] > transcribed_bases[to])
		swap(from, to); // <from> must point to the base with the lower genomic coordinate

	// for each transcript, calculate score reflecting how well the transcribed bases match annotated exons
	unordered_map<transcript_t,unsigned int> score;
	unordered_map<transcript_t,unsigned int> peak_score;
	unordered_map<transcript_t,bool> is_coding_at_breakpoint;
	size_t position = from;
	for (auto exon_set = exon_annotation_index[gene->contig].lower_bound(transcribed_bases[from]); exon_set != exon_annotation_index[gene->contig].end() && position >= min(from, to) && position <= max(from, to); ++exon_set) {

		// increase score by 1 for each base in an exon that is transcribed
		position_t last_transcribed_base = transcribed_bases[to];
		while (position >= min(from, to) && position <= max(from, to) && transcribed_bases[position] <= exon_set->first) {
			for (auto exon = exon_set->second.begin(); exon != exon_set->second.end(); ++exon) {
				if ((**exon).gene == gene && transcribed_bases[position] >= (**exon).start && transcribed_bases[position] <= (**exon).end) {
					score[(**exon).transcript]++;
					last_transcribed_base = transcribed_bases[position]; // remember for later to compute the number of exonic bases that are NOT transcribed
					if (position == breakpoint) {
						if (transcribed_bases[position] >= (**exon).coding_region_start && transcribed_bases[position] <= (**exon).coding_region_end)
							is_coding_at_breakpoint[(**exon).transcript] = true;
						if (transcribed_bases[position] == (**exon).start && *exon != (**exon).transcript->first_exon ||
						    transcribed_bases[position] == (**exon).end   && *exon != (**exon).transcript->last_exon)
							score[(**exon).transcript] += 10; // give a bonus when the breakpoint coincides with a splice site
					}
				}
			}
			position += (from <= to) ? +1 : -1;
		}

		for (auto exon = exon_set->second.begin(); exon != exon_set->second.end(); ++exon) {
			if ((**exon).gene == gene) {

				// note down peak score for each transcript
				peak_score[(**exon).transcript] = max(score[(**exon).transcript], peak_score[(**exon).transcript]);

				// decrease score by 1 for each base in an exon that is not transcribed
				position_t exon_start = (exon_set != exon_annotation_index[gene->contig].begin()) ? prev(exon_set)->first : (**exon).start-1;
				unsigned int exon_length = min(exon_set->first, transcribed_bases[to]) - max(last_transcribed_base+1, exon_start) + 1;
				score[(**exon).transcript] -= min(exon_length, score[(**exon).transcript]); // score must not go below 0
			}
		}

	}
	if (peak_score.empty())
		return NULL; // transcribed region does not overlap any annotated transcript

	// with which annotated transcript does the transcribed region have the best overlap?
	// search for the transcript with the best ratio of transcribed length/total length
	transcript_t best_transcript = peak_score.begin()->first;
	for (auto transcript = peak_score.begin(); transcript != peak_score.end(); ++transcript) {
		int transcript_length = transcript->first->last_exon->end - transcript->first->first_exon->start;
		int best_transcript_length = best_transcript->last_exon->end - best_transcript->first_exon->start;
		if (transcript->second > peak_score[best_transcript] ||
		    transcript->second == peak_score[best_transcript] && !is_coding_at_breakpoint[best_transcript] && is_coding_at_breakpoint[transcript->first] ||
		    transcript->second == peak_score[best_transcript] && is_coding_at_breakpoint[best_transcript] == is_coding_at_breakpoint[transcript->first] && transcript_length > best_transcript_length ||
		    transcript->second == peak_score[best_transcript] && is_coding_at_breakpoint[best_transcript] == is_coding_at_breakpoint[transcript->first] && transcript_length == best_transcript_length && transcript->first->id < best_transcript->id) { // IDs as tie breaker ensures deterministic behavior
			best_transcript = transcript->first;
		}
	}
	if (peak_score[best_transcript] == 0)
		return NULL;

	return best_transcript;
}

char dna_to_protein(const string& triplet) {
	string t = triplet;
	std::transform(t.begin(), t.end(), t.begin(), (int (*)(int))std::toupper);
	string d = t.substr(0, 2);
	if (d == "GC") { return 'A'; }
	else if (t == "TGT" || t == "TGC") { return 'C'; }
	else if (t == "GAT" || t == "GAC") { return 'D'; }
	else if (t == "GAA" || t == "GAG") { return 'E'; }
	else if (t == "TTT" || t == "TTC") { return 'F'; }
	else if (d == "GG") { return 'G'; }
	else if (t == "CAT" || t == "CAC") { return 'H'; }
	else if (t == "ATT" || t == "ATC" || t == "ATA") { return 'I'; }
	else if (t == "AAA" || t == "AAG") { return 'K'; }
	else if (d == "CT" || t == "TTA" || t == "TTG") { return 'L'; }
	else if (t == "ATG") { return 'M'; }
	else if (t == "AAT" || t == "AAC") { return 'N'; }
	else if (d == "CC") { return 'P'; }
	else if (t == "CAA" || t == "CAG") { return 'Q'; }
	else if (d == "CG" || t == "AGA" || t == "AGG") { return 'R'; }
	else if (d == "TC" || t == "AGT" || t == "AGC") { return 'S'; }
	else if (d == "AC") { return 'T'; }
	else if (d == "GT") { return 'V'; }
	else if (t == "TGG") { return 'W'; }
	else if (t == "TAT" || t == "TAC") { return 'Y'; }
	else if (t == "TAA" || t == "TAG" || t == "TGA") { return '*'; }
	else { return '?'; }
}

// determines reading frame of first base of given transcript based on coding exons overlapping the transcript
int get_reading_frame(const vector<position_t>& transcribed_bases, const int from, const int to, const transcript_t transcript, const gene_t gene, const assembly_t& assembly, exon_t& exon_with_start_codon) {

	// find exon containing start codon
	if (transcript == NULL)
		exon_with_start_codon = NULL;
	else
		exon_with_start_codon = (gene->strand == FORWARD) ? transcript->first_exon : transcript->last_exon;
	while (exon_with_start_codon != NULL && exon_with_start_codon->coding_region_start == -1)
		exon_with_start_codon = (gene->strand == FORWARD) ? exon_with_start_codon->next_exon : exon_with_start_codon->previous_exon;
	if (exon_with_start_codon == NULL)
		return -1; // no coding exon found => is a non-coding transcript

	// check if there is a start codon at the start of first coding exon, otherwise the annotation is incorrect
	string first_codon;
	if (gene->strand == FORWARD)
		first_codon = assembly.at(gene->contig).substr(exon_with_start_codon->coding_region_start, 3);
	else // gene->strand == REVERSE
		first_codon = dna_to_reverse_complement(assembly.at(gene->contig).substr(exon_with_start_codon->coding_region_end - 2, 3));
	if (first_codon != "ATG")
		return -1;

	// find a coding exon that is transcribed to determine reading frame
	int reading_frame = -1;
	position_t transcribed_coding_base = -1;
	for (exon_t exon = exon_with_start_codon; exon != NULL && transcribed_coding_base == -1; exon = (gene->strand == FORWARD) ? exon->next_exon : exon->previous_exon) {
		for (int position = from; position <= to && transcribed_coding_base == -1; position++)
			if (exon->coding_region_start <= transcribed_bases[position] && exon->coding_region_end >= transcribed_bases[position])
				transcribed_coding_base = position;
		if (transcribed_coding_base == -1) {
			reading_frame = (reading_frame + exon->coding_region_end - exon->coding_region_start + 1) % 3;
		} else {
			if (gene->strand == FORWARD)
				reading_frame += transcribed_bases[transcribed_coding_base] - exon->coding_region_start;
			else
				reading_frame += exon->coding_region_end - transcribed_bases[transcribed_coding_base];
			reading_frame = (reading_frame + 1) % 3;
		}
	}
	if (transcribed_coding_base == -1) // fusion transcript does not overlap any coding region
		return -1;

	// compute reading frame at start of fusion transcript
	for (int position = transcribed_coding_base - 1; position >= from; --position)
		if (transcribed_bases[position] != -1) // skip control characters and insertions
			reading_frame = (reading_frame == 0) ? 2 : reading_frame-1;

	return reading_frame;
}

string get_fusion_peptide_sequence(const string& transcript_sequence, const vector<position_t>& positions, const gene_t gene_5, const gene_t gene_3, const transcript_t transcript_5, const transcript_t transcript_3, const strand_t predicted_strand_3, const exon_annotation_index_t& exon_annotation_index, const assembly_t& assembly) {

	// abort if there is uncertainty in the transcript sequence around the junction or if the transcript sequence is unknown
	if (transcript_sequence.empty() || transcript_sequence == "." ||
	    transcript_sequence.find("...|") != string::npos || transcript_sequence.find("|...") != string::npos)
		return ".";

	if (assembly.find(gene_5->contig) == assembly.end() || assembly.find(gene_3->contig) == assembly.end())
		return "."; // we need the assembly to search for the start codon

	// split transcript into 5' and 3' parts and (possibly) non-template bases
	// moreover, remove sequences beyond "...", i.e., regions with unclear sequence
	size_t transcription_5_end = transcript_sequence.find('|') - 1;
	size_t transcription_5_start = transcript_sequence.rfind("...", transcription_5_end);
	if (transcription_5_start == string::npos)
		transcription_5_start = 0;
	else
		transcription_5_start += 3; // skip "..."
	size_t non_template_bases_length = transcript_sequence.find('|', transcription_5_end + 2);
	if (non_template_bases_length == string::npos)
		non_template_bases_length = 0;
	else
		non_template_bases_length -= transcription_5_end + 2;
	size_t transcription_3_start = transcription_5_end + 2;
	if (non_template_bases_length > 0)
		transcription_3_start += non_template_bases_length + 1;
	size_t transcription_3_end = transcript_sequence.find("...", transcription_3_start);
	if (transcription_3_end == string::npos)
		transcription_3_end = transcript_sequence.size() - 1;
	else
		transcription_3_end--;

	// determine reading frame of 5' gene
	exon_t start_exon_5 = NULL;
	int reading_frame_5 = get_reading_frame(positions, transcription_5_start, transcription_5_end, transcript_5, gene_5, assembly, start_exon_5);
	if (reading_frame_5 == -1)
		return "."; // 5' gene has no coding exons overlapping the transcribed region
	else if (reading_frame_5 != 0)
		reading_frame_5 = 3 - reading_frame_5;

	// determine reading frame of 3' gene
	exon_t start_exon_3 = NULL;
	int reading_frame_3 = -1;
	if (gene_3->strand == predicted_strand_3) // it makes no sense to determine the reading frame in case of anti-sense transcription
		reading_frame_3 = get_reading_frame(positions, transcription_3_start, transcription_3_end, transcript_3, gene_3, assembly, start_exon_3);

	// translate DNA to protein
	string peptide_sequence;
	peptide_sequence.reserve(transcript_sequence.size()/3+2);
	bool inside_indel = false; // needed to keep track of frame shifts (=> if so, convert to lowercase)
	bool inside_intron = false; // keep track of whether one or more bases of the current codon are in an intron (=> convert to lowercase)
	int frame_shift = 0; // keeps track of whether the current reading frame is shifted (=> convert to lowercase)
	int codon_5_bases = 0; // keeps track of how many bases of a codon come from the 5' gene to check if the codon spans the breakpoint
	int codon_3_bases = 0; // keeps track of how many bases of a codon come from the 3' gene to check if the codon spans the breakpoint
	bool found_start_codon = false;
	string codon;
	string reference_codon;
	for (size_t position = transcription_5_start + reading_frame_5; position < transcription_3_end; ++position) {

		// don't begin before start codon
		if (!found_start_codon) {
			if (positions[position] != -1 &&
			    (gene_5->strand == FORWARD && positions[position] >= start_exon_5->coding_region_start ||
			     gene_5->strand == REVERSE && positions[position] <= start_exon_5->coding_region_end))
				found_start_codon = true;
			else
				continue;
		}

		if (transcript_sequence[position] == 'A' || transcript_sequence[position] == 'T' || transcript_sequence[position] == 'C' || transcript_sequence[position] == 'G' ||
		    transcript_sequence[position] == 'a' || transcript_sequence[position] == 't' || transcript_sequence[position] == 'c' || transcript_sequence[position] == 'g' ||
		    transcript_sequence[position] == '?') {

			// count how many bases of the codon come from the 5' gene and how many from the 3' gene
			// to determine if a codon overlaps the breakpoint
			if (codon.size() == 0) {
				codon_5_bases = 0;
				codon_3_bases = 0;
			}
			if (position <= transcription_5_end)
				codon_5_bases++;
			else if (position >= transcription_3_start)
				codon_3_bases++;

			codon += transcript_sequence[position];

			// compare reference base at given position to check for non-silent SNPs/somatic SNVs
			if (positions[position] != -1) { // is not a control character or an insertion
				reference_codon += assembly.at((position <= transcription_5_end) ? gene_5->contig : gene_3->contig)[positions[position]];
				if (position <= transcription_5_end && gene_5->strand == REVERSE || position >= transcription_3_start && gene_3->strand == REVERSE)
					reference_codon[reference_codon.size()-1] = dna_to_complement(reference_codon[reference_codon.size()-1]);
			}

			if (inside_indel)
				frame_shift = (frame_shift + 1) % 3;

		} else if (transcript_sequence[position] == '[') {
			inside_indel = true;
		} else if (transcript_sequence[position] == ']') {
			inside_indel = false;
		} else if (transcript_sequence[position] == '-') {
			frame_shift = (frame_shift == 0) ? 2 : frame_shift - 1;
		}

		// check if we are in an intron/intergenic region or an exon
		if (positions[position] != -1) {
			inside_intron = true;
			exon_t next_exon = (position <= transcription_5_end) ? start_exon_5 : start_exon_3;
			while (next_exon != NULL) {
				if (positions[position] >= next_exon->start && positions[position] <= next_exon->end) {
					if (positions[position] >= next_exon->coding_region_start && positions[position] <= next_exon->coding_region_end)
						inside_intron = false;
					break;
				}
				if (position <= transcription_5_end && gene_5->strand == FORWARD || position >= transcription_3_start && gene_3->strand == FORWARD)
					next_exon = next_exon->next_exon;
				else
					next_exon = next_exon->previous_exon;
			}
		}

		if (codon.size() == 3) { // codon completed => translate to amino acid

			// translate codon to amino acid and check whether it differs from the reference assembly
			char amino_acid = dna_to_protein(codon);
			char reference_amino_acid = dna_to_protein(reference_codon);

			// convert aberrant amino acids to lowercase
			if (position > transcription_5_end && position < transcription_3_start || // non-template base
			    amino_acid != reference_amino_acid || // non-silent mutation
			    inside_indel || // indel
			    inside_intron || // intron
			    codon_5_bases != 3 && position <= transcription_5_end || // codon overlaps 5' breakpoint
			    codon_3_bases != 3 && position >= transcription_3_start || // codon overlaps 3' breakpoint
			    frame_shift != 0 || // out-of-frame base
			    position >= transcription_3_start && reading_frame_3 == -1) // 3' end is not a coding region
				amino_acid = tolower(amino_acid);

			peptide_sequence += amino_acid;
			codon.clear();
			reference_codon.clear();

			// terminate, if we hit a stop codon in the 3' gene
			if (codon_3_bases >= 2 && amino_acid == '*')
				break;
		}


		// mark end of 5' end as pipe
		if (position == transcription_5_end && codon.size() <= 1 ||
		    codon_5_bases == 2 && codon.size() == 0)
			if (peptide_sequence.empty() || peptide_sequence[peptide_sequence.size()-1] != '|')
				peptide_sequence += '|';

		// mark beginning of 3' end as pipe
		if (non_template_bases_length > 0)
			if (position + 2 == transcription_3_start && codon.size() <= 1 ||
			    codon_3_bases == 1 && codon.size() == 0)
				if (peptide_sequence.empty() || peptide_sequence[peptide_sequence.size()-1] != '|')
					peptide_sequence += '|';

		// check if fusion shifts frame of 3' gene
		if (position == transcription_3_start && reading_frame_3 != -1)
			frame_shift = (3 + reading_frame_3 + 1 - codon.size()) % 3;
	}

	return peptide_sequence;
}

string is_in_frame(const string& fusion_peptide_sequence) {

	if (fusion_peptide_sequence == ".")
		return ".";

	// declare fusion as out-of-frame, if there is a stop codon before the junction
	// unless there are in-frame codons between this stop codon and the junction
	size_t fusion_junction = fusion_peptide_sequence.rfind('|');
	size_t stop_codon_before_junction = fusion_peptide_sequence.rfind('*', fusion_junction);
	bool in_frame_codons_before_junction = false;
	for (size_t amino_acid = (stop_codon_before_junction < fusion_junction) ? stop_codon_before_junction+1 : 0; amino_acid < fusion_junction; ++amino_acid)
		if (fusion_peptide_sequence[amino_acid] >= 'A' && fusion_peptide_sequence[amino_acid] <= 'Z') {
			in_frame_codons_before_junction = true;
			break;
		}
	if (!in_frame_codons_before_junction)
		if (stop_codon_before_junction < fusion_junction)
			return "stop-codon";
		else
			return "out-of-frame";

	// a fusion is in-frame, if there is at least one amino acid in the 3' end of the fusion
	// that matches an amino acid of the 3' gene
	// such amino acids can easily be identified, because they are uppercase in the fusion sequence
	for (size_t amino_acid = fusion_peptide_sequence.size() - 1; amino_acid > fusion_junction; --amino_acid)
		if (fusion_peptide_sequence[amino_acid] >= 'A' && fusion_peptide_sequence[amino_acid] <= 'Z')
			return "in-frame";

	return "out-of-frame";
}

void write_fusions_to_file(fusions_t& fusions, const string& output_file, const coverage_t& coverage, const assembly_t& assembly, gene_annotation_index_t& gene_annotation_index, exon_annotation_index_t& exon_annotation_index, vector<string> contigs_by_id, const bool print_supporting_reads, const bool print_fusion_sequence, const bool print_peptide_sequence, const bool write_discarded_fusions) {
//TODO add "chr", if necessary

	// make a vector of pointers to all fusions
	// the vector will hold the fusions in sorted order
	vector<fusion_t*> sorted_fusions;
	if (write_discarded_fusions)
		sorted_fusions.reserve(fusions.size());
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (write_discarded_fusions != (fusion->second.filter == FILTER_none)) // either write filtered or unfiltered fusions
			sorted_fusions.push_back(&(fusion->second));

	// don't sort the discarded fusions
	if (!write_discarded_fusions) {

		// sometimes we report multiple breakpoints for the same pair of genes
		// in such cases, we want to make sure that all breakpoints of the same pair of genes
		// are grouped together in the output file, even if some of them have very low support
		// all such breakpoints will get the rank of the best ranking breakpoints
		// => find out what the best ranking breakpoints are for each pair of genes
		//    we store these breakpoints in a sorting object (sort_fusions_by_rank_of_best),
		//    which we use as a parameter to the sort function later
		unordered_map< tuple<gene_t/*gene1*/,gene_t/*gene2*/>, fusion_t* > best_fusion_by_gene_pair;
		for (auto fusion = sorted_fusions.begin(); fusion != sorted_fusions.end(); ++fusion) {
			fusion_t*& current_best = best_fusion_by_gene_pair[make_tuple((**fusion).gene1, (**fusion).gene2)];
			if (current_best == NULL || sort_fusions_by_support(*fusion, current_best))
				current_best = *fusion;
		}

		// sort all gene pairs by the rank of the best scoring breakpoints of a given gene pair
		sort_fusions_by_rank_of_best_t sort_fusions_by_rank_of_best;
		sort_fusions_by_rank_of_best.best = &best_fusion_by_gene_pair;
		sort(sorted_fusions.begin(), sorted_fusions.end(), sort_fusions_by_rank_of_best);
	}

	// write sorted list to file
	ofstream out(output_file);
	if (!out.is_open()) {
		cerr << "ERROR: failed to open output file: " << output_file << endl;
		exit(1);
	}
	out << "#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\tsite1\tsite2\ttype\tdirection1\tdirection2\tsplit_reads1\tsplit_reads2\tdiscordant_mates\tcoverage1\tcoverage2\tconfidence\tclosest_genomic_breakpoint1\tclosest_genomic_breakpoint2\tfilters\tfusion_transcript\tgene_id1\tgene_id2\treading_frame\ttranscript_id1\ttranscript_id2\tpeptide_sequence\tread_identifiers" << endl;
	for (auto fusion = sorted_fusions.begin(); fusion != sorted_fusions.end(); ++fusion) {

		// describe site of breakpoint
		string site_5 = get_fusion_site((**fusion).gene1, (**fusion).spliced1, (**fusion).exonic1, (**fusion).contig1, (**fusion).breakpoint1, exon_annotation_index);
		string site_3 = get_fusion_site((**fusion).gene2, (**fusion).spliced2, (**fusion).exonic2, (**fusion).contig2, (**fusion).breakpoint2, exon_annotation_index);
		
		// convert closest genomic breakpoints to strings of the format <chr>:<position>(<distance to transcriptomic breakpoint>)
		string closest_genomic_breakpoint_5, closest_genomic_breakpoint_3;
		if ((**fusion).closest_genomic_breakpoint1 >= 0) {
			closest_genomic_breakpoint_5 = contigs_by_id[(**fusion).contig1] + ":" + to_string(static_cast<long long int>((**fusion).closest_genomic_breakpoint1+1)) + "(" + to_string(static_cast<long long int>(abs((**fusion).breakpoint1 - (**fusion).closest_genomic_breakpoint1))) + ")";
		} else {
			closest_genomic_breakpoint_5 = ".";
		}
		if ((**fusion).closest_genomic_breakpoint2 >= 0) {
			closest_genomic_breakpoint_3 = contigs_by_id[(**fusion).contig2] + ":" + to_string(static_cast<long long int>((**fusion).closest_genomic_breakpoint2+1)) + "(" + to_string(static_cast<long long int>(abs((**fusion).breakpoint2 - (**fusion).closest_genomic_breakpoint2))) + ")";
		} else {
			closest_genomic_breakpoint_3 = ".";
		}

		// assign confidence scores
		string confidence;
		switch ((**fusion).confidence) {
			case CONFIDENCE_LOW:
				confidence = "low";
				break;
			case CONFIDENCE_MEDIUM:
				confidence = "medium";
				break;
			case CONFIDENCE_HIGH:
				confidence = "high";
				break;
		}

		// the 5' gene should always come first => swap columns, if necessary
		gene_t gene_5 = (**fusion).gene1; gene_t gene_3 = (**fusion).gene2;
		contig_t contig_5 = (**fusion).contig1; contig_t contig_3 = (**fusion).contig2;
		position_t breakpoint_5 = (**fusion).breakpoint1; position_t breakpoint_3 = (**fusion).breakpoint2;
		direction_t direction_5 = (**fusion).direction1; direction_t direction_3 = (**fusion).direction2;
		strand_t strand_5 = (**fusion).predicted_strand1; strand_t strand_3 = (**fusion).predicted_strand2;
		unsigned int split_reads_5 = (**fusion).split_reads1; unsigned int split_reads_3 = (**fusion).split_reads2;
		if ((**fusion).transcript_start == TRANSCRIPT_START_GENE2) {
			swap(gene_5, gene_3);
			swap(direction_5, direction_3);
			swap(contig_5, contig_3);
			swap(breakpoint_5, breakpoint_3);
			swap(site_5, site_3);
			swap(split_reads_5, split_reads_3);
			swap(closest_genomic_breakpoint_5, closest_genomic_breakpoint_3);
			swap(strand_5, strand_3);
		}

		int coverage_5 = coverage.get_coverage(contig_5, breakpoint_5, (direction_5 == UPSTREAM) ? DOWNSTREAM : UPSTREAM);
		int coverage_3 = coverage.get_coverage(contig_3, breakpoint_3, (direction_3 == UPSTREAM) ? DOWNSTREAM : UPSTREAM);

		// write line to output file
		out << gene_to_name(gene_5, contig_5, breakpoint_5, gene_annotation_index) << "\t" << gene_to_name(gene_3, contig_3, breakpoint_3, gene_annotation_index) << "\t"
		    << get_fusion_strand(strand_5, gene_5, (**fusion).predicted_strands_ambiguous) << "\t" << get_fusion_strand(strand_3, gene_3, (**fusion).predicted_strands_ambiguous) << "\t"
		    << contigs_by_id[contig_5] << ":" << (breakpoint_5+1) << "\t" << contigs_by_id[contig_3] << ":" << (breakpoint_3+1) << "\t"
		    << site_5 << "\t" << site_3 << "\t"
		    << get_fusion_type(**fusion) << "\t" << ((direction_5 == UPSTREAM) ? "upstream" : "downstream") << "\t" << ((direction_3 == UPSTREAM) ? "upstream" : "downstream") << "\t"
		    << split_reads_5 << "\t" << split_reads_3 << "\t" << (**fusion).discordant_mates << "\t"
		    << ((coverage_5 >= 0) ? to_string(static_cast<long long int>(coverage_5)) : ".") << "\t" << ((coverage_3 >= 0) ? to_string(static_cast<long long int>(coverage_3)) : ".") << "\t"
		    << confidence << "\t"
		    << closest_genomic_breakpoint_5 << "\t" << closest_genomic_breakpoint_3;

		// count the number of reads discarded by a given filter
		map<string,unsigned int> filters;
		if ((**fusion).filter != FILTER_none)
			filters[FILTERS[(**fusion).filter]] = 0;
		vector<chimeric_alignments_t::iterator> all_supporting_reads;
		all_supporting_reads.insert(all_supporting_reads.end(), (**fusion).split_read1_list.begin(), (**fusion).split_read1_list.end());
		all_supporting_reads.insert(all_supporting_reads.end(), (**fusion).split_read2_list.begin(), (**fusion).split_read2_list.end());
		all_supporting_reads.insert(all_supporting_reads.end(), (**fusion).discordant_mate_list.begin(), (**fusion).discordant_mate_list.end());
		for (auto chimeric_alignment = all_supporting_reads.begin(); chimeric_alignment != all_supporting_reads.end(); ++chimeric_alignment)
			if ((**chimeric_alignment).second.filter != FILTER_none)
				filters[FILTERS[(**chimeric_alignment).second.filter]]++;

		// output filters
		out << "\t";
		if (filters.empty()) {
			out << ".";
		} else {
			for (auto filter = filters.begin(); filter != filters.end(); ++filter) {
				if (filter != filters.begin())
					out << ",";
				out << filter->first;
				if (filter->second != 0)
					out << "(" << filter->second << ")";
			}
		}

		// print a fusion-spanning sequence
		out << "\t";
		string transcript_sequence;
		vector<position_t> positions;
		if (print_fusion_sequence || print_peptide_sequence)
			get_fusion_transcript_sequence(**fusion, assembly, transcript_sequence, positions);
		if (print_fusion_sequence) {
			out << transcript_sequence;
		} else {
			out << ".";
		}

		// print gene IDs
		out << "\t" << ((gene_5->is_dummy) ? "." : gene_5->gene_id) << "\t" << ((gene_3->is_dummy) ? "." : gene_3->gene_id);

		// print the translated protein sequence
		out << "\t";
		if (print_peptide_sequence) {
			transcript_t transcript_5 = get_transcript(transcript_sequence, positions, gene_5, strand_5, (**fusion).predicted_strands_ambiguous, 5, exon_annotation_index);
			transcript_t transcript_3 = get_transcript(transcript_sequence, positions, gene_3, strand_3, (**fusion).predicted_strands_ambiguous, 3, exon_annotation_index);
			string fusion_peptide_sequence = get_fusion_peptide_sequence(transcript_sequence, positions, gene_5, gene_3, transcript_5, transcript_3, strand_3, exon_annotation_index, assembly);
			out << is_in_frame(fusion_peptide_sequence) << "\t"
			    << ((transcript_5 == NULL) ? "." : transcript_5->name) << "\t"
			    << ((transcript_3 == NULL) ? "." : transcript_3->name) << "\t"
			    << fusion_peptide_sequence;
		} else {
			out << ".\t.\t.\t.";
		}

		// if requested, print identifiers of supporting reads
		out << "\t";
		if (print_supporting_reads && !all_supporting_reads.empty()) {
			for (auto read = all_supporting_reads.begin(); read != all_supporting_reads.end(); ++read) {
				if (read != all_supporting_reads.begin())
					out << ",";
				out << strip_hi_tag_from_read_name((**read).first);
			}
		} else {
			out << ".";
		}

		out << endl;
	}
	out.close();
	if (out.bad()) {
		cerr << "ERROR: failed to write to file" << endl;
		exit(1);
	}
}

