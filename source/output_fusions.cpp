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

using namespace std;

typedef map< position_t, map<string/*base*/,unsigned int/*frequency*/> > pileup_t;

void pileup_chimeric_alignments(vector<chimeric_alignments_t::iterator>& chimeric_alignments, const unsigned int mate, const bool reverse_complement, const direction_t direction, const position_t breakpoint, pileup_t& pileup) {

	unordered_map< tuple<position_t,position_t>/*intron boundaries*/, unsigned int/*frequency*/> introns;

	for (auto chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if ((**chimeric_alignment).second.filter == FILTERS.at("duplicates"))
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
		for (int cigar_element = 0; cigar_element < read.cigar.size(); cigar_element++) {
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
					for (position_t base = 0; base < read.cigar.op_length(cigar_element) - subtract_from_next_element; ++base, ++reference_offset)
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
					for (position_t base = 0; base < read.cigar.op_length(cigar_element) - subtract_from_next_element; ++base, ++read_offset, ++reference_offset)
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

void get_sequence_from_pileup(const pileup_t& pileup, const position_t breakpoint, const direction_t direction, const gene_t gene, const assembly_t& assembly, string& sequence, string& clipped_sequence) {

	// for each position, find the most frequent allele in the pileup
	bool intron_open = false; // keep track of whether the current position is in an intron
	bool intron_closed = true; // keep track of whether the current position is in an intron
	for (pileup_t::const_iterator position = pileup.begin(); position != pileup.end(); ++position) {

		if (position != pileup.begin() && prev(position)->first < position->first - 1 && !intron_open)
			sequence += "..."; // indicate uncovered stretches with an ellipsis

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
		                              most_frequent_base->first == reference_base) ? most_frequent_base->first : "n";


		if (most_frequent_base2 == "_") { // in the middle of an intron

			if (!intron_open) { // we are in an intron without ever having seen the start of it
				sequence += "...___"; // indicate missing stretch up until exon boundary
				intron_open = true;
				intron_closed = false;
			}

		} else if (most_frequent_base2 == ">") { // start of intron

			if (!intron_open) {
				sequence += "___";
				intron_open = true;
				intron_closed = false;
			}

		} else if (most_frequent_base2 == "<") { // end of intron

			if (!intron_open) // we are in an intron without ever having seen the start of it
				sequence += "...___"; // indicate missing stretch up until exon boundary
			intron_open = true;
			intron_closed = true;

		} else { // not an intron

			if (!intron_closed) // we are not in an intron anymore without ever having seen the end of it
				sequence += "..."; // indicate missing stretch up until exon boundary
			intron_open = false;
			intron_closed = true;

			// mark SNPs/SNVs via lowercase letters
			if (most_frequent_base2.size() > 1 /*insertion*/ || most_frequent_base2 != reference_base && reference_base != "N")
				std::transform(most_frequent_base2.begin(), most_frequent_base2.end(), most_frequent_base2.begin(), (int (*)(int))std::tolower);

			// mark insertions via brackets
			if (most_frequent_base2.size() > 1) {
				most_frequent_base2 = "[" + most_frequent_base2.substr(0, most_frequent_base2.size()-1) + "]" + most_frequent_base2[most_frequent_base2.size()-1];
				if (toupper(most_frequent_base2[most_frequent_base2.size()-1]) == reference_base[0])
					most_frequent_base2[most_frequent_base2.size()-1] = toupper(most_frequent_base2[most_frequent_base2.size()-1]);
			}

			if (direction == UPSTREAM && position->first < breakpoint || direction == DOWNSTREAM && position->first >= breakpoint)
				clipped_sequence += most_frequent_base2;
			else
				sequence += most_frequent_base2;

		}

	}
}

string get_fusion_transcript_sequence(fusion_t& fusion, const assembly_t& assembly) {

	if (fusion.predicted_strands_ambiguous || fusion.transcript_start_ambiguous)
		return "."; // sequence is unknown, because the strands cannot be determined

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
	int non_template_bases = 0;
	if (!fusion.spliced1 && !fusion.spliced2) {

		map<int, unsigned int> non_template_bases_count;
		for (auto read = fusion.split_read1_list.begin(); read != fusion.split_read2_list.end(); ++read) {

			// continue with split_read2_list if we have processed the split_read1_list
			if (read == fusion.split_read1_list.end()) {
				read = fusion.split_read2_list.begin();
				if (read == fusion.split_read2_list.end())
					break;
			}

			// there are non-template bases, if the sum of the clipped bases of split read and supplementary alignment are greater than the read length
			int clipped_split_read = ((**read).second[SPLIT_READ].strand == FORWARD) ? (**read).second[SPLIT_READ].cigar.op_length(0) : (**read).second[SPLIT_READ].cigar.op_length((**read).second[SPLIT_READ].cigar.size()-1);
			int clipped_supplementary = ((**read).second[SUPPLEMENTARY].strand == FORWARD) ? (**read).second[SUPPLEMENTARY].cigar.op_length((**read).second[SUPPLEMENTARY].cigar.size()-1) : (**read).second[SUPPLEMENTARY].cigar.op_length(0);
			int unmapped_bases = clipped_split_read + clipped_supplementary - (**read).second[SPLIT_READ].sequence.size();
			if (++non_template_bases_count[unmapped_bases] > non_template_bases_count[non_template_bases])
				non_template_bases = unmapped_bases;
		}

	}

	// determine most frequent bases in pileup
	string sequence1, sequence2, clipped_sequence1, clipped_sequence2;
	get_sequence_from_pileup(pileup1, fusion.breakpoint1, fusion.direction1, fusion.gene1, assembly, sequence1, clipped_sequence1);
	get_sequence_from_pileup(pileup2, fusion.breakpoint2, fusion.direction2, fusion.gene2, assembly, sequence2, clipped_sequence2);

	// if we have no split reads, the exact breakpoints are unknown => use ellipsis to indicate uncertainty
	if (fusion.split_read1_list.size() + fusion.split_read2_list.size() == 0) {
		if (fusion.direction1 == DOWNSTREAM)
			sequence1 += "...";
		else // fusion.direction1 == UPSTREAM
			sequence1 = "..." + sequence1;
		if (fusion.direction2 == DOWNSTREAM)
			sequence2 += "...";
		else // fusion.direction2 == UPSTREAM
			sequence2 = "..." + sequence2;
	}

	// add non-template bases (if there are any)
	if (non_template_bases > 0) {
		if (clipped_sequence1.size() >= non_template_bases) {
			std::transform(clipped_sequence1.begin(), clipped_sequence1.end(), clipped_sequence1.begin(), (int (*)(int))std::tolower);
			if (fusion.direction1 == UPSTREAM)
				sequence1 = clipped_sequence1.substr(clipped_sequence1.size() - non_template_bases) + sequence1;
			else
				sequence1 += clipped_sequence1.substr(0, non_template_bases);
		} else if (clipped_sequence2.size() >= non_template_bases) {
			std::transform(clipped_sequence2.begin(), clipped_sequence2.end(), clipped_sequence2.begin(), (int (*)(int))std::tolower);
			if (fusion.direction2 == UPSTREAM)
				sequence2 = clipped_sequence2.substr(clipped_sequence2.size() - non_template_bases) + sequence2;
			else
				sequence2 += clipped_sequence2.substr(0, non_template_bases);
		}
	}

	// look for mismatched bases (i.e., lowercase letters) next to the breakpoints, these are usually non-template bases
	bool sequence1_has_non_template_bases = false;
	bool sequence2_has_non_template_bases = false;
	if (fusion.direction1 == UPSTREAM) {
		int base = 0;
		while (base < sequence1.size() && (sequence1[base] == 'a' || sequence1[base] == 't' || sequence1[base] == 'c' || sequence1[base] == 'g'))
			++base;
		if (base > 0 && base < sequence1.size()) {
			sequence1 = sequence1.substr(0, base) + "|" + sequence1.substr(base);
			sequence1_has_non_template_bases = true;
		}
	} else if (fusion.direction1 == DOWNSTREAM) {
		int base = sequence1.size()-1;
		while (base >= 0 && (sequence1[base] == 'a' || sequence1[base] == 't' || sequence1[base] == 'c' || sequence1[base] == 'g'))
			--base;
		if (base < sequence1.size()-1 && base >= 0) {
			sequence1 = sequence1.substr(0, base+1) + "|" + sequence1.substr(base+1);
			sequence1_has_non_template_bases = true;
		}
	}
	if (fusion.direction2 == UPSTREAM) {
		int base = 0;
		while (base < sequence2.size() && (sequence2[base] == 'a' || sequence2[base] == 't' || sequence2[base] == 'c' || sequence2[base] == 'g'))
			++base;
		if (base > 0 && base < sequence2.size()) {
			sequence2 = sequence2.substr(0, base) + "|" + sequence2.substr(base);
			sequence2_has_non_template_bases = true;
		}
	} else if (fusion.direction2 == DOWNSTREAM) {
		int base = sequence2.size()-1;
		while (base >= 0 && (sequence2[base] == 'a' || sequence2[base] == 't' || sequence2[base] == 'c' || sequence2[base] == 'g'))
			--base;
		if (base < sequence2.size()-1 && base >= 0) {
			sequence2 = sequence2.substr(0, base+1) + "|" + sequence2.substr(base+1);
			sequence2_has_non_template_bases = true;
		}
	}

	string concatenated_sequence;
	if (fusion.transcript_start == TRANSCRIPT_START_GENE1) {
		concatenated_sequence = ((fusion.predicted_strand1 == FORWARD) ? sequence1 : dna_to_reverse_complement(sequence1));
		if (!sequence1_has_non_template_bases || !sequence2_has_non_template_bases) // otherwise the concatenated sequence might have three pipes
			concatenated_sequence += "|";
		concatenated_sequence += ((fusion.direction2 == UPSTREAM) ? sequence2 : dna_to_reverse_complement(sequence2));
	} else { // fusion.transcript_start == TRANSCRIPT_START_GENE2)
		concatenated_sequence = ((fusion.predicted_strand2 == FORWARD) ? sequence2 : dna_to_reverse_complement(sequence2));
		if (!sequence2_has_non_template_bases || !sequence1_has_non_template_bases) // otherwise the concatenated sequence might have three pipes
			concatenated_sequence += "|";
		concatenated_sequence += ((fusion.direction1 == UPSTREAM) ? sequence1 : dna_to_reverse_complement(sequence1));
	}

	return concatenated_sequence;
}

bool sort_fusions_by_support(const fusion_t* x, const fusion_t* y) {
	if (x->confidence != y->confidence)
		return x->confidence > y->confidence;
	else if (x->supporting_reads() != y->supporting_reads())
		return x->supporting_reads() > y->supporting_reads();
	else if (x->evalue != y->evalue)
		return x->evalue < y->evalue;
	else
		return x->gene1->start + x->gene2->start < y->gene1->start + y->gene2->start; // this does not really sort, it only ensures that fusions between the
		                                                                              // same pair of genes are grouped together, if e-value and supporting reads are equal
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
			bool is_utr = false;
			unsigned int is_3_end = 0;
			unsigned int is_5_end = 0;
			for (exon_set_t::iterator exon = exons.begin(); exon != exons.end(); ++exon) {
				if ((**exon).gene == gene) {
					has_overlapping_exon = true;
					if ((**exon).is_utr)
						is_utr = true;
					if ((**exon).is_transcript_end)
						++is_3_end;
					if ((**exon).is_transcript_start)
						++is_5_end;
				}
			}
			if (!has_overlapping_exon) {
				site = "intron";
			} else if (is_utr) {
				if (is_3_end > is_5_end) {
					site = "3'UTR";
				} else if (is_3_end < is_5_end) {
					site = "5'UTR";
				} else {
					site = "UTR";
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

void write_fusions_to_file(fusions_t& fusions, const string& output_file, const assembly_t& assembly, gene_annotation_index_t& gene_annotation_index, exon_annotation_index_t& exon_annotation_index, vector<string> contigs_by_id, const bool print_supporting_reads, const bool print_fusion_sequence, const bool write_discarded_fusions) {
//TODO add "chr", if necessary

	// make a vector of pointers to all fusions
	// the vector will hold the fusions in sorted order
	vector<fusion_t*> sorted_fusions;
	if (write_discarded_fusions)
		sorted_fusions.reserve(fusions.size());
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion)
		if (write_discarded_fusions != (fusion->second.filter == NULL)) // either write filtered or unfiltered fusions
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
		cerr << "ERROR: Failed to open output file '" << output_file << "'." << endl;
		exit(1);
	}
	out << "#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\tsite1\tsite2\ttype\tdirection1\tdirection2\tsplit_reads1\tsplit_reads2\tdiscordant_mates\tconfidence\tclosest_genomic_breakpoint1\tclosest_genomic_breakpoint2\tfilters\tfusion_transcript\tread_identifiers" << endl;
	for (auto fusion = sorted_fusions.begin(); fusion != sorted_fusions.end(); ++fusion) {

		// describe site of breakpoint
		string site1 = get_fusion_site((**fusion).gene1, (**fusion).spliced1, (**fusion).exonic1, (**fusion).contig1, (**fusion).breakpoint1, exon_annotation_index);
		string site2 = get_fusion_site((**fusion).gene2, (**fusion).spliced2, (**fusion).exonic2, (**fusion).contig2, (**fusion).breakpoint2, exon_annotation_index);
		
		// convert closest genomic breakpoints to strings of the format <chr>:<position>(<distance to transcriptomic breakpoint>)
		string closest_genomic_breakpoint1, closest_genomic_breakpoint2;
		if ((**fusion).closest_genomic_breakpoint1 >= 0) {
			closest_genomic_breakpoint1 = contigs_by_id[(**fusion).contig1] + ":" + to_string(static_cast<long long int>((**fusion).closest_genomic_breakpoint1+1)) + "(" + to_string(static_cast<long long int>(abs((**fusion).breakpoint1 - (**fusion).closest_genomic_breakpoint1))) + ")";
		} else {
			closest_genomic_breakpoint1 = ".";
		}
		if ((**fusion).closest_genomic_breakpoint2 >= 0) {
			closest_genomic_breakpoint2 = contigs_by_id[(**fusion).contig2] + ":" + to_string(static_cast<long long int>((**fusion).closest_genomic_breakpoint2+1)) + "(" + to_string(static_cast<long long int>(abs((**fusion).breakpoint2 - (**fusion).closest_genomic_breakpoint2))) + ")";
		} else {
			closest_genomic_breakpoint2 = ".";
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
		gene_t gene1 = (**fusion).gene1; gene_t gene2 = (**fusion).gene2;
		contig_t contig1 = (**fusion).contig1; contig_t contig2 = (**fusion).contig2;
		position_t breakpoint1 = (**fusion).breakpoint1; position_t breakpoint2 = (**fusion).breakpoint2;
		direction_t direction1 = (**fusion).direction1; direction_t direction2 = (**fusion).direction2;
		strand_t strand1 = (**fusion).predicted_strand1; strand_t strand2 = (**fusion).predicted_strand2;
		unsigned int split_reads1 = (**fusion).split_reads1; unsigned int split_reads2 = (**fusion).split_reads2;
		if ((**fusion).transcript_start == TRANSCRIPT_START_GENE2) {
			swap(gene1, gene2);
			swap(direction1, direction2);
			swap(contig1, contig2);
			swap(breakpoint1, breakpoint2);
			swap(site1, site2);
			swap(split_reads1, split_reads2);
			swap(closest_genomic_breakpoint1, closest_genomic_breakpoint2);
			swap(strand1, strand2);
		}

		// write line to output file
		out << gene_to_name(gene1, contig1, breakpoint1, gene_annotation_index) << "\t" << gene_to_name(gene2, contig2, breakpoint2, gene_annotation_index) << "\t"
		    << get_fusion_strand(strand1, gene1, (**fusion).predicted_strands_ambiguous) << "\t" << get_fusion_strand(strand2, gene2, (**fusion).predicted_strands_ambiguous) << "\t"
		    << contigs_by_id[contig1] << ":" << (breakpoint1+1) << "\t" << contigs_by_id[contig2] << ":" << (breakpoint2+1) << "\t"
		    << site1 << "\t" << site2 << "\t"
		    << get_fusion_type(**fusion) << "\t" << ((direction1 == UPSTREAM) ? "upstream" : "downstream") << "\t" << ((direction2 == UPSTREAM) ? "upstream" : "downstream") << "\t"
		    << split_reads1 << "\t" << split_reads2 << "\t" << (**fusion).discordant_mates << "\t"
		    << confidence << "\t"
		    << closest_genomic_breakpoint1 << "\t" << closest_genomic_breakpoint2;

		// count the number of reads discarded by a given filter
		map<string,unsigned int> filters;
		if ((**fusion).filter != NULL)
			filters[*(**fusion).filter] = 0;
		vector<chimeric_alignments_t::iterator> all_supporting_reads;
		all_supporting_reads.insert(all_supporting_reads.end(), (**fusion).split_read1_list.begin(), (**fusion).split_read1_list.end());
		all_supporting_reads.insert(all_supporting_reads.end(), (**fusion).split_read2_list.begin(), (**fusion).split_read2_list.end());
		all_supporting_reads.insert(all_supporting_reads.end(), (**fusion).discordant_mate_list.begin(), (**fusion).discordant_mate_list.end());
		for (auto chimeric_alignment = all_supporting_reads.begin(); chimeric_alignment != all_supporting_reads.end(); ++chimeric_alignment)
			if ((**chimeric_alignment).second.filter != NULL)
				filters[*(**chimeric_alignment).second.filter]++;

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
		if (print_fusion_sequence) {
			out << get_fusion_transcript_sequence(**fusion, assembly);
		} else {
			out << ".";
		}

		// if requested, print identifiers of supporting reads
		out << "\t";
		if (print_supporting_reads && !all_supporting_reads.empty()) {
			for (auto read = all_supporting_reads.begin(); read != all_supporting_reads.end(); ++read) {
				if (read != all_supporting_reads.begin())
					out << ",";
				out << (**read).first;
			}
		} else {
			out << ".";
		}

		out << endl;
	}
	out.close();
	if (out.bad()) {
		cerr << "ERROR: Failed to write to file" << endl;
		exit(1);
	}
}

