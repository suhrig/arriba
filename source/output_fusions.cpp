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
#include "output_fusions.hpp"

using namespace std;

// try to determine which gene makes the 5' end of the transcript by
// looking at which promoter likely drives transcription
enum transcript_start_t { TRANSCRIPT_START_GENE1, TRANSCRIPT_START_GENE2, TRANSCRIPT_START_AMBIGUOUS };
transcript_start_t get_start_of_transcript(const fusion_t& fusion) {
	if (fusion.gene1->is_dummy) { // if gene1 is a dummy gene, then gene2 has priority; otherwise gene1 has priority
		// TODO this branch needs a make-over when we support strand-specific libraries
		if (fusion.gene2->strand == FORWARD && fusion.direction2 == DOWNSTREAM) { // transcript = gene2(+) -> gene1(+/-)
			return TRANSCRIPT_START_GENE2;
		} else if (fusion.gene2->strand == REVERSE && fusion.direction2 == UPSTREAM) { // transcript = gene2(-) -> gene1(+/-)
			return TRANSCRIPT_START_GENE2;
		} else if (fusion.gene2->strand == FORWARD && fusion.gene1->contig == fusion.gene2->contig && fusion.gene1->end < fusion.gene2->start) { // potential read-through fusion with upstream unannotated exon
			return TRANSCRIPT_START_GENE1;
		} else if (fusion.gene2->strand == REVERSE && fusion.gene1->contig == fusion.gene2->contig && fusion.gene1->start > fusion.gene2->end) { // potential read-through fusion with downstream unannotated exon
			return TRANSCRIPT_START_GENE1;
		} else { // ambiguous, since orientation of dummy gene is unclear
			return TRANSCRIPT_START_AMBIGUOUS;
		}
	} else if (!fusion.exonic1 && fusion.exonic2 || // if breakpoint1 is intronic and breakpoint2 isn't, then gene2 has priority
	           !fusion.spliced1 && fusion.spliced2) { // if both breakpoints are exonic, but only breakpoint2 is at a splice site, then gene2 has priority
		if (fusion.gene2->strand == FORWARD && fusion.direction2 == DOWNSTREAM) { // transcript = gene2(+) -> gene1(+/-)
			return TRANSCRIPT_START_GENE2;
		} else if (fusion.gene2->strand == REVERSE && fusion.direction2 == UPSTREAM) { // transcript = gene2(-) -> gene1(+/-)
			return TRANSCRIPT_START_GENE2;
		} else if (fusion.gene1->strand == FORWARD && fusion.direction1 == DOWNSTREAM) { // transcript = gene1(+) -> gene2(+/-)
			return TRANSCRIPT_START_GENE1;
		} else if (fusion.gene1->strand == REVERSE && fusion.direction1 == UPSTREAM) { // transcript = gene1(-) -> gene2(+/-)
			return TRANSCRIPT_START_GENE1;
		} else { // end-to-end-fused genes
			return TRANSCRIPT_START_AMBIGUOUS;
		}
	// in all other cases gene1 has priority
	} else if (fusion.gene1->strand == FORWARD && fusion.direction1 == DOWNSTREAM) { // transcript = gene1(+) -> gene2(+/-)
		return TRANSCRIPT_START_GENE1;
	} else if (fusion.gene1->strand == REVERSE && fusion.direction1 == UPSTREAM) { // transcript = gene1(-) -> gene2(+/-)
		return TRANSCRIPT_START_GENE1;
	} else if (fusion.gene2->strand == FORWARD && fusion.direction2 == DOWNSTREAM) { // transcript = gene2(+) -> gene1(+/-)
		return TRANSCRIPT_START_GENE2;
	} else if (fusion.gene2->strand == REVERSE && fusion.direction2 == UPSTREAM) { // transcript = gene2(-) -> gene1(+/-)
		return TRANSCRIPT_START_GENE2;
	} else { // end-to-end-fused genes
		return TRANSCRIPT_START_AMBIGUOUS;
	}
}

typedef map< position_t, map<string/*base*/,unsigned int/*frequency*/> > pileup_t;

void pileup_chimeric_alignments(vector<mates_t*>& chimeric_alignments, const unsigned int mate, const bool reverse_complement, const direction_t direction, const position_t breakpoint, pileup_t& pileup) {
	for (auto i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {

		if ((**i).filter == FILTERS.at("duplicates"))
			continue; // skip duplicates

		alignment_t& read = (**i)[mate]; // introduce alias for cleaner code

		if ((**i).size() == 2) // discordant mate
			if (!(direction == DOWNSTREAM && read.strand == FORWARD && read.end   <= breakpoint+2 && read.end   >= breakpoint-200 ||
			      direction == UPSTREAM   && read.strand == REVERSE && read.start >= breakpoint-2 && read.start <= breakpoint+200)) // only consider discordant mates close to the breakpoints (we don't care about the ones in other exons)
				continue;

		string read_sequence = (mate == SUPPLEMENTARY) ? (**i)[SPLIT_READ].sequence : read.sequence;
		if (reverse_complement)
			read_sequence = dna_to_reverse_complement(read_sequence);

		position_t read_offset = 0;
		position_t reference_offset = read.start;
		int subtract_from_next_element = 0;
		for (int cigar_element = 0; cigar_element < read.cigar.size(); cigar_element++) {
			switch (read.cigar.operation(cigar_element)) {
				case BAM_CINS:
					pileup[reference_offset][read_sequence.substr(read_offset, read.cigar.op_length(cigar_element)+1)]++;
					read_offset += read.cigar.op_length(cigar_element) + 1; // +1, because we take one base from the next element
					++reference_offset; // +1, because we take one base from the next element
					subtract_from_next_element = 1; // because we took one base from the next element
					break;
				case BAM_CREF_SKIP:
					reference_offset += read.cigar.op_length(cigar_element) - subtract_from_next_element;
					subtract_from_next_element = 0;
					break;
				case BAM_CDEL:
					for (position_t base = 0; base < read.cigar.op_length(cigar_element) - subtract_from_next_element; ++base, ++reference_offset)
						pileup[reference_offset]["-"]++; // indicate deletion by dash
					subtract_from_next_element = 0;
					break;
				case BAM_CSOFT_CLIP:
					if (mate == SPLIT_READ &&
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
}

void get_sequence_from_pileup(const pileup_t& pileup, const position_t breakpoint, const direction_t direction, const gene_t gene,  string& sequence, string& clipped_sequence) {

	// for each position, find the most frequent allele in the pileup
	position_t previous_position;
	for (pileup_t::const_iterator position = pileup.begin(); position != pileup.end(); ++position) {

		if (position != pileup.begin() && previous_position < position->first - 1)
			sequence += "..."; // indicate introns/gaps with an ellipsis
		previous_position = position->first;

		// find out base in reference to mark SNPs/SNVs
		string reference_base = "N";
		if (position->first >= gene->start && // check if we have the sequence for the given position
		    gene->sequence.size() > position->first - gene->start)
			reference_base = gene->sequence[position->first - gene->start];

		// find most frequent allele at current position and compute coverage
		auto base = position->second.begin();
		auto most_frequent_base = base;
		unsigned int coverage = most_frequent_base->second;
		for (++base; base != position->second.end(); ++base) {
			if (base->second > most_frequent_base->second || base->first == reference_base && base->second+1 >= most_frequent_base->second)
				most_frequent_base = base;
			coverage += base->second;
		}

		// we trust the base, if it has a frequency of >= 75%
		string most_frequent_base2 = (most_frequent_base->second >= 0.75 * coverage || most_frequent_base->first == reference_base) ? most_frequent_base->first : "n";

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

string get_fusion_transcript_sequence(fusion_t& fusion, const transcript_start_t transcript_start) {

	if (transcript_start == TRANSCRIPT_START_AMBIGUOUS)
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
			int clipped_split_read = ((**read)[SPLIT_READ].strand == FORWARD) ? (**read)[SPLIT_READ].cigar.op_length(0) : (**read)[SPLIT_READ].cigar.op_length((**read)[SPLIT_READ].cigar.size()-1);
			int clipped_supplementary = ((**read)[SUPPLEMENTARY].strand == FORWARD) ? (**read)[SUPPLEMENTARY].cigar.op_length((**read)[SUPPLEMENTARY].cigar.size()-1) : (**read)[SUPPLEMENTARY].cigar.op_length(0);
			int unmapped_bases = clipped_split_read + clipped_supplementary - (**read)[SPLIT_READ].sequence.size();
			if (++non_template_bases_count[unmapped_bases] > non_template_bases_count[non_template_bases])
				non_template_bases = unmapped_bases;
		}

	}

	// determine most frequent bases in pileup
	string sequence1, sequence2, clipped_sequence1, clipped_sequence2;
	get_sequence_from_pileup(pileup1, fusion.breakpoint1, fusion.direction1, fusion.gene1, sequence1, clipped_sequence1);
	get_sequence_from_pileup(pileup2, fusion.breakpoint2, fusion.direction2, fusion.gene2, sequence2, clipped_sequence2);

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

	string concatenated_sequence = ".";
	if (transcript_start == TRANSCRIPT_START_GENE1) {
		concatenated_sequence = ((fusion.gene1->strand == FORWARD) ? sequence1 : dna_to_reverse_complement(sequence1));
		if (!sequence1_has_non_template_bases || !sequence2_has_non_template_bases) // otherwise the concatenated sequence might have three pipes
			concatenated_sequence += "|";
		concatenated_sequence += ((fusion.direction2 == UPSTREAM) ? sequence2 : dna_to_reverse_complement(sequence2));
	} else { // transcript_start == TRANSCRIPT_START_GENE2)
		concatenated_sequence = ((fusion.gene2->strand == FORWARD) ? sequence2 : dna_to_reverse_complement(sequence2));
		if (!sequence2_has_non_template_bases || !sequence1_has_non_template_bases) // otherwise the concatenated sequence might have three pipes
			concatenated_sequence += "|";
		concatenated_sequence += ((fusion.direction1 == UPSTREAM) ? sequence1 : dna_to_reverse_complement(sequence1));
	}

	return concatenated_sequence;
}

bool sort_fusions_by_support(const fusion_t* x, const fusion_t* y) {
	if (x->supporting_reads() != y->supporting_reads())
		return x->supporting_reads() > y->supporting_reads();
	else if (x->evalue != y->evalue)
		return x->evalue < y->evalue;
	else
		return x->gene1->start + x->gene2->start < y->gene1->start + y->gene2->start; // this does not really sort, it only ensures that fusions between the
		                                                                              // same pair of genes are grouped together, if e-value and supporting reads are equal
}

typedef tuple<gene_t /*gene1*/, gene_t /*gene2*/> gene_pair_t;
class sort_fusions_by_rank_of_best_t {
public:
	unordered_map< gene_pair_t, fusion_t* > best;
	bool operator()(const fusion_t* x, const fusion_t* y) {
		fusion_t* best_x = best.at(make_tuple(x->gene1, x->gene2));
		fusion_t* best_y = best.at(make_tuple(y->gene1, y->gene2));
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
			for (gene_set_t::iterator gene = index_hit1->second.begin(); gene != index_hit1->second.end(); gene = index_hit1->second.upper_bound(*gene)) {
				if (!(**gene).is_dummy) {
					if (!result.empty())
						result += ",";
					result += (**gene).name + "(" + to_string(breakpoint - (**gene).end) + ")";
				}
			}
		}

		// go downstream until we find a non-dummy gene
		while (index_hit2 != gene_annotation_index[contig].end() && (index_hit2->second.empty() || (**(index_hit2->second.begin())).is_dummy))
			++index_hit2;

		// append downstream flanking genes with distances to gene name
		if (index_hit2 != gene_annotation_index[contig].end()) {
			for (gene_set_t::iterator gene = index_hit2->second.begin(); gene != index_hit2->second.end(); gene = index_hit2->second.upper_bound(*gene)) {
				if (!(**gene).is_dummy) {
					if (!result.empty())
						result += ",";
					result += (**gene).name + "(" + to_string((**gene).start - breakpoint) + ")";
				}
			}
		}

		return result;
	}
}

string get_fusion_type(const fusion_t& fusion) {
	if (fusion.contig1 != fusion.contig2) {
		if (fusion.gene1->is_dummy || fusion.gene2->is_dummy) {
			return "translocation/truncation"; // fusion ends up in intergenic region => gene is truncated
		} else if (fusion.direction1 == fusion.direction2 && fusion.gene1->strand != fusion.gene2->strand ||
		         fusion.direction1 != fusion.direction2 && fusion.gene1->strand == fusion.gene2->strand) {
			return "translocation"; // orderly fusion yielding a hybrid protein
		} else {
			if (fusion.direction1 == UPSTREAM && fusion.direction2 == UPSTREAM && fusion.gene1->strand == FORWARD && fusion.gene2->strand == FORWARD ||
			    fusion.direction1 == DOWNSTREAM && fusion.direction2 == DOWNSTREAM && fusion.gene1->strand == REVERSE && fusion.gene2->strand == REVERSE) {
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
					return "deletion/read-through";
				} else {
					return "deletion/5'-5'";
				}
			} else {
					return "deletion/3'-3'"; // tail-to-tail fusion
			}
		} else if (fusion.direction1 == fusion.direction2) {
			if (fusion.gene1->strand != fusion.gene2->strand) {
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
			if (fusion.gene1->strand == fusion.gene2->strand) {
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

string get_fusion_strand(const gene_t gene1, const gene_t gene2, const direction_t direction1, const direction_t direction2, const bool is_start, const transcript_start_t& transcript_start) {
	string result;

	// determine strand of gene as per annotation
	if (gene1->is_dummy)
		result = ".";
	else
		result = (gene1->strand == FORWARD) ? "+" : "-";

	// separator between gene strand and fusion strand
	result += "/";

	// determine most likely strand of fusion from the directions of the fusion
	if (transcript_start == TRANSCRIPT_START_AMBIGUOUS) {
		result += ".";
	} else {
		if (is_start) {
			if (gene1->is_dummy) {
				if (gene2->is_dummy) {
					result += ".";
				} else {
					if (direction1 != direction2) {
						result += (gene2->strand == FORWARD) ? "+" : "-";
					} else {
						result += (gene2->strand == FORWARD) ? "-" : "+";
					}
				}
			} else {
				result += (gene1->strand == FORWARD) ? "+" : "-";
			}
		} else { // !is_start
			if (gene2->is_dummy) {
				if (gene1->is_dummy) {
					result += ".";
				} else {
					result += (gene1->strand == FORWARD) ? "+" : "-";
				}
			} else {
				if (direction1 != direction2) {
					result += (gene2->strand == FORWARD) ? "+" : "-";
				} else {
					result += (gene2->strand == FORWARD) ? "-" : "+";
				}
			}
		}
	}

	return result;
}

enum confidence_t { LOW_CONFIDENCE, MEDIUM_CONFIDENCE, HIGH_CONFIDENCE };
confidence_t get_confidence(const fusion_t& fusion, const map< gene_t, vector<fusion_t*> >& fusions_by_gene) {
	confidence_t result;

	if (fusion.filter != NULL) {
		result = LOW_CONFIDENCE;
	} else {
		if (fusion.evalue > 0.3) {
			result = LOW_CONFIDENCE;
			if (fusion.spliced1 && fusion.spliced2 && !fusion.is_read_through()) {
				// look for multiple spliced translocations
				unsigned int number_of_spliced_breakpoints = 0;
				auto fusions_of_gene = fusions_by_gene.find(fusion.gene1);
				for (auto fusion_of_gene = fusions_of_gene->second.begin(); fusion_of_gene != fusions_of_gene->second.end(); ++fusion_of_gene) {
					if ((**fusion_of_gene).gene2 == fusion.gene2 && (**fusion_of_gene).gene2 == fusion.gene2 && (**fusion_of_gene).spliced1 && (**fusion_of_gene).spliced2)
						++number_of_spliced_breakpoints;
				}
				fusions_of_gene = fusions_by_gene.find(fusion.gene2);
				for (auto fusion_of_gene = fusions_of_gene->second.begin(); fusion_of_gene != fusions_of_gene->second.end(); ++fusion_of_gene) {
					if ((**fusion_of_gene).gene1 == fusion.gene1 && (**fusion_of_gene).gene2 == fusion.gene2 && (**fusion_of_gene).spliced1 && (**fusion_of_gene).spliced2)
						++number_of_spliced_breakpoints;
				}
				if (number_of_spliced_breakpoints >= 2)
					result = MEDIUM_CONFIDENCE;
			}
		} else if (fusion.is_read_through()) {
			result = LOW_CONFIDENCE;
			if ((fusion.split_reads1 > 0 && fusion.split_reads2 > 0 || fusion.split_reads1 > 0 && fusion.discordant_mates > 0 || fusion.split_reads2 > 0 && fusion.discordant_mates > 0) && fusion.supporting_reads() > 9) {
				result = MEDIUM_CONFIDENCE;
			} else {
				// look for multiple deletions involving the same gene
				unsigned int number_of_deletions = 0;
				auto fusions_of_gene = fusions_by_gene.find(fusion.gene1);
				for (auto fusion_of_gene = fusions_of_gene->second.begin(); fusion_of_gene != fusions_of_gene->second.end(); ++fusion_of_gene) {
					if ((**fusion_of_gene).filter == NULL &&
					    (**fusion_of_gene).split_reads1 + (**fusion_of_gene).split_reads2 > 0 &&
					    (**fusion_of_gene).direction1 == DOWNSTREAM && (**fusion_of_gene).direction2 == UPSTREAM &&
					    ((**fusion_of_gene).breakpoint1 != fusion.breakpoint1 || (**fusion_of_gene).breakpoint2 != fusion.breakpoint2) &&
					    (**fusion_of_gene).breakpoint2 > fusion.breakpoint1 && (**fusion_of_gene).breakpoint1 < fusion.breakpoint2) {
						++number_of_deletions;
					}
				}
				fusions_of_gene = fusions_by_gene.find(fusion.gene2);
				for (auto fusion_of_gene = fusions_of_gene->second.begin(); fusion_of_gene != fusions_of_gene->second.end(); ++fusion_of_gene) {
					if ((**fusion_of_gene).filter == NULL &&
					    (**fusion_of_gene).split_reads1 + (**fusion_of_gene).split_reads2 > 0 &&
					    (**fusion_of_gene).direction1 == DOWNSTREAM && (**fusion_of_gene).direction2 == UPSTREAM &&
					    ((**fusion_of_gene).breakpoint1 != fusion.breakpoint1 || (**fusion_of_gene).breakpoint2 != fusion.breakpoint2) &&
					    (**fusion_of_gene).breakpoint2 > fusion.breakpoint1 && (**fusion_of_gene).breakpoint1 < fusion.breakpoint2) {
						++number_of_deletions;
						}
				}
				if (number_of_deletions >= 1)
					result = MEDIUM_CONFIDENCE;
			}

		} else if (fusion.breakpoint_overlaps_both_genes() || fusion.gene1 == fusion.gene2) { // intragenic event
			if (fusion.split_reads1 + fusion.split_reads2 == 0) {
				result = LOW_CONFIDENCE;
			} else if (!fusion.exonic1 && !fusion.exonic2) {
				if (fusion.split_reads1 > 0 && fusion.split_reads2 > 0) {
					result = HIGH_CONFIDENCE;
				} else {
					result = MEDIUM_CONFIDENCE;
				}
			} else if (!fusion.exonic1 || !fusion.exonic2) {
				if (fusion.split_reads1 > 3 && fusion.split_reads2 > 3) {
					result = HIGH_CONFIDENCE;
				} else {
					result = MEDIUM_CONFIDENCE;
				}
			} else {
				result = LOW_CONFIDENCE;
			}
		} else if (fusion.split_reads1 + fusion.split_reads2 == 0 || fusion.split_reads1 + fusion.discordant_mates == 0 || fusion.split_reads2 + fusion.discordant_mates == 0) {
			result = MEDIUM_CONFIDENCE;
		} else {
			result = HIGH_CONFIDENCE;
		}

		if (fusion.evalue > 0.2) {
			if (result == HIGH_CONFIDENCE)
				result = MEDIUM_CONFIDENCE;
			else if (result == MEDIUM_CONFIDENCE)
				result = LOW_CONFIDENCE;
		}

		if (fusion.closest_genomic_breakpoint1 >= 0) { // has genomic support
			if (result == LOW_CONFIDENCE)
				result = MEDIUM_CONFIDENCE;
			else if (result == MEDIUM_CONFIDENCE)
				result = HIGH_CONFIDENCE;
		}
	}

	return result;
}

void write_fusions_to_file(fusions_t& fusions, const string& output_file, gene_annotation_index_t& gene_annotation_index, exon_annotation_index_t& exon_annotation_index, vector<string> contigs_by_id, const bool print_supporting_reads, const bool print_fusion_sequence, const bool write_discarded_fusions) {
//TODO add "chr", if necessary

	// make a vector of pointers to all fusions
	// the vector will hold the fusions in sorted order
	vector<fusion_t*> sorted_fusions;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
		if (write_discarded_fusions || i->second.filter == NULL)
			sorted_fusions.push_back(&(i->second));
	}

	// don't sort the discarded fusions
	if (!write_discarded_fusions) {

		// sometimes we report multiple breakpoints for the same pair of genes
		// in such cases, we want to make sure that all breakpoints of the same pair of genes
		// are grouped together in the output file, even if some of them have very low support
		// all such breakpoints will get the rank of the best ranking breakpoints
		// => find out what the best ranking breakpoints are for each pair of genes
		//    we store these breakpoints in a sorting object (sort_fusions_by_rank_of_best),
		//    which we use as a parameter to the sort function later
		sort_fusions_by_rank_of_best_t sort_fusions_by_rank_of_best;
		for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
			if (i->second.filter != NULL)
				continue; // skip discarded fusions

			gene_pair_t gene_pair = make_tuple(i->second.gene1, i->second.gene2);
			if (sort_fusions_by_rank_of_best.best.find(gene_pair) == sort_fusions_by_rank_of_best.best.end()) {
				sort_fusions_by_rank_of_best.best[gene_pair] = &(i->second);
			} else {
				fusion_t*& current_best = sort_fusions_by_rank_of_best.best[gene_pair];
				if (sort_fusions_by_support(&(i->second), current_best))
					current_best = &(i->second);
			}
		}

		// sort all gene pairs by the rank of the best scoring breakpoints of a given gene pair
		sort(sorted_fusions.begin(), sorted_fusions.end(), ref(sort_fusions_by_rank_of_best));
	}

	// get all fusions by gene, we need this for better confidence scoring
	map< gene_t, vector<fusion_t*> > fusions_by_gene;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		fusions_by_gene[fusion->second.gene1].push_back(&fusion->second);
		fusions_by_gene[fusion->second.gene2].push_back(&fusion->second);
	}

	// write sorted list to file
	ofstream out(output_file);
	if (!out.is_open()) {
		cerr << "ERROR: Failed to open output file '" << output_file << "'." << endl;
		exit(1);
	}
	out << "#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\tsite1\tsite2\ttype\tdirection1\tdirection2\tsplit_reads1\tsplit_reads2\tdiscordant_mates\tconfidence\tclosest_genomic_breakpoint1\tclosest_genomic_breakpoint2\tfilters\tfusion_transcript\tread_identifiers" << endl;
	for (vector<fusion_t*>::iterator i = sorted_fusions.begin(); i != sorted_fusions.end(); ++i) {

		if (((**i).filter == NULL) == write_discarded_fusions) // write either filtered or unfiltered fusions, but not both
			continue;

		// get the gene which likely makes the 5' end of the transcript first
		transcript_start_t transcript_start = get_start_of_transcript(**i);

		// prepare columns
		gene_t gene1 = (**i).gene1; gene_t gene2 = (**i).gene2;
		contig_t contig1 = (**i).contig1; contig_t contig2 = (**i).contig2;
		position_t breakpoint1 = (**i).breakpoint1; position_t breakpoint2 = (**i).breakpoint2;
		direction_t direction1 = (**i).direction1; direction_t direction2 = (**i).direction2;
		unsigned int split_reads1 = (**i).split_reads1; unsigned int split_reads2 = (**i).split_reads2;
		string site1, site2;
		if ((**i).gene1->is_dummy) {
			site1 = "intergenic";
		} else if ((**i).spliced1) {
			site1 = "splice-site";
		} else if ((**i).exonic1) {
			site1 = "exonic";
		} else {
			site1 = "intronic";
		}
		if ((**i).gene2->is_dummy) {
			site2 = "intergenic";
		} else if ((**i).spliced2) {
			site2 = "splice-site";
		} else if ((**i).exonic2) {
			site2 = "exonic";
		} else {
			site2 = "intronic";
		}
		string closest_genomic_breakpoint1, closest_genomic_breakpoint2;
		if ((**i).closest_genomic_breakpoint1 >= 0) {
			closest_genomic_breakpoint1 = contigs_by_id[(**i).contig1] + ":" + to_string((**i).closest_genomic_breakpoint1+1) + "(" + to_string(abs((**i).breakpoint1 - (**i).closest_genomic_breakpoint1)) + ")";
		} else {
			closest_genomic_breakpoint1 = ".";
		}
		if ((**i).closest_genomic_breakpoint2 >= 0) {
			closest_genomic_breakpoint2 = contigs_by_id[(**i).contig2] + ":" + to_string((**i).closest_genomic_breakpoint2+1) + "(" + to_string(abs((**i).breakpoint2 - (**i).closest_genomic_breakpoint2)) + ")";
		} else {
			closest_genomic_breakpoint2 = ".";
		}
		string confidence;
		switch (get_confidence(**i, fusions_by_gene)) {
			case LOW_CONFIDENCE:
				confidence = "low";
				break;
			case MEDIUM_CONFIDENCE:
				confidence = "medium";
				break;
			case HIGH_CONFIDENCE:
				confidence = "high";
				break;
		}

		// the 5' gene should always come first => swap columns, if necessary
		if (transcript_start == TRANSCRIPT_START_GENE2) {
			swap(gene1, gene2);
			swap(direction1, direction2);
			swap(contig1, contig2);
			swap(breakpoint1, breakpoint2);
			swap(site1, site2);
			swap(split_reads1, split_reads2);
			swap(closest_genomic_breakpoint1, closest_genomic_breakpoint2);
		}

		// write line to output file
		out << gene_to_name(gene1, contig1, breakpoint1, gene_annotation_index) << "\t" << gene_to_name(gene2, contig2, breakpoint2, gene_annotation_index) << "\t"
		    << get_fusion_strand(gene1, gene2, direction1, direction2, true, transcript_start) << "\t" << get_fusion_strand(gene2, gene1, direction2, direction1, false, transcript_start) << "\t"
		    << contigs_by_id[contig1] << ":" << (breakpoint1+1) << "\t" << contigs_by_id[contig2] << ":" << (breakpoint2+1) << "\t"
		    << site1 << "\t" << site2 << "\t"
		    << get_fusion_type(**i) << "\t" << ((direction1 == UPSTREAM) ? "upstream" : "downstream") << "\t" << ((direction2 == UPSTREAM) ? "upstream" : "downstream") << "\t"
		    << split_reads1 << "\t" << split_reads2 << "\t" << (**i).discordant_mates << "\t"
		    << confidence << "\t"
		    << closest_genomic_breakpoint1 << "\t" << closest_genomic_breakpoint2;

		// count the number of reads discarded by a given filter
		map<string,unsigned int> filters;
		if ((**i).filter != NULL)
			filters[*(**i).filter] = 0;
		vector<mates_t*> all_supporting_reads = (**i).split_read1_list;
		all_supporting_reads.insert(all_supporting_reads.end(), (**i).split_read1_list.begin(), (**i).split_read1_list.end());
		all_supporting_reads.insert(all_supporting_reads.end(), (**i).split_read2_list.begin(), (**i).split_read2_list.end());
		all_supporting_reads.insert(all_supporting_reads.end(), (**i).discordant_mate_list.begin(), (**i).discordant_mate_list.end());
		for (auto chimeric_alignment = all_supporting_reads.begin(); chimeric_alignment != all_supporting_reads.end(); ++chimeric_alignment)
			if ((**chimeric_alignment).filter != NULL)
				filters[*(**chimeric_alignment).filter]++;

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
			out << get_fusion_transcript_sequence(**i, transcript_start);
		} else {
			out << ".";
		}

		// if requested, print identifiers of supporting reads
		out << "\t";
		if (print_supporting_reads) {
			for (auto j = all_supporting_reads.begin(); j != all_supporting_reads.end(); ++j) {
				if (j != all_supporting_reads.begin())
					out << ",";
				out << (**j).name;
			}
		} else {
			out << ".";
		}

		out << endl;
	}
	out.close();
}

