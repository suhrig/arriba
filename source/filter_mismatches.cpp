#include <cmath>
#include <string>
#include "sam.h"
#include "annotation.hpp"
#include "assembly.hpp"
#include "common.hpp"
#include "filter_mismatches.hpp"

using namespace std;

void count_mismatches(const alignment_t& alignment, const string& sequence, const assembly_t& assembly, unsigned int& mismatches, unsigned int& alignment_length) {

	// calculate template length and the number of mismatches
	mismatches = 0;
	alignment_length = 0;
	position_t reference_position = alignment.start;
	position_t read_position = 0;
	for (unsigned int i = 0; i < alignment.cigar.size(); ++i) {
		switch (alignment.cigar.operation(i)) {
			case BAM_CSOFT_CLIP:
			case BAM_CHARD_CLIP:
				read_position += alignment.cigar.op_length(i);
				// clipping which might result from overlapping with the breakpoint is not counted as mismatch
				if (!(i == 0 && alignment.strand == REVERSE || i == alignment.cigar.size()-1 && alignment.strand == FORWARD))
					mismatches++;
				break;
			case BAM_CDEL:
				mismatches++;
				// fall through
			case BAM_CREF_SKIP:
				reference_position += alignment.cigar.op_length(i);
				break;
			case BAM_CINS:
				mismatches++;
				read_position += alignment.cigar.op_length(i);
				break;
			case BAM_CMATCH:
			case BAM_CEQUAL:
			case BAM_CDIFF:
				for (unsigned int operation_i = 1; operation_i <= alignment.cigar.op_length(i); ++operation_i) {
					if (sequence[read_position] != 'N') {
						if (sequence[read_position] != assembly.at(alignment.contig)[reference_position])
							mismatches++;
						alignment_length++;
					}
					reference_position++;
					read_position++;
				}
				break;
		}
	}
}

double calculate_binomial_coefficient(const unsigned int k, const unsigned int n) {
	double result = 1;

	// calculate n over m
	for (unsigned int i = n - k + 1; i <= n; ++i)
		result *= i;
	for (unsigned int i = 1; i <= k; ++i)
		result /= i;

	return result;
}

float calculate_binomial_distribution(const unsigned int k, const unsigned int n, const float p) {
	return calculate_binomial_coefficient(k, n) * pow(p, k) * pow(1-p, n-k);
}

bool test_mismatch_probability(const alignment_t& alignment, const string& sequence, const assembly_t& assembly, const float mismatch_probability, long unsigned int genome_size, const float pvalue_cutoff) {

	// Alignment artifacts with many mismatches arise from two sources:
	// 1. read incorrectly aligned to homologous sequence
	// 2. randomly generated sequence (e.g., somatic hypermutation) aligns somewhere in the genome by chance
	// The latter is more probable for short sequences, but becomes irrelevant for longer sequences.
	// The former is mostly relevant for medium-sized sequences.

	// estimate probability of observing the given number of mismatches by chance using a binomial model
	unsigned int mismatches, alignment_length;
	count_mismatches(alignment, sequence, assembly, mismatches, alignment_length);
	if (calculate_binomial_distribution(mismatches, alignment_length, mismatch_probability) < pvalue_cutoff) {
		return true;
	} else if (mismatches > 0) {
		// estimate probability that a random sequence aligns somewhere in the genome
		long double number_of_permutations_of_bases = pow(4/*#bases*/, alignment_length - mismatches);
		if (genome_size >= number_of_permutations_of_bases) { // short sequences are practically guaranteed to have a hit in the genome by random chance
			return true;
		} else {
			// discard read, if there is a >1% chance that it is a random sequence that just happens to have a match in the genome
			return (1 - pow(1 - genome_size/number_of_permutations_of_bases, calculate_binomial_coefficient(mismatches, alignment_length))) > 0.01;
		}
	} else
		return false;
}

unsigned int filter_mismatches(chimeric_alignments_t& chimeric_alignments, const assembly_t& assembly, const vector<bool>& interesting_contigs, const float mismatch_probability, const float pvalue_cutoff) {

	// calculate size of genome
	// we'll need this to calculate the probability of finding a match in the genome given a random sequence of bases
	long unsigned int genome_size = 0;
	for (assembly_t::const_iterator contig = assembly.begin(); contig != assembly.end(); ++contig)
		if (interesting_contigs[contig->first])
			genome_size += contig->second.size();

	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != NULL)
			continue; // read has already been filtered

		// discard chimeric alignments which have too many mismatches
		if (chimeric_alignment->second.size() == 2) { // discordant mates
			
			if (test_mismatch_probability(chimeric_alignment->second[MATE1], chimeric_alignment->second[MATE1].sequence, assembly, mismatch_probability, genome_size, pvalue_cutoff) ||
			    test_mismatch_probability(chimeric_alignment->second[MATE2], chimeric_alignment->second[MATE2].sequence, assembly, mismatch_probability, genome_size, pvalue_cutoff)) {
				chimeric_alignment->second.filter = FILTERS.at("mismatches");
				continue;
			}
		} else { // split read
			if (test_mismatch_probability(chimeric_alignment->second[MATE1], chimeric_alignment->second[MATE1].sequence, assembly, mismatch_probability, genome_size, pvalue_cutoff) ||
			    test_mismatch_probability(chimeric_alignment->second[SUPPLEMENTARY], (chimeric_alignment->second[SUPPLEMENTARY].strand == chimeric_alignment->second[SPLIT_READ].strand) ? chimeric_alignment->second[SPLIT_READ].sequence : dna_to_reverse_complement(chimeric_alignment->second[SPLIT_READ].sequence), assembly, mismatch_probability, genome_size, pvalue_cutoff)) {
				chimeric_alignment->second.filter = FILTERS.at("mismatches");
				continue;
			}
		}

		++remaining;
	}
	return remaining;
}

