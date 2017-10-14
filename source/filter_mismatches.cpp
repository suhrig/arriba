#include <cmath>
#include <string>
#include "sam.h"
#include "annotation.hpp"
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
				read_position += alignment.cigar.op_length(i);
				// fall through
			case BAM_CHARD_CLIP:
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

// calculate the binomial distribution
float calculate_binomial(const unsigned int k, const unsigned int n, const float p) {
	float result = 1;

	// calculate n over m
	for (int i = n - k + 1; i <= n; ++i)
		result *= i;
	for (int i = 1; i <= k; ++i)
		result /= i;

	// calculate p^k * (1-p)^(n-k)
	result *= pow(p, k) * pow(1-p, n-k);

	return result;
}

bool test_mismatch_probability(const alignment_t& alignment, const string& sequence, const assembly_t& assembly, const float mismatch_probability, const float pvalue_cutoff) {
	// estimate probability of observing the given number of mismatches by chance using a binomial model
	unsigned int mismatches, alignment_length;
	count_mismatches(alignment, sequence, assembly, mismatches, alignment_length);
	return calculate_binomial(mismatches, alignment_length, mismatch_probability) < pvalue_cutoff;
}

unsigned int filter_mismatches(chimeric_alignments_t& chimeric_alignments, const assembly_t& assembly, const float mismatch_probability, const float pvalue_cutoff) {
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != NULL)
			continue; // read has already been filtered

		// discard chimeric alignments which have too many mismatches
		if (chimeric_alignment->second.size() == 2) { // discordant mates
			
			if (test_mismatch_probability(chimeric_alignment->second[MATE1], chimeric_alignment->second[MATE1].sequence, assembly, mismatch_probability, pvalue_cutoff) ||
			    test_mismatch_probability(chimeric_alignment->second[MATE2], chimeric_alignment->second[MATE2].sequence, assembly, mismatch_probability, pvalue_cutoff)) {
				chimeric_alignment->second.filter = FILTERS.at("mismatches");
				continue;
			}
		} else { // split read
			if (test_mismatch_probability(chimeric_alignment->second[MATE1], chimeric_alignment->second[MATE1].sequence, assembly, mismatch_probability, pvalue_cutoff) ||
			    test_mismatch_probability(chimeric_alignment->second[SUPPLEMENTARY], (chimeric_alignment->second[SUPPLEMENTARY].strand == chimeric_alignment->second[SPLIT_READ].strand) ? chimeric_alignment->second[SPLIT_READ].sequence : dna_to_reverse_complement(chimeric_alignment->second[SPLIT_READ].sequence), assembly, mismatch_probability, pvalue_cutoff)) {
				chimeric_alignment->second.filter = FILTERS.at("mismatches");
				continue;
			}
		}

		++remaining;
	}
	return remaining;
}

