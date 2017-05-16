#include <cmath>
#include <set>
#include <string>
#include <tuple>
#include "sam.h"
#include "common.hpp"
#include "filter_mismatches.hpp"

using namespace std;

typedef set< tuple<contig_t,position_t,char> > snp_list_t;

void count_mismatches(const alignment_t& alignment, const assembly_t& assembly, const snp_list_t& snps, unsigned int& mismatch_count, snp_list_t& mismatches, unsigned int& alignment_length) {
	// calculate template length and the number of mismatches
	mismatch_count = 0;
	mismatches.clear();
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
					mismatch_count++;
				break;
			case BAM_CDEL:
				mismatch_count++;
				// fall through
			case BAM_CREF_SKIP:
				reference_position += alignment.cigar.op_length(i);
				break;
			case BAM_CINS:
				mismatch_count++;
				read_position += alignment.cigar.op_length(i);
				break;
			case BAM_CMATCH:
				for (unsigned int operation_i = 1; operation_i <= alignment.cigar.op_length(i); ++operation_i) {
					if (alignment.sequence[read_position] != 'N') {
						if (snps.find(make_tuple(alignment.contig, reference_position, alignment.sequence[read_position])) == snps.end()) {
							if (alignment.sequence[read_position] != assembly.at(alignment.contig)[reference_position]) {
								mismatch_count++;
								mismatches.insert(make_tuple(alignment.contig, reference_position, alignment.sequence[read_position]));
							}
						}
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

void get_snps(const chimeric_alignments_t& chimeric_alignments, const assembly_t& assembly, const float mismatch_probability, const float pvalue_cutoff, snp_list_t& snps) {
	unsigned int mismatch_count, alignment_length;
	snp_list_t empty_list, potential_snps;
	for (chimeric_alignments_t::const_iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {
		if (chimeric_alignment->second.filter != NULL)
			continue; // read has already been filtered
		for (mates_t::const_iterator mate = chimeric_alignment->second.begin(); mate != chimeric_alignment->second.end(); ++mate) {
			count_mismatches(*mate, assembly, empty_list, mismatch_count, potential_snps, alignment_length);
			if (calculate_binomial(mismatch_count, alignment_length, mismatch_probability) >= pvalue_cutoff)
				snps.insert(potential_snps.begin(), potential_snps.end());
		}
	}
}

bool test_mismatch_probability(const alignment_t& alignment, const assembly_t& assembly, const snp_list_t& snps, const float mismatch_probability, const float pvalue_cutoff) {
	// estimate probability of observing the given number of mismatches by chance using a binomial model
	unsigned int mismatch_count, alignment_length;
	snp_list_t mismatches;
	count_mismatches(alignment, assembly, snps, mismatch_count, mismatches, alignment_length);
	return calculate_binomial(mismatches.size(), alignment_length, mismatch_probability) < pvalue_cutoff;
}

unsigned int filter_mismatches(chimeric_alignments_t& chimeric_alignments, const assembly_t& assembly, const float mismatch_probability, const float pvalue_cutoff) {
	snp_list_t snps;
	get_snps(chimeric_alignments, assembly, mismatch_probability, pvalue_cutoff, snps);
	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != NULL)
			continue; // read has already been filtered

		// discard chimeric alignments which have too many mismatches
		if (chimeric_alignment->second.size() == 2) { // discordant mates
			
			if (test_mismatch_probability(chimeric_alignment->second[MATE1], assembly, snps, mismatch_probability, pvalue_cutoff) ||
			    test_mismatch_probability(chimeric_alignment->second[MATE2], assembly, snps, mismatch_probability, pvalue_cutoff)) {
				chimeric_alignment->second.filter = FILTERS.at("mismatches");
				continue;
			}
		} else { // split read
			if (test_mismatch_probability(chimeric_alignment->second[MATE1], assembly, snps, mismatch_probability, pvalue_cutoff) ||
			    test_mismatch_probability(chimeric_alignment->second[SUPPLEMENTARY], assembly, snps, mismatch_probability, pvalue_cutoff)) {
				chimeric_alignment->second.filter = FILTERS.at("mismatches");
				continue;
			}
		}

		++remaining;
	}
	return remaining;
}

