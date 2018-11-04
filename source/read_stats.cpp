#include <climits>
#include <iostream>
#include <list>
#include <vector>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "read_stats.hpp"

bool estimate_mate_gap_distribution(const chimeric_alignments_t& chimeric_alignments, float& mate_gap_mean, float& mate_gap_stddev, const gene_annotation_index_t& gene_annotation_index, const exon_annotation_index_t& exon_annotation_index) {

	// use MATE1 and MATE2 from split reads to calculate insert size distribution
	list<int> mate_gaps;
	unsigned int count = 0;
	for (chimeric_alignments_t::const_iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {
		if (chimeric_alignment->second.filter != NULL || chimeric_alignment->second.single_end)
			continue;
		if (chimeric_alignment->second.size() == 3) { // split read

			alignment_t const* forward_mate = &(chimeric_alignment->second[MATE1]);
			alignment_t const* reverse_mate = &(chimeric_alignment->second[SPLIT_READ]);
			if (forward_mate->strand == REVERSE)
				swap(forward_mate, reverse_mate);

			mate_gaps.push_back(get_spliced_distance(forward_mate->contig, forward_mate->end, reverse_mate->start, DOWNSTREAM, UPSTREAM, *forward_mate->genes.begin(), exon_annotation_index));

			count++;
			if (count > 100000)
				break; // the sample size should be big enough
		}
	}

	if (count < 10000) {
		cerr << "WARNING: not enough chimeric reads to estimate mate gap distribution, using default values" << endl;
		return false;
	}

	bool no_more_outliers = false;
	while (true) {
		// calculate mean
		mate_gap_mean = 0;
		for (list<int>::iterator i = mate_gaps.begin(); i != mate_gaps.end(); i++)
			mate_gap_mean += *i;
		mate_gap_mean /= count;

		// calculate standard deviation
		mate_gap_stddev = 0;
		for (list<int>::iterator i = mate_gaps.begin(); i != mate_gaps.end(); i++)
			mate_gap_stddev += (*i - mate_gap_mean) * (*i - mate_gap_mean);
		mate_gap_stddev = sqrt(1.0/(count-1) * mate_gap_stddev);

		// due to alterantive splicing the mate gap distribution is not distributed normally
		// there are usually many outliers which inflate the standard deviation
		// => remove outliers until the distribution resembles a normal distribution
		//    (i.e., 68.24% are inside the range: mean +/- 1*stddev)
		unsigned int within_range = 0;
		for (list<int>::iterator i = mate_gaps.begin(); i != mate_gaps.end(); ++i)
			if (*i > mate_gap_mean - mate_gap_stddev || *i < mate_gap_mean + mate_gap_stddev)
				within_range++;
		if (1.0*within_range/count < 0.683 || no_more_outliers)
			break; // all outliers have been removed

		// remove outliers, if the mate gap distribution is not yet normally distributed
		no_more_outliers = true;
		for (list<int>::iterator i = mate_gaps.begin(); i != mate_gaps.end();) {
			if (*i < mate_gap_mean - 3*mate_gap_stddev || *i > mate_gap_mean + 3*mate_gap_stddev) {
				i = mate_gaps.erase(i); // remove outlier
				count--;
				no_more_outliers = false;
			} else {
				++i;
			}
		}
	}

	return true;
}

strandedness_t detect_strandedness(const chimeric_alignments_t& chimeric_alignments, const gene_annotation_index_t& gene_annotation_index) {
	unsigned int count = 0;
	unsigned int matching_strand = 0;
	for (chimeric_alignments_t::const_iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {
		if (chimeric_alignment->second.filter == NULL) {
			gene_set_t genes;
			alignment_t const* first_in_pair = (chimeric_alignment->second[MATE1].first_in_pair) ? &(chimeric_alignment->second[MATE1]) : &(chimeric_alignment->second[MATE2]);
			get_annotation_by_coordinate(first_in_pair->contig, first_in_pair->start, first_in_pair->end, genes, gene_annotation_index);
			if (genes.size() == 1) {
				if (first_in_pair->strand == (**genes.begin()).strand)
					matching_strand++;
				count++;
				if (count >= 100000)
					break; // the sample size should be sufficient
			}
		}
	}

	if (count < 10000) {
		return STRANDEDNESS_NO; // not enough statistical power => assume no
	} else if (matching_strand < 0.3 * count) {
		return STRANDEDNESS_REVERSE;
	} else if (matching_strand > 0.7 * count) {
		return STRANDEDNESS_YES;
	} else
		return STRANDEDNESS_NO; // not enough signal => assume no
}

// initialize data structure to compute coverage for windows of size <COVERAGE_RESOLUTION>
coverage_t::coverage_t(const contigs_t& contigs, const assembly_t& assembly) {
	fragment_starts.resize(contigs.size());
	fragment_ends.resize(contigs.size());
	coverage.resize(contigs.size());
	for (assembly_t::const_iterator contig = assembly.begin(); contig != assembly.end(); ++contig) {
		if (!contig->second.empty()) {
			unsigned int windows = contig->second.size() / COVERAGE_RESOLUTION + 2; //+2 to avoid array-out-of-bounds errors
			fragment_starts[contig->first].resize(windows);
			fragment_ends[contig->first].resize(windows);
			coverage[contig->first].resize(windows);
		}
	}
}

// add alignment to coverage
void coverage_t::add_fragment(bam1_t* mate1, bam1_t* mate2, const bool is_read_through_alignment) {

	// fake paired-end data, if single-end data given, to avoid NULL pointer exceptions
	if (mate2 == NULL)
		mate2 = mate1;

	if (mate1->core.tid >= fragment_starts.size() || fragment_starts[mate1->core.tid].empty() ||
	    mate2->core.tid >= fragment_starts.size() || fragment_starts[mate2->core.tid].empty())
		return; // ignore reads on uninteresting contigs

	bool is_chimeric = is_read_through_alignment;
	if (mate1->core.flag & BAM_FPAIRED) { // paired-end data
		if (!(mate1->core.flag & BAM_FPROPER_PAIR) || // discordant mates
	            !(mate1->core.flag & BAM_FREVERSE) && bam_cigar_type(bam_cigar_op(bam_get_cigar(mate1)[0])) == BAM_CSOFT_CLIP || // clipped at start => likely split-read
	            !(mate2->core.flag & BAM_FREVERSE) && bam_cigar_type(bam_cigar_op(bam_get_cigar(mate2)[0])) == BAM_CSOFT_CLIP ||
	            (mate1->core.flag & BAM_FREVERSE) && bam_cigar_type(bam_cigar_op(bam_get_cigar(mate1)[mate1->core.n_cigar-1])) == BAM_CSOFT_CLIP || // clipped at end => likely split-read
	            (mate2->core.flag & BAM_FREVERSE) && bam_cigar_type(bam_cigar_op(bam_get_cigar(mate2)[mate2->core.n_cigar-1])) == BAM_CSOFT_CLIP)
			is_chimeric = true;
	} else { // single-end data
		if (bam_cigar_type(bam_cigar_op(bam_get_cigar(mate1)[0])) == BAM_CSOFT_CLIP ||
		    bam_cigar_type(bam_cigar_op(bam_get_cigar(mate1)[mate1->core.n_cigar-1])) == BAM_CSOFT_CLIP)
			is_chimeric = true; // alignment is clipped => likely split-read
	}

	// store start of fragment
	if (!is_chimeric) { // the 'no_coverage' filter should only consider non-chimeric reads
		if (!(mate1->core.flag & BAM_FREVERSE) || !(mate1->core.flag & BAM_FPAIRED))
			fragment_starts[mate1->core.tid][mate1->core.pos/COVERAGE_RESOLUTION] = true;
		else
			fragment_starts[mate2->core.tid][mate2->core.pos/COVERAGE_RESOLUTION] = true;
	}

	// compute coverage from CIGAR string
	position_t position1 = mate1->core.pos;
	position_t position2 = mate2->core.pos;
	position_t position = min(position1, position2);
	int window = position/COVERAGE_RESOLUTION;
	unsigned int i1 = 0;
	unsigned int i2 = 0;
	while (true) {

		// go through CIGAR operations of both mates in parallel,
		// always choosing the one that consumes the smallest amount of the reference
		uint32_t cigar_op1, cigar_op2;
		unsigned int op_length1;
		unsigned int op_length2;
		// find out how much of the reference the next CIGAR element of each mate would consume
		if (i1 < mate1->core.n_cigar) {
			cigar_op1 = bam_get_cigar(mate1)[i1];
			op_length1 = (bam_cigar_type(bam_cigar_op(cigar_op1)) & 2/*consume reference*/) ? bam_cigar_oplen(cigar_op1) : 0;
		} else { // CIGAR elements of mate1 completely processed
			op_length1 = 0;
			window = max(window, position2/COVERAGE_RESOLUTION);
		}
		if (i2 < mate2->core.n_cigar) {
			cigar_op2 = bam_get_cigar(mate2)[i2];
			op_length2 = (bam_cigar_type(bam_cigar_op(cigar_op2)) & 2/*consume reference*/) ? bam_cigar_oplen(cigar_op2) : 0;
		} else { // CIGAR elements of mate2 completely processed
			op_length2 = 0;
			window = max(window, position1/COVERAGE_RESOLUTION);
		}
		// pick the mate whose next CIGAR element consumes the least amount of reference
		contig_t contig;
		uint32_t cigar_op;
		if (i1 < mate1->core.n_cigar && (position1 + op_length1 < position2 + op_length2 || i2 >= mate2->core.n_cigar)) {
			i1++;
			if (op_length1 == 0)
				continue;
			cigar_op = cigar_op1;
			contig = mate1->core.tid;
			position1 += op_length1;
			position = position1;
		} else if (i2 < mate2->core.n_cigar) {
			i2++;
			if (op_length2 == 0)
				continue;
			cigar_op = cigar_op2;
			contig = mate2->core.tid;
			position2 += op_length2;
			position = position2;
		} else {
			break; // all CIGAR elements of both mates have been processed
		}

		// increase coverage counter of windows that CIGAR element overlaps with
		if (bam_cigar_type(bam_cigar_op(cigar_op)) & 1/*consume query*/) {
			while (window <= position/COVERAGE_RESOLUTION) {
				if (coverage[contig][window] < USHRT_MAX)
					if (position - window * COVERAGE_RESOLUTION >= COVERAGE_RESOLUTION/2) // read must overlap at least half of the window
						coverage[contig][window]++;
				++window;
			}
		} else {
			window = position/COVERAGE_RESOLUTION;
		}

	}

	// store end of fragment
	if (!is_chimeric) { // the 'no_coverage' filter should only consider non-chimeric reads
		if ((mate1->core.flag & BAM_FREVERSE) || !(mate1->core.flag & BAM_FPAIRED))
			fragment_ends[mate1->core.tid][(position1-1)/COVERAGE_RESOLUTION] = true;
		else
			fragment_ends[mate2->core.tid][(position2-1)/COVERAGE_RESOLUTION] = true;
	}
}

// returns true, if a fragment begins at the given position
bool coverage_t::fragment_starts_here(const contig_t contig, const position_t start, const position_t end) const {
	if (contig >= fragment_starts.size())
		return false;
	for (unsigned int window = start/COVERAGE_RESOLUTION + 1; window <= end/COVERAGE_RESOLUTION; ++window) {
		if (window >= fragment_starts[contig].size())
			return false;
		if (fragment_starts[contig][window])
			return true;
	}
	return false;
}

// returns true, if a fragment ends at the given position
bool coverage_t::fragment_ends_here(const contig_t contig, const position_t start, const position_t end) const {
	if (contig >= fragment_ends.size())
		return false;
	for (unsigned int window = start/COVERAGE_RESOLUTION; window < end/COVERAGE_RESOLUTION; ++window) {
		if (window >= fragment_ends[contig].size())
			return false;
		if (fragment_ends[contig][window])
			return true;
	}
	return false;
}

// get coverage within a window of <COVERAGE_RESOLUTION> upstream or downstream of given position
int coverage_t::get_coverage(const contig_t contig, const position_t position, const direction_t direction) const {
	if (contig >= coverage.size() || coverage[contig].empty())
		return -1;
	if (direction == UPSTREAM) {
		if (position < COVERAGE_RESOLUTION)
			return 0;
		else
			return coverage[contig][position/COVERAGE_RESOLUTION-1];
	} else { // direction == DOWNSTREAM
		return coverage[contig][position/COVERAGE_RESOLUTION+1];
	}
}

