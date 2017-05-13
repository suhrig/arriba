#include <list>
#include "common.hpp"
#include "annotation.hpp"
#include "read_stats.hpp"

unsigned long int count_mapped_reads(const string& bam_file_path, const vector<bool>& interesting_contigs) {
	unsigned long int result = 0;
	bam_index_t* bam_index = bam_index_load(bam_file_path.c_str());
	for (unsigned int i = 0; i < interesting_contigs.size(); ++i) {
        	unsigned long int mapped, unmapped;
		hts_idx_get_stat(bam_index, i, &mapped, &unmapped);
		if (interesting_contigs[i]) // only count reads on interesting contigs
			result += mapped;
	}
	return result;
}

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
		if (1.0*within_range/mate_gaps.size() < 0.683 || no_more_outliers)
			break; // all outliers have been removed

		// remove outliers, if the mate gap distribution is not yet normally distributed
		no_more_outliers = true;
		for (list<int>::iterator i = mate_gaps.begin(); i != mate_gaps.end();) {
			if (*i < mate_gap_mean - 3*mate_gap_stddev || *i > mate_gap_mean + 3*mate_gap_stddev) {
				i = mate_gaps.erase(i); // remove outlier
				no_more_outliers = false;
			} else {
				++i;
			}
		}
	}

	return true;
}

