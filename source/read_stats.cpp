#include <list>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "read_stats.hpp"

void estimate_mate_gap_distribution(const string& bam_file_path, float& mate_gap_mean, float& mate_gap_stddev, const gene_annotation_index_t& gene_annotation_index, const exon_annotation_index_t& exon_annotation_index, const unsigned int sample_size) {

	// open BAM file
	BGZF* bam_file = bam_open(bam_file_path.c_str(), "rb");
	bam_header_t* bam_header = bam_header_read(bam_file);

	// cycle through BAM file and note down mate gaps
	bam1_t* bam_record = bam_init1();
	unsigned int count = 0;
	list<int> mate_gaps;
	while (
	       count <= sample_size &&		   // read at most <sample_size> reads
	       bam_read1(bam_file, bam_record) > 0 // or until the end of the file is reached
	      ) {

		// ignore secondary alignments and mates which aren't paired properly
		// also, we only use mates on the forward strand in order to not count the same mate gap twice
		if ((bam_record->core.flag & BAM_FPROPER_PAIR) && !(bam_record->core.flag & BAM_FSECONDARY) && !(bam_record->core.flag & BAM_FREVERSE)) {

			position_t mate1_end = bam_endpos(bam_record);

			// check if both mates map to exons of the same gene
			// (if not then the alignment might be strange and mate gap calculation would be skewed)
			gene_set_t mate1_genes, mate2_genes, combined_genes;
			get_annotation_by_coordinate(bam_record->core.tid, bam_record->core.pos, mate1_end, mate1_genes, gene_annotation_index);
			get_annotation_by_coordinate(bam_record->core.mtid, bam_record->core.mpos, bam_record->core.mpos, mate2_genes, gene_annotation_index);
			combine_annotations(mate1_genes, mate2_genes, combined_genes, false);

			if (combined_genes.size() == 1 && // both mates map to exons of one (and only one) gene
			    !(**combined_genes.begin()).is_dummy) { // ignore intergenic regions
				mate_gaps.push_back(get_spliced_distance(bam_record->core.tid, mate1_end, bam_record->core.mpos, DOWNSTREAM, UPSTREAM, *combined_genes.begin(), exon_annotation_index));
				count++;
			}
		}
	}

	// close BAM file
	bam_destroy1(bam_record);
	bam_header_destroy(bam_header);
	bam_close(bam_file);

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
}

