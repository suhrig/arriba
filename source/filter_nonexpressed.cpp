#include <string>
#include <unordered_map>
#include "cram.h"
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "filter_nonexpressed.hpp"

using namespace std;

bool region_has_non_chimeric_reads(samFile* bam_file, hts_idx_t* bam_index, const contig_t contig, const position_t start, const position_t end, const direction_t direction, const chimeric_alignments_t& chimeric_alignments, const bool exonic) {
	bam1_t* bam_record = bam_init1();
	bool result = false;

	hts_itr_t* bam_record_iterator = sam_itr_queryi(bam_index, contig, start, end);
	while (sam_itr_next(bam_file, bam_record_iterator, bam_record) >= 0) {
		if (chimeric_alignments.find((char*) bam_get_qname(bam_record)) == chimeric_alignments.end()) { // ignore chimeric reads
			if (exonic || // in the case of exonic breakpoints, we accept any overlapping read
				      // in the case of intronic/intergenic breakpoints, we are more strict:
				      // ignore reads which span OVER the given region due to splicing
				      // we only count those that start or end near the breakpoint

			    (direction == UPSTREAM   && !(bam_record->core.flag & BAM_FREVERSE) && bam_record->core.pos   >= start && bam_record->core.pos   <  end && (bam_record->core.mpos                         >= start || !(bam_record->core.flag & BAM_FPAIRED))) ||
			    (direction == DOWNSTREAM &&  (bam_record->core.flag & BAM_FREVERSE) && bam_endpos(bam_record) >  start && bam_endpos(bam_record) <= end && (bam_record->core.mpos+bam_record->core.l_qseq <= end   || !(bam_record->core.flag & BAM_FPAIRED)))) {
			                                                                                                                                            // ^quick and dirty way of calculating the end of the mate
				result = true;
				break;
			}
		}
	}

	sam_itr_destroy(bam_record_iterator);
	bam_destroy1(bam_record);

	return result;
}

unsigned int filter_nonexpressed(fusions_t& fusions, const string& bam_file_path, const string& assembly_file_path, const chimeric_alignments_t& chimeric_alignments, const exon_annotation_index_t& exon_annotation_index, const int max_mate_gap) {

	// open BAM file and load index
	samFile* bam_file = sam_open(bam_file_path.c_str(), "rb");
	if (bam_file->is_cram)
		cram_set_option(bam_file->fp.cram, CRAM_OPT_REFERENCE, assembly_file_path.c_str());
	hts_idx_t* bam_index = sam_index_load(bam_file, bam_file_path.c_str());

	// for each fusion, check if there is any expression around the breakpoint
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != NULL)
			continue; // fusion has already been filtered

		if (!fusion->second.is_read_through()) {

			if (fusion->second.split_reads1 + fusion->second.split_reads2 != 0 &&
			    fusion->second.split_reads1 + fusion->second.discordant_mates != 0 &&
			    fusion->second.split_reads2 + fusion->second.discordant_mates != 0) {
				++remaining;
				continue; // don't filter fusions with high support
			}

			if (fusion->second.spliced1 || fusion->second.spliced2) {
				++remaining;
				continue; // don't filter spliced breakpoints (they are more credible)
			}
		}

		position_t start, end;
		bool is_in_terminal_exon;

		// check if breakpoint1 is in a terminal exon
		exon_set_t exons;
		get_annotation_by_coordinate(fusion->second.contig1, fusion->second.breakpoint1, fusion->second.breakpoint1, exons, exon_annotation_index);
		is_in_terminal_exon = false;
		for (auto exon = exons.begin(); exon != exons.end() && !is_in_terminal_exon; ++exon)
			if ((**exon).gene == fusion->second.gene1 && ((**exon).is_transcript_start || (**exon).is_transcript_end))
				is_in_terminal_exon = true;

		if (!is_in_terminal_exon) {
			// check if there is coverage around breakpoint1
			if (fusion->second.direction1 == UPSTREAM) {
				start = fusion->second.breakpoint1;
				if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0)
					start -= max_mate_gap;
				end = max(fusion->second.breakpoint1 + max_mate_gap, fusion->second.anchor_start1);
			} else {
				start = min(fusion->second.breakpoint1 - max_mate_gap, fusion->second.anchor_start1);
				end = fusion->second.breakpoint1;
				if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0)
					end += max_mate_gap;
			}
			if (!region_has_non_chimeric_reads(bam_file, bam_index, fusion->second.contig1, start, end, fusion->second.direction1, chimeric_alignments, fusion->second.exonic1)) {
				fusion->second.filter = FILTERS.at("non_expressed");
				continue;
			}
		}

		// check if breakpoint2 is in a terminal exon
		exons.clear();
		get_annotation_by_coordinate(fusion->second.contig2, fusion->second.breakpoint2, fusion->second.breakpoint2, exons, exon_annotation_index);
		is_in_terminal_exon = false;
		for (auto exon = exons.begin(); exon != exons.end() && !is_in_terminal_exon; ++exon)
			if ((**exon).gene == fusion->second.gene2 && ((**exon).is_transcript_start || (**exon).is_transcript_end))
				is_in_terminal_exon = true;

		if (!is_in_terminal_exon) {
			// check if there is coverage around breakpoint2
			if (fusion->second.direction2 == UPSTREAM) {
				start = fusion->second.breakpoint2;
				if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0)
					start -= max_mate_gap;
				end = max(fusion->second.breakpoint2 + max_mate_gap, fusion->second.anchor_start2);
			} else {
				start = min(fusion->second.breakpoint2 - max_mate_gap, fusion->second.anchor_start2);
				end = fusion->second.breakpoint2;
				if (fusion->second.split_reads1 + fusion->second.split_reads2 == 0)
					end += max_mate_gap;
			}
			if (!region_has_non_chimeric_reads(bam_file, bam_index, fusion->second.contig2, start, end, fusion->second.direction2, chimeric_alignments, fusion->second.exonic2)) {
				fusion->second.filter = FILTERS.at("non_expressed");
				continue;
			}
		}

		remaining++;
	}

	// close BAM file and index
	hts_idx_destroy(bam_index);
	sam_close(bam_file);

	return remaining;
}

