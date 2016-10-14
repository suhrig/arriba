#include <string>
#include <unordered_map>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "filter_nonexpressed.hpp"

using namespace std;

bool region_has_non_chimeric_reads(BGZF* bam_file, bam_index_t* bam_index, const contig_t contig, const position_t start, const position_t end, direction_t direction, chimeric_alignments_t& chimeric_alignments, bool exonic) {
	bam1_t* bam_record = bam_init1();
	bool result = false;

	bam_iter_t bam_record_iterator = bam_iter_query(bam_index, contig, start, end);
	while (bam_iter_read(bam_file, bam_record_iterator, bam_record) >= 0) {
		if (chimeric_alignments.find((char*) bam1_qname(bam_record)) == chimeric_alignments.end()) { // ignore chimeric reads
			if (exonic || // in the case of exonic breakpoints, we accept any overlapping read
				      // in the case of intronic/intergenic breakpoints, we are more strict:
				      // ignore reads which span OVER the given region due to splicing
				      // we only count those that start or end near the breakpoint

			    (direction == UPSTREAM   && !(bam_record->core.flag & BAM_FREVERSE) && bam_record->core.pos   >= start && bam_record->core.pos   <  end && bam_record->core.mpos                         >= start) ||
			    (direction == DOWNSTREAM &&  (bam_record->core.flag & BAM_FREVERSE) && bam_endpos(bam_record) >  start && bam_endpos(bam_record) <= end && bam_record->core.mpos+bam_record->core.l_qseq <= end  )) {
			                                                                                                                                            // ^quick and dirty way of calculating the end of the mate
				result = true;
				break;
			}
		}
	}

	bam_iter_destroy(bam_record_iterator);
	bam_destroy1(bam_record);

	return result;
}

unsigned int filter_nonexpressed(fusions_t& fusions, const string& bam_file_path, chimeric_alignments_t& chimeric_alignments, exon_annotation_t& exon_annotation) {

	// When the breakpoint is close to the start/end of a transcript, then
	// there might not be any reads in the normal rna.bam file around the breakpoint
	// other than the chimeric reads.
	// For this reason, we do not filter these breakpoints, when there is no expression
	// in their vicinity. But we need to determine which exons are the first/last
	// exons of a gene.
	unordered_map< gene_t,vector<exon_annotation_record_t*> > terminal_exons_by_gene;
	for (exon_annotation_t::iterator exon = exon_annotation.begin(); exon != exon_annotation.end(); ++exon) {

		if (terminal_exons_by_gene[exon->gene].size() == 0) { // we encounter the first exon of this gene => simply add it
			terminal_exons_by_gene[exon->gene].push_back(&(*exon));
		} else if (terminal_exons_by_gene[exon->gene].size() == 1) { // we encounter the second exon of this gene => simply add it
			if (exon->start < terminal_exons_by_gene[exon->gene][0]->start || exon->end > terminal_exons_by_gene[exon->gene][0]->end) { // ignore exons which are encompassed by other exons
				if (exon->start <= terminal_exons_by_gene[exon->gene][0]->start && exon->end >= terminal_exons_by_gene[exon->gene][0]->end) {
					terminal_exons_by_gene[exon->gene][0] = &(*exon);
				} else {
					terminal_exons_by_gene[exon->gene].push_back(&(*exon));
					if (terminal_exons_by_gene[exon->gene][0]->start > terminal_exons_by_gene[exon->gene][1]->start ||
					    terminal_exons_by_gene[exon->gene][0]->end > terminal_exons_by_gene[exon->gene][1]->end)
						swap(terminal_exons_by_gene[exon->gene][0], terminal_exons_by_gene[exon->gene][1]); // make sure exons are sorted by coordinate
				}
			}

		} else { // we encounter the third (or more) exon of this gene => check if it is a terminal exon
			if (terminal_exons_by_gene[exon->gene][0]->start > exon->start ||
			    terminal_exons_by_gene[exon->gene][0]->start == exon->start && terminal_exons_by_gene[exon->gene][0]->end < exon->end) {
				terminal_exons_by_gene[exon->gene][0] = &(*exon);
			} else if (terminal_exons_by_gene[exon->gene][1]->end < exon->end ||
			           terminal_exons_by_gene[exon->gene][1]->end == exon->end && terminal_exons_by_gene[exon->gene][1]->start > exon->start) {
				terminal_exons_by_gene[exon->gene][1] = &(*exon);
			}
		}

	}

	// open BAM file and load index
	BGZF* bam_file = bam_open(bam_file_path.c_str(), "rb");
	bam_index_t* bam_index = bam_index_load(bam_file_path.c_str());

	// for each fusion, check if there is any expression around the breakpoint
	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (!i->second.filters.empty())
			continue; // fusion has already been filtered


		if (!(i->second.contig1 == i->second.contig2 && i->second.breakpoint2 - i->second.breakpoint1 < 400000 && i->second.direction1 == DOWNSTREAM && i->second.direction2 == UPSTREAM)) { // fusion is not read-through

			if (i->second.split_reads1 + i->second.split_reads2 != 0 &&
			    i->second.split_reads1 + i->second.discordant_mates != 0 &&
			    i->second.split_reads2 + i->second.discordant_mates != 0) {
				++remaining;
				continue; // don't filter fusions with high support
			}

			if (i->second.spliced1 || i->second.spliced2) {
				++remaining;
				continue; // don't filter spliced breakpoints (they are more credible)
			}
		}

		position_t start, end;
		bool is_in_terminal_exon;

		// check if breakpoint1 is in a terminal exon
		is_in_terminal_exon = false;
		for (auto exon = terminal_exons_by_gene[i->second.gene1].begin(); exon != terminal_exons_by_gene[i->second.gene1].end() && !is_in_terminal_exon; ++exon)
			if (i->second.breakpoint1 >= (**exon).start && i->second.breakpoint1 <= (**exon).end)
				is_in_terminal_exon = true;

		if (!is_in_terminal_exon) {
			// check if there is coverage around breakpoint1
			if (i->second.direction1 == UPSTREAM) {
				start = i->second.breakpoint1;
				if (i->second.split_reads1 + i->second.split_reads2 == 0)
					start -= 200;
				end = max(i->second.breakpoint1 + 200, i->second.anchor_start1);
			} else {
				start = min(i->second.breakpoint1 - 200, i->second.anchor_start1);
				end = i->second.breakpoint1;
				if (i->second.split_reads1 + i->second.split_reads2 == 0)
					end += 200;
			}
			if (!region_has_non_chimeric_reads(bam_file, bam_index, i->second.contig1, start, end, i->second.direction1, chimeric_alignments, i->second.exonic1)) {
				i->second.filters.insert(FILTERS.at("non_expressed"));
				continue;
			}
		}

		// check if breakpoint2 is in a terminal exon
		is_in_terminal_exon = false;
		for (auto exon = terminal_exons_by_gene[i->second.gene2].begin(); exon != terminal_exons_by_gene[i->second.gene2].end() && !is_in_terminal_exon; ++exon)
			if (i->second.breakpoint2 >= (**exon).start && i->second.breakpoint2 <= (**exon).end)
				is_in_terminal_exon = true;

		if (!is_in_terminal_exon) {
			// check if there is coverage around breakpoint2
			if (i->second.direction2 == UPSTREAM) {
				start = i->second.breakpoint2;
				if (i->second.split_reads1 + i->second.split_reads2 == 0)
					start -= 200;
				end = max(i->second.breakpoint2 + 200, i->second.anchor_start2);
			} else {
				start = min(i->second.breakpoint2 - 200, i->second.anchor_start2);
				end = i->second.breakpoint2;
				if (i->second.split_reads1 + i->second.split_reads2 == 0)
					end += 200;
			}
			if (!region_has_non_chimeric_reads(bam_file, bam_index, i->second.contig2, start, end, i->second.direction2, chimeric_alignments, i->second.exonic2)) {
				i->second.filters.insert(FILTERS.at("non_expressed"));
				continue;
			}
		}

		remaining++;
	}

	// close BAM file and index
	bam_close(bam_file);
	bam_index_destroy(bam_index);

	return remaining;
}

