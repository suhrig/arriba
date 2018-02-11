#include <algorithm>
#include <iostream>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "extract_reads_read_through_fusion.hpp"

using namespace std;

void write_read_through_alignment(BGZF* bam_file, const bam1_t* const bam_record) {
	if (bam_write1(bam_file, bam_record) < 0) {
		cerr << "ERROR: failed to write BAM record to read-through fusions file." << endl;
		exit(1);
	}
}

bool find_spanning_intron(const bam1_t* read, const position_t gene1_end, const position_t gene2_start, unsigned int& cigar_op, position_t& read_pos, position_t& gene2_pos) {

	if (read->core.n_cigar < 3)
		return false; // there cannot be any introns, when there are less than 3 CIGAR operations

	position_t before_cigar_op = read->core.pos;
	position_t after_cigar_op;
	for (unsigned int i = 0; i < read->core.n_cigar; ++i) {
		unsigned int op_length = bam_cigar2rlen(1, bam1_cigar(read)+i);
		after_cigar_op = before_cigar_op + op_length;
		if ((bam1_cigar(read)[i] & 15) == BAM_CREF_SKIP && // this is an intron
		    (before_cigar_op <= gene1_end   && after_cigar_op >  gene1_end || // the intron spans over the end of gene1
		     before_cigar_op <  gene2_start && after_cigar_op >= gene2_start)) { // the intron spans over the start of gene2
			cigar_op = i;
			read_pos = bam_cigar2qlen(i, bam1_cigar(read));
			gene2_pos = after_cigar_op;
			return true;
		}
		before_cigar_op = after_cigar_op;
	}

	return false;
}

void clip_end(BGZF* read_through_fusions_file, bam1_t* read, unsigned int cigar_op, position_t read_pos, bool is_supplementary) {

	// clip read after intron given in cigar_op
	bam1_t* clipped_read = bam_dup1(read);
	if (clipped_read == NULL) {
		cerr << "ERROR: failed to allocate memory." << endl;
		exit(1);
	}
	copy( // throw away end of CIGAR string by shifting the end of the BAM record a bit to the front
		(char*) (bam1_cigar(read) + read->core.n_cigar),
		((char*) read->data) + read->data_len,
		(char*) (bam1_cigar(clipped_read) + cigar_op + 1)
	);
	clipped_read->core.n_cigar = cigar_op + 1; // correct length of CIGAR string after throwing away the end
	bam1_cigar(clipped_read)[cigar_op] = ((read->core.l_qseq - read_pos)<<4) + BAM_CSOFT_CLIP; // soft-clip spliced part of read
	clipped_read->data_len -= (read->core.n_cigar - cigar_op - 1) * 4; // correct length of variable data of BAM record
	if (is_supplementary)
		clipped_read->core.flag |= BAM_FSECONDARY; // mark supplementary read as secondary alignment
	write_read_through_alignment(read_through_fusions_file, clipped_read); // write to output
	bam_destroy1(clipped_read); // free memory
}

void clip_start(BGZF* read_through_fusions_file, bam1_t* read, unsigned int cigar_op, position_t read_pos, position_t gene2_pos, bool is_supplementary) {

	// clip read before intron given in cigar_op
	bam1_t* clipped_read = bam_dup1(read);
	if (clipped_read == NULL) {
		cerr << "ERROR: failed to allocate memory." << endl;
		exit(1);
	}
	for (unsigned int i = cigar_op+1; i < read->core.n_cigar; ++i) // copy end of CIGAR string to beginning
		bam1_cigar(clipped_read)[i - cigar_op] = bam1_cigar(clipped_read)[i];
	copy( // throw away end of CIGAR string by shifting the end of the BAM record a bit to the front
		(char*) (bam1_cigar(read) + read->core.n_cigar),
		((char*) read->data) + read->data_len,
		(char*) (bam1_cigar(clipped_read) + read->core.n_cigar - cigar_op)
	);
	clipped_read->core.n_cigar = read->core.n_cigar - cigar_op; // correct length of CIGAR string after throwing away the end
	bam1_cigar(clipped_read)[0] = (read_pos<<4) + BAM_CSOFT_CLIP; // soft-clip spliced part of read
	clipped_read->data_len -= cigar_op * 4; // correct length of variable data of BAM record
	clipped_read->core.pos = gene2_pos; // set start of split read to coordinate in 2nd gene
	if (is_supplementary)
		clipped_read->core.flag |= BAM_FSECONDARY; // mark supplementary read as secondary alignment
	write_read_through_alignment(read_through_fusions_file, clipped_read); // write to output
	bam_destroy1(clipped_read); // free memory
}

void extract_read_through_fusion(BGZF* read_through_fusions_file, bam1_t* forward_mate, bam1_t* reverse_mate, const gene_annotation_index_t& gene_annotation_index) {

	if (forward_mate->core.flag & BAM_FUNMAP) // ignore unmapped reads
		return;

	if (forward_mate->core.flag & BAM_FPAIRED) // paired-end data
		if ((reverse_mate->core.flag & BAM_FUNMAP) || // ignore unmapped reads
		    !(forward_mate->core.flag & BAM_FPROPER_PAIR)) // ignore discordant mates, they are in the chimeric.bam file already
		return;

	// find out which read is on the forward strand and which on the reverse
	if (forward_mate->core.flag & BAM_FREVERSE)
		swap(forward_mate, reverse_mate);

	// check if one mate maps inside the gene and the other outside
	gene_set_t forward_mate_genes, reverse_mate_genes;
	if (forward_mate != NULL)
		get_annotation_by_coordinate(forward_mate->core.tid, forward_mate->core.pos, forward_mate->core.pos, forward_mate_genes, gene_annotation_index);
	else
		get_annotation_by_coordinate(reverse_mate->core.tid, reverse_mate->core.pos, reverse_mate->core.pos, forward_mate_genes, gene_annotation_index);
	if (reverse_mate != NULL)
		get_annotation_by_coordinate(reverse_mate->core.tid, bam_endpos(reverse_mate), bam_endpos(reverse_mate), reverse_mate_genes, gene_annotation_index);
	else
		get_annotation_by_coordinate(forward_mate->core.tid, bam_endpos(forward_mate), bam_endpos(forward_mate), reverse_mate_genes, gene_annotation_index);
	gene_set_t common_genes;
	combine_annotations(forward_mate_genes, reverse_mate_genes, common_genes, false);
	if (common_genes.empty() && !(forward_mate_genes.empty() && reverse_mate_genes.empty())) { // mate1 and mate2 map to different genes => potential read-through fusion

		// there are three possibilites to get here:
		// (1) we have a split read (one mate has an intron and part of this mate maps inside the gene and the other part outside)
		// (2) we have a discordant mate (one mate maps far away from the gene of the other mate)
		// (3) we have normally mapped mates, but the end of one of them reaches just outside the gene => ignore those

		// find the biggest genes that the mates overlap with
		position_t forward_gene_start, forward_gene_end, reverse_gene_start, reverse_gene_end;
		get_boundaries_of_biggest_gene(forward_mate_genes, forward_gene_start, forward_gene_end);
		get_boundaries_of_biggest_gene(reverse_mate_genes, reverse_gene_start, reverse_gene_end);

		// if a mate does not overlap with any gene,
		// use the boundaries of the mate or of the gene of the other mate
		if (forward_gene_end == -1)
			forward_gene_end = reverse_gene_start - 1;
		if (reverse_gene_start == -1)
			reverse_gene_start = forward_gene_end + 1;

		// check for possibility (1) => find intron in CIGAR string which spans the fusion genes
		unsigned int forward_cigar_op, reverse_cigar_op;
		unsigned int forward_matches_before_intron, reverse_matches_before_intron;
		position_t forward_read_pos, reverse_read_pos;
		position_t forward_gene2_pos, reverse_gene2_pos;
		bool forward_mate_has_intron = (forward_mate == NULL) ? false : find_spanning_intron(forward_mate, forward_gene_end, reverse_gene_start, forward_cigar_op, forward_read_pos, forward_gene2_pos);
		bool reverse_mate_has_intron = (reverse_mate == NULL) ? false : find_spanning_intron(reverse_mate, forward_gene_end, reverse_gene_start, reverse_cigar_op, reverse_read_pos, reverse_gene2_pos);
		if (forward_mate_has_intron &&
		    (!reverse_mate_has_intron || forward_read_pos < reverse_mate->core.l_qseq - reverse_read_pos)) { // if both mates are clipped, use the one with the longer segment as anchor

			// make split read and supplementary from forward mate
			clip_start(read_through_fusions_file, forward_mate, forward_cigar_op, forward_read_pos, forward_gene2_pos, false); // split read
			clip_end(read_through_fusions_file, forward_mate, forward_cigar_op, forward_read_pos, true); // supplementary
					
			if (reverse_mate != NULL) {
				// set mpos of reverse mate to start of primary alignment of forward mate
				reverse_mate->core.mpos = forward_gene2_pos;
				if (reverse_mate_has_intron) { // reverse mate overlaps with breakpoint => clip it
					clip_start(read_through_fusions_file, reverse_mate, reverse_cigar_op, reverse_read_pos, reverse_gene2_pos, false);
				} else { // reverse mate overlaps with forward mate, but not with breakpoint
					// write reverse mate to output as is
					write_read_through_alignment(read_through_fusions_file, reverse_mate);
				}
			}

		} else if (reverse_mate_has_intron) {

			// make split read and supplementary from reverse mate
			clip_start(read_through_fusions_file, reverse_mate, reverse_cigar_op, reverse_read_pos, reverse_gene2_pos, true); // supplementary
			clip_end(read_through_fusions_file, reverse_mate, reverse_cigar_op, reverse_read_pos, false); // split read

			if (forward_mate != NULL) {
				if (forward_mate_has_intron) { // forward mate overlaps with breakpoints => clip it
					clip_end(read_through_fusions_file, forward_mate, forward_cigar_op, forward_read_pos, false);
				} else { // forward mate overlaps with reverse mate, but not with breakpoint
					// write forward mate to output as is
					write_read_through_alignment(read_through_fusions_file, forward_mate);
				}
			}

		// check for possibility (2) => the mates must be contained within different genes
		} else if (forward_mate != NULL && reverse_mate != NULL &&
		           reverse_mate->core.pos   >= reverse_gene_start &&
		           bam_endpos(forward_mate) <= forward_gene_end) {

			// write discordant mates to output file
			write_read_through_alignment(read_through_fusions_file, forward_mate);
			write_read_through_alignment(read_through_fusions_file, reverse_mate);

		} // else possibility (3)
	}
}

