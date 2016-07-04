#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <string>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "options.hpp"
#include "options_extract_read-through_fusions.hpp"

using namespace std;

bool find_spanning_intron(const bam1_t* read, const position_t gene1_end, const position_t gene2_start, unsigned int& cigar_op, position_t& read_pos, position_t& gene2_pos) {

	if (read->core.n_cigar < 3)
		return false; // there cannot be any introns, when there are less than 3 CIGAR operations

	bool found = false; // return value of function

	position_t before_cigar_op = read->core.pos;
	position_t after_cigar_op;
	for (unsigned int i = 0; i < read->core.n_cigar; ++i) {
		after_cigar_op = before_cigar_op + bam_cigar2rlen(1, bam1_cigar(read)+i);
		if ((bam1_cigar(read)[i] & BAM_CREF_SKIP) == BAM_CREF_SKIP && // this is an intron
		    (before_cigar_op <= gene1_end   && after_cigar_op >  gene1_end || // the intron spans over the end of gene1
		     before_cigar_op <  gene2_start && after_cigar_op >= gene2_start)) { // the intron spans over the start of gene2
			found = true;
			cigar_op = i;
			read_pos = bam_cigar2qlen(i, bam1_cigar(read));
			gene2_pos = after_cigar_op;
			break;
		}
		before_cigar_op = after_cigar_op;
	}

	return found;
}

void clip_end(BGZF* output_bam_file, bam1_t* read, unsigned int cigar_op, position_t read_pos, bool is_supplementary) {

	// clip read after intron given in cigar_op
	bam1_t* clipped_read = bam_dup1(read);
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
	bam_write1(output_bam_file, clipped_read); // write to output
	bam_destroy1(clipped_read); // free memory
}

void clip_start(BGZF* output_bam_file, bam1_t* read, unsigned int cigar_op, position_t read_pos, position_t gene2_pos, bool is_supplementary) {

	// clip read before intron given in cigar_op
	bam1_t* clipped_read = bam_dup1(read);
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
	bam_write1(output_bam_file, clipped_read); // write to output
	bam_destroy1(clipped_read); // free memory
}

typedef map<string,bam1_t*> buffered_bam_records_t;

int main(int argc, char **argv) {

	options_t options = parse_arguments(argc, argv);

	// open input BAM file
	BGZF* input_bam_file = bam_open(options.input_bam_file.c_str(), "rb");
	bam_header_t* bam_header = bam_header_read(input_bam_file);

	// make a mapping of contig name -> contig ID
	contigs_t contigs;
	for (int32_t i = 0; i < bam_header->n_targets; ++i)
		contigs[removeChr(bam_header->target_name[i])] = i;

	// read gene annotation
	annotation_t gene_annotation;
	read_annotation_bed(options.gene_annotation_file, gene_annotation, contigs);
	annotation_index_t gene_annotation_index;
	make_annotation_index(gene_annotation, gene_annotation_index, contigs);

	// open output BAM file
	BGZF* output_bam_file = bam_open(options.output_bam_file.c_str(), "w");
	bam_header_write(output_bam_file, bam_header);

	// read BAM records
	bam1_t* bam_record = bam_init1();
	buffered_bam_records_t buffered_bam_records; // holds the first mate until we have found the second
	while (bam_read1(input_bam_file, bam_record) > 0) {

		if (!(bam_record->core.flag & BAM_FPROPER_PAIR) || // ignore discordant mates, they are in the chimeric.bam file already
		    (bam_record->core.flag & BAM_FUNMAP) || (bam_record->core.flag & BAM_FMUNMAP) || // ignore single mates
		    (bam_record->core.flag & BAM_FSECONDARY)) // ignore secondary alignments
			continue;
	
		// try to insert the mate into the buffered BAM records
		// if there was already a record with the same read name, insertion will fail (set to false) and
		// existing_element will point to the mate which was already in the buffered BAM records
		pair<buffered_bam_records_t::iterator,bool> existing_element = buffered_bam_records.insert(pair<string,bam1_t*>((char*) bam1_qname(bam_record), bam_record));

		if (get<1>(existing_element)) { // insertion succeeded, i.e., this is the first mate with the given read name, which we encounter
			
			bam_record = bam_init1(); // allocate memory for the next record

		} else { // insertion failed, i.e., we have already read the first mate previously

			// get the mate that already existed in the buffered BAM records
			buffered_bam_records_t::iterator buffered_bam_record = get<0>(existing_element);

			// find out which read is on the forward strand and which on the reverse
			bam1_t* forward_mate = buffered_bam_record->second;
			bam1_t* reverse_mate = bam_record;
			if (forward_mate->core.flag & BAM_FREVERSE)
				swap(forward_mate, reverse_mate);

			// check if one mate maps inside the gene and the other outside
			gene_set_t forward_mate_genes, reverse_mate_genes;
			get_annotation_by_coordinate(forward_mate->core.tid, forward_mate->core.pos, forward_mate->core.pos, forward_mate_genes, gene_annotation_index);
			get_annotation_by_coordinate(bam_record->core.tid, bam_endpos(reverse_mate), bam_endpos(reverse_mate), reverse_mate_genes, gene_annotation_index);
			gene_set_t common_genes;
			combine_annotations(forward_mate_genes, reverse_mate_genes, common_genes, false);
			if (common_genes.empty() && !(forward_mate_genes.empty() && reverse_mate_genes.empty())) { // mate1 and mate2 map to different genes => potential read-through fusion

				// there are three possibilites to get here:
				// (1) we have a split read (one mate has an intron and part of this mate maps inside the gene and the other part outside)
				// (2) we have a discordant mate (one mate maps far away from the gene of the other mate)
				// (3) we have normally mapped mates, but the end of one of them reaches just outside the gene => ignore those

				// find the biggest genes that the mates overlap with
				position_t forward_gene_start, forward_gene_end, reverse_gene_start, reverse_gene_end;
				get_boundaries_of_biggest_gene(forward_mate_genes, gene_annotation, forward_gene_start, forward_gene_end);
				get_boundaries_of_biggest_gene(reverse_mate_genes, gene_annotation, reverse_gene_start, reverse_gene_end);

				// if a mate does not overlap with any gene,
				// use the boundaries of the mate or of the gene of the other mate
				if (forward_gene_start == -1)
					forward_gene_start = forward_mate->core.pos;
				if (forward_gene_end == -1)
					forward_gene_end = reverse_gene_start - 1;
				if (reverse_gene_start == -1)
					reverse_gene_start = forward_gene_end + 1;
				if (reverse_gene_end == -1)
					reverse_gene_end = bam_endpos(reverse_mate);

				// check for possibility (1) => find intron in CIGAR string which spans the fusion genes
				unsigned int cigar_op;
				position_t read_pos;
				position_t gene2_pos;
				if (find_spanning_intron(forward_mate, forward_gene_end, reverse_gene_start, cigar_op, read_pos, gene2_pos)) {

					// make split read and supplementary from forward mate
					clip_start(output_bam_file, forward_mate, cigar_op, read_pos, gene2_pos, false); // split read
					clip_end(output_bam_file, forward_mate, cigar_op, read_pos, true); // supplementary
					
					if (reverse_mate->core.pos >= bam_endpos(forward_mate)) {
						// write reverse mate to output as is
						bam_write1(output_bam_file, reverse_mate);
					} else { // reverse mate overlaps with forward mate
						if (find_spanning_intron(reverse_mate, forward_gene_end, reverse_gene_start, cigar_op, read_pos, gene2_pos)) { // reverse mate overlaps with breakpoint => clip it
							clip_start(output_bam_file, reverse_mate, cigar_op, read_pos, gene2_pos, false);
						} else { // reverse mate overlaps with forward mate, but not with breakpoint
							// write reverse mate to output as is
							bam_write1(output_bam_file, reverse_mate);
						}
					}

				} else if (find_spanning_intron(reverse_mate, forward_gene_end, reverse_gene_start, cigar_op, read_pos, gene2_pos)) {

					// make split read and supplementary from reverse mate
					clip_start(output_bam_file, reverse_mate, cigar_op, read_pos, gene2_pos, true); // supplementary
					clip_end(output_bam_file, reverse_mate, cigar_op, read_pos, false); // split read

					if (bam_endpos(forward_mate) <= reverse_mate->core.pos) {
						// write forward mate to output as is
						bam_write1(output_bam_file, forward_mate);
					} else { // forward mate overlaps with reverse mate
						if (find_spanning_intron(forward_mate, forward_gene_end, reverse_gene_start, cigar_op, read_pos, gene2_pos)) { // forward mate overlaps with breakpoints => clip it
							clip_end(output_bam_file, forward_mate, cigar_op, read_pos, false);
						} else { // forward mate overlaps with reverse mate, but not with breakpoint
							// write forward mate to output as is
							bam_write1(output_bam_file, forward_mate);
						}
					}

				} else
				// check for possibility (2) => the mates must be contained within different genes
				if (reverse_mate->core.pos   >= reverse_gene_start &&
				    bam_endpos(forward_mate) <= forward_gene_end) {

					// write discordant mates to output file
					bam_write1(output_bam_file, forward_mate);
					bam_write1(output_bam_file, reverse_mate);

				} // else possibility (3)
			}

			// free memory
			bam_destroy1(buffered_bam_record->second);
			buffered_bam_records.erase(buffered_bam_record);
		}

	}

	// close BAM files
	bam_destroy1(bam_record);
	bam_header_destroy(bam_header);
	bam_close(input_bam_file);
	bam_close(output_bam_file);

	return 0;
}

