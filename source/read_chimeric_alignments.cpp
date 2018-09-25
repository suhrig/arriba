#include <algorithm>
#include <iostream>
#include <string>
#include "cram.h"
#include "sam.h"
#include "annotation.hpp"
#include "common.hpp"
#include "read_chimeric_alignments.hpp"
#include "read_stats.hpp"

using namespace std;

typedef map<string,bam1_t*> buffered_bam_records_t;
typedef vector<contig_t> tid_to_contig_t;

bool find_spanning_intron(const bam1_t* bam_record, const position_t gene1_end, const position_t gene2_start, unsigned int& cigar_op, position_t& read_pos) {

	if (bam_record->core.n_cigar < 3)
		return false; // there cannot be any introns, when there are less than 3 CIGAR operations

	position_t before_cigar_op = bam_record->core.pos;
	position_t after_cigar_op;
	for (unsigned int i = 0; i < bam_record->core.n_cigar; ++i) {
		uint32_t current_op = bam_get_cigar(bam_record)[i];
		unsigned int op_length = (bam_cigar_type(bam_cigar_op(current_op)) & 2) ? bam_cigar_oplen(current_op) : 0;
		after_cigar_op = before_cigar_op + op_length;
		if ((current_op & 15) == BAM_CREF_SKIP && // this is an intron
		    (before_cigar_op <= gene1_end   && after_cigar_op >  gene1_end || // the intron spans over the end of gene1
		     before_cigar_op <  gene2_start && after_cigar_op >= gene2_start)) { // the intron spans over the start of gene2
			cigar_op = i;
			read_pos = bam_cigar2qlen(i, bam_get_cigar(bam_record));
			return true;
		}
		before_cigar_op = after_cigar_op;
	}

	return false;
}

void add_chimeric_alignment(chimeric_alignments_t& chimeric_alignments, const bam1_t* bam_record, unsigned int cigar_op = 0, const position_t read_pos = 0, const bool clip_start = false, const bool clip_end = false, const bool is_supplementary = false) {

	// convert bam1_t structure into our own structure and discard information we don't need
	string name = (char*) bam_get_qname(bam_record);
	mates_t* mates = &chimeric_alignments[name];
	mates->single_end = !(bam_record->core.flag & BAM_FPAIRED);
	mates->resize(mates->size()+1);
	alignment_t& alignment = (*mates)[mates->size()-1];
	alignment.strand = (bam_record->core.flag & BAM_FREVERSE) ? REVERSE : FORWARD;
	alignment.first_in_pair = bam_record->core.flag & BAM_FREAD1;
	alignment.contig = bam_record->core.tid;
	alignment.supplementary = is_supplementary;
	if (!is_supplementary) { // only keep sequence in memory, if this is not the supplementary alignment (because then it's already stored in the split-read)
		alignment.sequence.resize(bam_record->core.l_qseq);
		for (unsigned int i = 0; i < bam_record->core.l_qseq; ++i)
			alignment.sequence[i] = seq_nt16_str[bam_seqi(bam_get_seq(bam_record), i)];
	}

	// read-through alignments need to be split into a split-read and a supplementary alignment
	if (clip_start) {
		alignment.start = bam_record->core.pos + bam_cigar2rlen(cigar_op + 1, bam_get_cigar(bam_record));
		alignment.end = bam_endpos(bam_record) - 1;
		alignment.cigar.resize(bam_record->core.n_cigar - cigar_op);
		alignment.cigar[0] = (read_pos<<4) + BAM_CSOFT_CLIP; // soft-clip start of read
		for (unsigned int i = cigar_op+1; i < bam_record->core.n_cigar; ++i) // copy cigar operations from <cigar_op> onwards
			alignment.cigar[i-cigar_op] = bam_get_cigar(bam_record)[i];
	} else if (clip_end) {
		alignment.start = bam_record->core.pos;
		alignment.end = bam_record->core.pos + bam_cigar2rlen(cigar_op, bam_get_cigar(bam_record)) - 1;
		alignment.cigar.resize(cigar_op + 1);
		for (unsigned int i = 0; i < cigar_op; ++i) // copy cigar operations up until <cigar_op>
			alignment.cigar[i] = bam_get_cigar(bam_record)[i];
		alignment.cigar[cigar_op] = ((bam_record->core.l_qseq - read_pos)<<4) + BAM_CSOFT_CLIP; // soft-clip end of read
	} else { // do not clip, i.e., alignment is already split into a split-read and a supplementary alignment
		alignment.start = bam_record->core.pos;
		alignment.end = bam_endpos(bam_record) - 1;
		alignment.cigar.resize(bam_record->core.n_cigar);
		for (unsigned int i = 0; i < alignment.cigar.size(); ++i)
			alignment.cigar[i] = bam_get_cigar(bam_record)[i];
	}
}

bool extract_read_through_alignment(chimeric_alignments_t& chimeric_alignments, bam1_t* forward_mate, bam1_t* reverse_mate, const gene_annotation_index_t& gene_annotation_index, const bool separate_chimeric_bam_file) {

	if (forward_mate->core.flag & BAM_FUNMAP) // ignore unmapped reads
		return false;

	if (forward_mate->core.flag & BAM_FPAIRED) // paired-end data
		if ((reverse_mate->core.flag & BAM_FUNMAP) || // ignore unmapped reads
		    !(forward_mate->core.flag & BAM_FPROPER_PAIR)) // ignore discordant mates, they are in the chimeric.bam file already
		return false;

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
		bool forward_mate_has_intron = (forward_mate == NULL) ? false : find_spanning_intron(forward_mate, forward_gene_end, reverse_gene_start, forward_cigar_op, forward_read_pos);
		bool reverse_mate_has_intron = (reverse_mate == NULL) ? false : find_spanning_intron(reverse_mate, forward_gene_end, reverse_gene_start, reverse_cigar_op, reverse_read_pos);
		if (forward_mate_has_intron &&
		    (!reverse_mate_has_intron || forward_read_pos < reverse_mate->core.l_qseq - reverse_read_pos)) { // if both mates are clipped, use the one with the longer segment as anchor

			// only add the read-through alignment, if it is not also a chimeric alignment
			if (!separate_chimeric_bam_file || chimeric_alignments.find((char*) bam_get_qname(forward_mate)) == chimeric_alignments.end()) {
				// make split read and supplementary from forward mate
				add_chimeric_alignment(chimeric_alignments, forward_mate, forward_cigar_op, forward_read_pos, true/*clip start*/, false, false/*split-read*/);
				add_chimeric_alignment(chimeric_alignments, forward_mate, forward_cigar_op, forward_read_pos, false, true/*clip end*/, true/*supplementary*/);

				if (reverse_mate != NULL) { // paired-end
					if (reverse_mate_has_intron) // reverse mate overlaps with breakpoint => clip it
						add_chimeric_alignment(chimeric_alignments, reverse_mate, reverse_cigar_op, reverse_read_pos, true, false, false);
					else // reverse mate overlaps with forward mate, but not with breakpoint => add it as is
						add_chimeric_alignment(chimeric_alignments, reverse_mate);
				}
				return true;
			}

		} else if (reverse_mate_has_intron) {

			// only add the read-through alignment, if it is not also a chimeric alignment
			if (!separate_chimeric_bam_file || chimeric_alignments.find((char*) bam_get_qname(reverse_mate)) == chimeric_alignments.end()) {
				// make split read and supplementary from reverse mate
				add_chimeric_alignment(chimeric_alignments, reverse_mate, reverse_cigar_op, reverse_read_pos, true/*clip start*/, false, true/*supplementary*/);
				add_chimeric_alignment(chimeric_alignments, reverse_mate, reverse_cigar_op, reverse_read_pos, false, true/*clip end*/, false/*split-read*/);

				if (forward_mate != NULL) { // paired-end
					if (forward_mate_has_intron) // forward mate overlaps with breakpoints => clip it at the end
						add_chimeric_alignment(chimeric_alignments, forward_mate, forward_cigar_op, forward_read_pos, false, true, false);
					else // forward mate overlaps with reverse mate, but not with breakpoint => add it as is
						add_chimeric_alignment(chimeric_alignments, forward_mate);
				}
				return true;
			}

		// check for possibility (2) => the mates must be contained within different genes
		} else if (forward_mate != NULL && reverse_mate != NULL && // paired-end
		           reverse_mate->core.pos   >= reverse_gene_start &&
		           bam_endpos(forward_mate) <= forward_gene_end) {

			// only add the read-through alignment, if it is not also a chimeric alignment
			if (!separate_chimeric_bam_file || chimeric_alignments.find((char*) bam_get_qname(forward_mate)) == chimeric_alignments.end()) {
				// add discordant mates to chimeric alignments file
				add_chimeric_alignment(chimeric_alignments, forward_mate);
				add_chimeric_alignment(chimeric_alignments, reverse_mate);
				return true;
			}

		} // else possibility (3)
	}
	return false;
}

unsigned int read_chimeric_alignments(const string& bam_file_path, const string& assembly_file_path, chimeric_alignments_t& chimeric_alignments, unsigned long int& mapped_reads, coverage_t& coverage, contigs_t& contigs, const contigs_t& interesting_contigs, const gene_annotation_index_t& gene_annotation_index, const bool separate_chimeric_bam_file, const bool is_rna_bam_file) {

	// open BAM file
	samFile* bam_file = sam_open(bam_file_path.c_str(), "rb");
	if (bam_file->is_cram)
		cram_set_option(bam_file->fp.cram, CRAM_OPT_REFERENCE, assembly_file_path.c_str());
	bam_hdr_t* bam_header = sam_hdr_read(bam_file);

	// add contigs which are not yet listed in <contigs>
	// and make a map tid -> contig, because the contig IDs in the BAM file need not necessarily match the contig IDs in the GTF file
	tid_to_contig_t tid_to_contig(bam_header->n_targets);
	vector<bool> interesting_tids(bam_header->n_targets);
	for (uint32_t target = 0; target < bam_header->n_targets; ++target) {
		string contig_name = removeChr(bam_header->target_name[target]);
		contigs.insert(pair<string,contig_t>(contig_name, contigs.size())); // this fails (i.e., nothing is inserted), if the contig already exists
		tid_to_contig[target] = contigs[contig_name];
		if (is_rna_bam_file) // only count reads of Aligned.out.bam, not of Chimeric.out.sam
			interesting_tids[target] = (interesting_contigs.find(contig_name) != interesting_contigs.end()) || interesting_contigs.empty();
	}

	// read BAM records
	bam1_t* bam_record = bam_init1();
	if (bam_record == NULL) {
		cerr << "ERROR: failed to allocate memory." << endl;
		exit(1);
	}
	buffered_bam_records_t buffered_bam_records; // holds the first mate until we have found the second
	while (sam_read1(bam_file, bam_header, bam_record) >= 0) {

		if (is_rna_bam_file)
			if ((bam_record->core.flag & (BAM_FSECONDARY | BAM_FUNMAP)) || (bam_record->core.flag & BAM_FPAIRED) && (bam_record->core.flag & BAM_FMUNMAP)) // ignore multi-mapping and unmapped reads
				continue;

		// fix contig number to match ours
		bam_record->core.tid = tid_to_contig[bam_record->core.tid];

		// when chimOutType==SeparateSAMold, then the chimeric alignments are in a separate file and we can load everything
		// when chimOutType==WithinBAM, then we need to extract them based on SAM attributes:
		if (separate_chimeric_bam_file && !is_rna_bam_file || // read everything from Chimeric.out.sam
		    !separate_chimeric_bam_file && is_rna_bam_file && // extract chimeric reads from Aligned.out.bam
		    ((bam_record->core.flag & BAM_FPAIRED) && !(bam_record->core.flag & BAM_FPROPER_PAIR) || // discordant mate
		     (bam_record->core.flag & BAM_FSUPPLEMENTARY))) { // supplementary
			add_chimeric_alignment(chimeric_alignments, bam_record, 0, 0, false, false, separate_chimeric_bam_file && (bam_record->core.flag & BAM_FSECONDARY) || (bam_record->core.flag & BAM_FSUPPLEMENTARY));
			continue; // supplementary alignments and discordant mates are added directly; only split reads and read-through fusions need to be buffered (see below)
		}

		// count mapped reads on interesting contigs
		if (interesting_tids[bam_record->core.tid])
			mapped_reads++;

		// for paired-end data we need to wait until we have read both mates
		bam1_t* previously_seen_mate = NULL;
		if (bam_record->core.flag & BAM_FPAIRED) {

			// try to insert the mate into the buffered BAM records
			// if there was already a record with the same read name, insertion will fail (->second set to false) and
			// previously_seen_mate->first will point to the mate which was already in the buffered BAM records
			pair<buffered_bam_records_t::iterator,bool> find_previously_seen_mate = buffered_bam_records.insert(pair<string,bam1_t*>((char*) bam_get_qname(bam_record), bam_record));
			if (!find_previously_seen_mate.second) { // this is the second mate we have seen
				previously_seen_mate = find_previously_seen_mate.first->second;
				buffered_bam_records.erase(find_previously_seen_mate.first); // remove from lookup buffer, we don't need it anymore
			}

		}

		if ((bam_record->core.flag & BAM_FPAIRED) && previously_seen_mate == NULL) { // this is the first mate with the given read name, which we encounter
			
			bam_record = bam_init1(); // allocate memory for the next record
			if (bam_record == NULL) {
				cerr << "ERROR: failed to allocate memory." << endl;
				exit(1);
			}


		} else { // single-end data or we have already read the first mate previously

			bool already_added = false;
			if (bam_aux_get(bam_record, "SA") != NULL || previously_seen_mate != NULL && bam_aux_get(previously_seen_mate, "SA") != NULL) { // split-read
				if (!separate_chimeric_bam_file) {
					add_chimeric_alignment(chimeric_alignments, bam_record);
					if (previously_seen_mate != NULL)
						add_chimeric_alignment(chimeric_alignments, previously_seen_mate);
				}
				already_added = true;
			}

			if (is_rna_bam_file && !already_added)
				already_added = extract_read_through_alignment(chimeric_alignments, bam_record, previously_seen_mate, gene_annotation_index, separate_chimeric_bam_file);

			if (!already_added && interesting_tids[bam_record->core.tid]) {
				coverage.add_alignment(bam_record);
				if (previously_seen_mate != NULL)
					coverage.add_alignment(previously_seen_mate);
			}

			if (previously_seen_mate != NULL)
				bam_destroy1(previously_seen_mate);
		}
	}

	// close BAM file
	bam_destroy1(bam_record);
	bam_hdr_destroy(bam_header);
	sam_close(bam_file);

	return chimeric_alignments.size();
}

void assign_strands_from_strandedness(chimeric_alignments_t& chimeric_alignments, const strandedness_t strandedness) {
	if (strandedness != STRANDEDNESS_NO) {
		for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {
			unsigned int first_in_pair = (chimeric_alignment->second[MATE1].first_in_pair) ? MATE1 : MATE2;
			unsigned int second_in_pair = (chimeric_alignment->second[MATE1].first_in_pair) ? MATE2 : MATE1;
			chimeric_alignment->second[first_in_pair].predicted_strand = complement_strand_if(chimeric_alignment->second[first_in_pair].strand, strandedness == STRANDEDNESS_REVERSE);
			chimeric_alignment->second[first_in_pair].predicted_strand_ambiguous = false;
			chimeric_alignment->second[second_in_pair].predicted_strand = complement_strand_if(chimeric_alignment->second[first_in_pair].predicted_strand, chimeric_alignment->second[first_in_pair].strand == chimeric_alignment->second[second_in_pair].strand);
			chimeric_alignment->second[second_in_pair].predicted_strand_ambiguous = false;
			if (chimeric_alignment->second.size() == 3) { // split read
				chimeric_alignment->second[SUPPLEMENTARY].predicted_strand = complement_strand_if(chimeric_alignment->second[SPLIT_READ].predicted_strand, chimeric_alignment->second[SUPPLEMENTARY].strand != chimeric_alignment->second[SPLIT_READ].strand);
				chimeric_alignment->second[SUPPLEMENTARY].predicted_strand_ambiguous = false;
			}
		}
	}
}

