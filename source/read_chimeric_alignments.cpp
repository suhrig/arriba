#include <algorithm>
#include <climits>
#include <iostream>
#include <string>
#include <vector>
#include "cram.h"
#include "htrie_map.h"
#include "sam.h"
#include "annotation.hpp"
#include "common.hpp"
#include "read_chimeric_alignments.hpp"
#include "read_stats.hpp"

using namespace std;

typedef tsl::htrie_map<char,bam1_t*> collated_bam_records_t;
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
		if (bam_cigar_op(current_op) == BAM_CREF_SKIP && // this is an intron
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

inline strand_t get_strand(const bam1_t* bam_record) {
	return (bam_record->core.flag & BAM_FREVERSE) ? REVERSE : FORWARD;
}

const unsigned char CLIP_NONE = 0;
const unsigned char CLIP_START = 1;
const unsigned char CLIP_END = 2;
void add_chimeric_alignment(mates_t& mates, const bam1_t* bam_record, const bool is_supplementary = false, const unsigned int cigar_op = 0, const unsigned char clip = CLIP_NONE) {

	// convert bam1_t structure into our own structure and discard information we don't need
	mates.single_end = !(bam_record->core.flag & BAM_FPAIRED);
	mates.duplicate = mates.duplicate || (bam_record->core.flag & BAM_FDUP);
	mates.resize(mates.size()+1);
	alignment_t& alignment = mates[mates.size()-1];
	alignment.strand = get_strand(bam_record);
	alignment.first_in_pair = bam_record->core.flag & BAM_FREAD1;
	alignment.contig = bam_record->core.tid;
	alignment.supplementary = is_supplementary;
	if (!is_supplementary) { // only keep sequence in memory, if this is not the supplementary alignment (because then it's already stored in the split-read)
		alignment.sequence.resize(bam_record->core.l_qseq);
		for (int i = 0; i < bam_record->core.l_qseq; ++i)
			alignment.sequence[i] = seq_nt16_str[bam_seqi(bam_get_seq(bam_record), i)];
	}

	// read-through alignments need to be split into a split-read and a supplementary alignment
	if (clip == CLIP_START) {
		alignment.start = bam_record->core.pos + bam_cigar2rlen(cigar_op, bam_get_cigar(bam_record));
		alignment.end = bam_endpos(bam_record) - 1;
		alignment.cigar.resize(bam_record->core.n_cigar - cigar_op + 1);
		alignment.cigar[0] = bam_cigar_gen(bam_cigar2qlen(cigar_op, bam_get_cigar(bam_record)), BAM_CSOFT_CLIP); // soft-clip start of read
		for (unsigned int i = cigar_op; i < bam_record->core.n_cigar; ++i) // copy cigar operations from <cigar_op> onwards
			alignment.cigar[i-cigar_op+1] = bam_get_cigar(bam_record)[i];
	} else if (clip == CLIP_END) {
		alignment.start = bam_record->core.pos;
		alignment.end = bam_record->core.pos + bam_cigar2rlen(cigar_op + 1, bam_get_cigar(bam_record)) - 1;
		alignment.cigar.resize(cigar_op + 2);
		for (unsigned int i = 0; i <= cigar_op; ++i) // copy cigar operations up until <cigar_op>
			alignment.cigar[i] = bam_get_cigar(bam_record)[i];
		alignment.cigar[cigar_op+1] = bam_cigar_gen(bam_record->core.l_qseq - bam_cigar2qlen(cigar_op+1, bam_get_cigar(bam_record)), BAM_CSOFT_CLIP); // soft-clip end of read
	} else { // do not clip, i.e., alignment is already split into a split-read and a supplementary alignment
		alignment.start = bam_record->core.pos;
		alignment.end = bam_endpos(bam_record) - 1;
		alignment.cigar.resize(bam_record->core.n_cigar);
		for (unsigned int i = 0; i < alignment.cigar.size(); ++i)
			alignment.cigar[i] = bam_get_cigar(bam_record)[i];
	}
}

bool extract_read_through_alignment(chimeric_alignments_t& chimeric_alignments, const string& read_name, bam1_t* forward_mate, bam1_t* reverse_mate, const gene_annotation_index_t& gene_annotation_index, const bool separate_chimeric_bam_file) {

	// find out which read is on the forward strand and which on the reverse
	if (get_strand(forward_mate) == REVERSE)
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
		position_t forward_read_pos, reverse_read_pos;
		bool forward_mate_has_intron = (forward_mate == NULL) ? false : find_spanning_intron(forward_mate, forward_gene_end, reverse_gene_start, forward_cigar_op, forward_read_pos);
		bool reverse_mate_has_intron = (reverse_mate == NULL) ? false : find_spanning_intron(reverse_mate, forward_gene_end, reverse_gene_start, reverse_cigar_op, reverse_read_pos);
		if (forward_mate_has_intron &&
		    (!reverse_mate_has_intron || forward_read_pos < reverse_mate->core.l_qseq - reverse_read_pos)) { // if both mates are clipped, use the one with the longer segment as anchor

			// store read-through alignments, unless they are already stored as chimeric alignments
			pair<chimeric_alignments_t::iterator,bool> mates = chimeric_alignments.insert(pair<string,mates_t>(read_name, mates_t()));
			if (mates.second) { // insertion succeeded => the alignments have not been stored previously

				// make split read and supplementary from forward mate
				add_chimeric_alignment(mates.first->second, forward_mate, false, forward_cigar_op+1, CLIP_START);
				add_chimeric_alignment(mates.first->second, forward_mate, true/*supplementary*/, forward_cigar_op-1, CLIP_END);

				if (reverse_mate != NULL) { // paired-end
					if (reverse_mate_has_intron) // reverse mate overlaps with breakpoint => clip it
						add_chimeric_alignment(mates.first->second, reverse_mate, false, reverse_cigar_op+1, CLIP_START);
					else // reverse mate overlaps with forward mate, but not with breakpoint => add it as is
						add_chimeric_alignment(mates.first->second, reverse_mate);
				}

				return true;
			}

		} else if (reverse_mate_has_intron) {

			// add read-through alignments, unless they are already stored as chimeric alignments
			pair<chimeric_alignments_t::iterator,bool> mates = chimeric_alignments.insert(pair<string,mates_t>(read_name, mates_t()));
			if (mates.second) { // insertion succeeded => the alignments have not been stored previously

				// make split read and supplementary from reverse mate
				add_chimeric_alignment(mates.first->second, reverse_mate, true/*supplementary*/, reverse_cigar_op+1, CLIP_START);
				add_chimeric_alignment(mates.first->second, reverse_mate, false, reverse_cigar_op-1, CLIP_END);
	
				if (forward_mate != NULL) { // paired-end
					if (forward_mate_has_intron) // forward mate overlaps with breakpoints => clip it at the end
						add_chimeric_alignment(mates.first->second, forward_mate, false, forward_cigar_op-1, CLIP_END);
					else // forward mate overlaps with reverse mate, but not with breakpoint => add it as is
						add_chimeric_alignment(mates.first->second, forward_mate);
				}

				return true;
			}

		// check for possibility (2) => the mates must be contained within different genes
		} else if (forward_mate != NULL && reverse_mate != NULL && // paired-end
		           reverse_mate->core.pos   >= reverse_gene_start &&
		           bam_endpos(forward_mate) <= forward_gene_end) {

			// add read-through alignments, unless they are already stored as chimeric alignments
			pair<chimeric_alignments_t::iterator,bool> mates = chimeric_alignments.insert(pair<string,mates_t>(read_name, mates_t()));
			if (mates.second) { // insertion succeeded => the alignments have not been stored previously
				add_chimeric_alignment(mates.first->second, forward_mate);
				add_chimeric_alignment(mates.first->second, reverse_mate);
			}
			return true;

		} // else possibility (3)
	}

	return false;
}

// when the insert size is small, two mates fully overlap and the ends are clipped, because they contain adapter sequence
// => ignore them, do not consider them for alignment of candidate internal tandem repeat reads
bool clipped_sequence_is_adapter(const bam1_t* mate1, const bam1_t* mate2) {
	if (mate1 == NULL || mate2 == NULL)
		return false;
	if (mate1->core.pos == mate2->core.pos) {
		if (get_strand(mate1) == REVERSE && bam_cigar_op(bam_get_cigar(mate1)[0]) == BAM_CSOFT_CLIP &&
		    get_strand(mate2) == FORWARD && bam_cigar_op(bam_get_cigar(mate2)[mate2->core.n_cigar-1]) == BAM_CSOFT_CLIP &&
		    bam_cigar_oplen(bam_get_cigar(mate1)[0]) == bam_cigar_oplen(bam_get_cigar(mate2)[mate2->core.n_cigar-1]))
			return true;
		if (get_strand(mate2) == REVERSE && bam_cigar_op(bam_get_cigar(mate2)[0]) == BAM_CSOFT_CLIP &&
		    get_strand(mate1) == FORWARD && bam_cigar_op(bam_get_cigar(mate1)[mate1->core.n_cigar-1]) == BAM_CSOFT_CLIP &&
		    bam_cigar_oplen(bam_get_cigar(mate2)[0]) == bam_cigar_oplen(bam_get_cigar(mate1)[mate1->core.n_cigar-1]))
			return true;
	}
	return false;
}

// STAR is bad at aligning internal tandem duplications
// => if we see a clipped read, check manually if it can be aligned as a tandem duplication
bool is_tandem_duplication(const bam1_t* bam_record, const assembly_t& assembly, const unsigned int max_itd_length, alignment_t& tandem_alignment) {

	const unsigned int min_clipped_length = 12; // ignore reads with less than this many clipped bases
	const unsigned int min_duplication_length = 9; // lowering this value is probably pointless, because STAR aligns as an indel instead of clipped read
	const unsigned int max_duplication_length = max_itd_length;
	const unsigned int max_mismatches = 1; // abort alignment after hitting this many mismatches
	const unsigned int max_non_template_bases = 6; // ignore mismatches at the beginning of the tandem alignment, because this is a frequent occurrence with ITDs
	const unsigned int min_alignment_length = 15; // at least this many bases must align to consider the alignment valid

	if (bam_record == NULL)
		return false;

	// check if read is clipped at start or end, and if so, find the longer clipped segment
	unsigned int clipped_sequence_length = 0;
	unsigned int clipped_sequence_position = 0;
	bool clipped_start = true;
	int alignment_direction = +1;
	int alignment_window_start = 0;
	int alignment_window_end = 0;
	if (bam_cigar_op(bam_get_cigar(bam_record)[0]) == BAM_CSOFT_CLIP && bam_cigar_oplen(bam_get_cigar(bam_record)[0]) >= min_clipped_length) {
		// read is clipped at start
		clipped_sequence_length = bam_cigar_oplen(bam_get_cigar(bam_record)[0]);
		clipped_sequence_position = 0;
		alignment_direction = -1;
		alignment_window_start = bam_record->core.pos + min_duplication_length - clipped_sequence_length;
		alignment_window_end = bam_record->core.pos + max_duplication_length - clipped_sequence_length;
		clipped_start = true;
	}
	if (bam_cigar_op(bam_get_cigar(bam_record)[bam_record->core.n_cigar-1]) == BAM_CSOFT_CLIP && bam_cigar_oplen(bam_get_cigar(bam_record)[bam_record->core.n_cigar-1]) >= max(min_clipped_length, clipped_sequence_length)) {
		// read is clipped at end
		clipped_sequence_length = bam_cigar_oplen(bam_get_cigar(bam_record)[bam_record->core.n_cigar-1]);
		clipped_sequence_position = bam_record->core.l_qseq - clipped_sequence_length;
		alignment_direction = +1;
		alignment_window_start = bam_endpos(bam_record) - max_duplication_length;
		alignment_window_end = bam_endpos(bam_record) - min_duplication_length;
		clipped_start = false;
	}
	if (clipped_sequence_length == 0)
		return false; // read is not clipped

	if (assembly.find(bam_record->core.tid) == assembly.end())
		return false; // contig sequence unavailable and thus no way to make an alignment
	const string& contig_sequence = assembly.at(bam_record->core.tid);
	if (alignment_window_end + max_duplication_length + clipped_sequence_length + 1 >= contig_sequence.size() ||
	    alignment_window_start <= (int) (max_duplication_length + clipped_sequence_length + 1))
		return false; // ignore alignments close to contig boundaries to avoid array out-of-bounds errors

	// try to align clipped sequence in a window of size <max_duplication_length>
	string clipped_sequence;
	clipped_sequence.resize(clipped_sequence_length);
	for (unsigned int i = 0; i < clipped_sequence_length; ++i)
		clipped_sequence[i] = seq_nt16_str[bam_seqi(bam_get_seq(bam_record), clipped_sequence_position + i)];
	for (int contig_pos = alignment_window_start; contig_pos <= alignment_window_end; ++contig_pos) {

		// align at given position and abort when too many mismatches have been encountered
		unsigned int matches = 0;
		unsigned int mismatches = 0;
		tandem_alignment.start = contig_sequence.size();
		tandem_alignment.end = -1;
		for (unsigned int i = 0; i < clipped_sequence_length && mismatches <= max_mismatches; i++) {
			int read_pos = (alignment_direction == +1) ? i : clipped_sequence_length - 1 - i;
			if (contig_sequence[contig_pos+read_pos] == clipped_sequence[read_pos]) {
				matches++;
				if (contig_pos + read_pos < tandem_alignment.start)
					tandem_alignment.start = contig_pos + read_pos;
				if (contig_pos + read_pos > tandem_alignment.end)
					tandem_alignment.end = contig_pos + read_pos;
			} else if (i >= max_non_template_bases) {
				mismatches++;
			}
		}

		// return tandem alignment as result if it has sufficient quality
		if (matches >= min_alignment_length || matches + mismatches == clipped_sequence_length) {
			tandem_alignment.strand = get_strand(bam_record);
			tandem_alignment.first_in_pair = bam_record->core.flag & BAM_FREAD1;
			tandem_alignment.contig = bam_record->core.tid;
			tandem_alignment.supplementary = !(bam_record->core.flag & BAM_FPAIRED) ||
			                                 clipped_start  && get_strand(bam_record) == FORWARD ||
			                                 !clipped_start && get_strand(bam_record) == REVERSE;
			if (!tandem_alignment.supplementary) { // only keep sequence in memory, if this is not the supplementary alignment (because then it's already stored in the split-read)
				tandem_alignment.sequence.resize(bam_record->core.l_qseq);
				for (int i = 0; i < bam_record->core.l_qseq; ++i)
					tandem_alignment.sequence[i] = seq_nt16_str[bam_seqi(bam_get_seq(bam_record), i)];
			}
			// construct CIGAR string
			uint32_t clip_left = (clipped_start) ? 0 : bam_record->core.l_qseq - clipped_sequence_length;
			uint32_t clip_right = (clipped_start) ? bam_record->core.l_qseq - clipped_sequence_length : 0;
			if (tandem_alignment.start > contig_pos)
				clip_left += tandem_alignment.start - contig_pos;
			if (tandem_alignment.end < contig_pos + (int) clipped_sequence_length - 1)
				clip_right += contig_pos + clipped_sequence_length - 1 - tandem_alignment.end;
			if (clip_left > 0)
				tandem_alignment.cigar.push_back(bam_cigar_gen(clip_left, BAM_CSOFT_CLIP));
			tandem_alignment.cigar.push_back(bam_cigar_gen(tandem_alignment.end - tandem_alignment.start + 1, BAM_CMATCH));
			if (clip_right > 0)
				tandem_alignment.cigar.push_back(bam_cigar_gen(clip_right, BAM_CSOFT_CLIP));
			return true;
		}
	}

	return false; // if we get here, no tandem duplication alignment was found
}

// remove alignments when supplementary flags are missing or when there are too many/few alignment records
// in addition, reformat single-end reads as if they were paired-end, such that the rest of Arriba does not need to care about single-end vs. paired-end
unsigned int remove_malformed_alignments(chimeric_alignments_t& chimeric_alignments) {

	unsigned int malformed_count = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end();) {

		bool malformed = false;

		if (chimeric_alignment->second.single_end) {
			if (chimeric_alignment->second.size() == 2 &&
			    (chimeric_alignment->second[MATE1].supplementary !=/*xor*/ chimeric_alignment->second[MATE2].supplementary)) { // there must be exactly one supplementary flag

				// use the alignment with the shorter anchor as the SUPPLEMENTARY and the longer one as the SPLIT_READ
				// and copy the split read in the MATE1 place (to simulate paired-end data)
				if (chimeric_alignment->second[MATE1].end - chimeric_alignment->second[MATE1].start > chimeric_alignment->second[MATE2].end - chimeric_alignment->second[MATE2].start) {
					chimeric_alignment->second.push_back(chimeric_alignment->second[MATE2]);
					chimeric_alignment->second[MATE2] = chimeric_alignment->second[MATE1];
				} else {
					chimeric_alignment->second.push_back(chimeric_alignment->second[MATE1]);
					chimeric_alignment->second[MATE1] = chimeric_alignment->second[MATE2];
				}

				// MATE1 and SPLIT_READ must have the sequence, SUPPLEMENTARY must not
				if (!chimeric_alignment->second[MATE1].supplementary) {
					chimeric_alignment->second[SPLIT_READ].sequence = chimeric_alignment->second[MATE1].sequence;
				} else if (!chimeric_alignment->second[SPLIT_READ].supplementary) {
					chimeric_alignment->second[MATE1].sequence = chimeric_alignment->second[SPLIT_READ].sequence;
				} else { // !chimeric_alignment->second[SUPPLEMENTARY].supplementary
					chimeric_alignment->second[MATE1].sequence = chimeric_alignment->second[SUPPLEMENTARY].sequence;
					chimeric_alignment->second[SPLIT_READ].sequence = chimeric_alignment->second[SUPPLEMENTARY].sequence;
				}
				chimeric_alignment->second[SUPPLEMENTARY].sequence.clear();

				// set supplementary flag like it would be set if we had paired-end data
				chimeric_alignment->second[SUPPLEMENTARY].supplementary = true;
				chimeric_alignment->second[MATE1].supplementary = false;
				chimeric_alignment->second[SPLIT_READ].supplementary = false;

				// set strands like they would be set if we had paired-end data
				if (chimeric_alignment->second[SPLIT_READ].sequence.length() - chimeric_alignment->second[SPLIT_READ].preclipping() - ((chimeric_alignment->second[SPLIT_READ].strand == chimeric_alignment->second[SUPPLEMENTARY].strand) ? chimeric_alignment->second[SUPPLEMENTARY].postclipping() : chimeric_alignment->second[SUPPLEMENTARY].preclipping()) <
				    chimeric_alignment->second[SPLIT_READ].sequence.length() - chimeric_alignment->second[SPLIT_READ].postclipping() - ((chimeric_alignment->second[SPLIT_READ].strand == chimeric_alignment->second[SUPPLEMENTARY].strand) ? chimeric_alignment->second[SUPPLEMENTARY].preclipping() : chimeric_alignment->second[SUPPLEMENTARY].postclipping())) {
					if (chimeric_alignment->second[SPLIT_READ].strand == FORWARD) {
						chimeric_alignment->second[MATE1].strand = complement_strand(chimeric_alignment->second[MATE1].strand);
					} else {
						chimeric_alignment->second[SPLIT_READ].strand = complement_strand(chimeric_alignment->second[SPLIT_READ].strand);
						chimeric_alignment->second[SUPPLEMENTARY].strand = complement_strand(chimeric_alignment->second[SUPPLEMENTARY].strand);
					}
				} else {
					if (chimeric_alignment->second[SPLIT_READ].strand == REVERSE) {
						chimeric_alignment->second[MATE1].strand = complement_strand(chimeric_alignment->second[MATE1].strand);
					} else {
						chimeric_alignment->second[SPLIT_READ].strand = complement_strand(chimeric_alignment->second[SPLIT_READ].strand);
						chimeric_alignment->second[SUPPLEMENTARY].strand = complement_strand(chimeric_alignment->second[SUPPLEMENTARY].strand);
					}
				}

				// the split read must be clipped at the end, the supplementary alignment at the beginning
				// and the clipped segment of the split read must encompass the supplementary and vice versa
				unsigned int clipped_bases_split_read = (chimeric_alignment->second[SPLIT_READ].strand == FORWARD) ? chimeric_alignment->second[SPLIT_READ].preclipping() : chimeric_alignment->second[SPLIT_READ].postclipping();
				unsigned int clipped_bases_supplementary = (chimeric_alignment->second[SUPPLEMENTARY].strand == FORWARD) ? chimeric_alignment->second[SUPPLEMENTARY].postclipping() : chimeric_alignment->second[SUPPLEMENTARY].preclipping();
				if (clipped_bases_split_read < chimeric_alignment->second[SPLIT_READ].sequence.size() - clipped_bases_supplementary ||
				    clipped_bases_supplementary < chimeric_alignment->second[SPLIT_READ].sequence.size() - clipped_bases_split_read)
					malformed = true;

			} else {
				// if we get here, there are either too many alignments with the same name or too few
				// or something is wrong with the supplementary flags
				malformed = true;
			}

		} else { // paired-end

			if (chimeric_alignment->second.size() == 3) { // split read

				// make sure supplementary alignment is in the SUPPLEMENTARY place
				if (chimeric_alignment->second[MATE1].supplementary) {
					swap(chimeric_alignment->second[MATE1], chimeric_alignment->second[SUPPLEMENTARY]);
				} else if (chimeric_alignment->second[MATE2].supplementary) {
					swap(chimeric_alignment->second[MATE2], chimeric_alignment->second[SUPPLEMENTARY]);
				}

				// make sure the split read is in the SPLIT_READ place
				if (chimeric_alignment->second[SPLIT_READ].first_in_pair != chimeric_alignment->second[SUPPLEMENTARY].first_in_pair)
					swap(chimeric_alignment->second[MATE1], chimeric_alignment->second[MATE2]);

				// make sure we have exactly one supplementary alignment, or else something is wrong
				if (chimeric_alignment->second[MATE1].supplementary ||
				    chimeric_alignment->second[SPLIT_READ].supplementary ||
				    !chimeric_alignment->second[SUPPLEMENTARY].supplementary)
					malformed = true;

				// mate1 and mate2 should align in a colinear fashion
				if (chimeric_alignment->second[MATE1].contig != chimeric_alignment->second[SPLIT_READ].contig ||
				    chimeric_alignment->second[MATE1].strand == chimeric_alignment->second[SPLIT_READ].strand)
					malformed = true;

				// the split read must be clipped at the end, the supplementary alignment at the beginning
				// and the clipped segment of the split read must encompass the supplementary and vice versa
				unsigned int clipped_bases_split_read = (chimeric_alignment->second[SPLIT_READ].strand == FORWARD) ? chimeric_alignment->second[SPLIT_READ].preclipping() : chimeric_alignment->second[SPLIT_READ].postclipping();
				unsigned int clipped_bases_supplementary = (chimeric_alignment->second[SUPPLEMENTARY].strand == FORWARD) ? chimeric_alignment->second[SUPPLEMENTARY].postclipping() : chimeric_alignment->second[SUPPLEMENTARY].preclipping();
				if (clipped_bases_split_read < chimeric_alignment->second[SPLIT_READ].sequence.size() - clipped_bases_supplementary ||
				    clipped_bases_supplementary < chimeric_alignment->second[SPLIT_READ].sequence.size() - clipped_bases_split_read)
					malformed = true;

			} else if (chimeric_alignment->second.size() == 2) { // discordant mate
				// none of the mates should have the supplementary bit set, or else something is wrong
				if (chimeric_alignment->second[MATE1].supplementary || chimeric_alignment->second[MATE2].supplementary)
					malformed = true;
			} else {
				// if we get here, there are either too many alignments with the same name or too few
				malformed = true;
			}

		}

		if (malformed) {
			malformed_count++;
			chimeric_alignment = chimeric_alignments.erase(chimeric_alignment);
		} else {
			++chimeric_alignment;
		}
	}

	return malformed_count;
}

// paired-end overlap alignment produces chimeric alignments where the wrong end is the clipped segment,
// when reads are longer than the insert size
// we must ignore those, because they are probably adapters
bool is_clipped_at_correct_end(const bam1_t* bam_record) {
	if (!(bam_record->core.flag & BAM_FPAIRED))
		return true; // single-end reads may be clipped at either end
	//paired-end read must be clipped at the fragment end
	unsigned int clipped_end;
	if (bam_record->core.flag & BAM_FSUPPLEMENTARY)
		clipped_end = (get_strand(bam_record) == FORWARD) ? bam_record->core.n_cigar-1 : 0;
	else
		clipped_end = (get_strand(bam_record) == FORWARD) ? 0 : bam_record->core.n_cigar-1;
	uint32_t clipped_cigar = bam_cigar_op(bam_get_cigar(bam_record)[clipped_end]);
	return clipped_cigar == BAM_CSOFT_CLIP || clipped_cigar == BAM_CHARD_CLIP;
}

unsigned int read_chimeric_alignments(const string& bam_file_path, const assembly_t& assembly, const string& assembly_file_path, chimeric_alignments_t& chimeric_alignments, unsigned long int& mapped_reads, vector<unsigned long int>& mapped_viral_reads_by_contig, coverage_t& coverage, contigs_t& contigs, vector<string>& original_contig_names, const string& interesting_contigs, const string& viral_contigs, const gene_annotation_index_t& gene_annotation_index, const bool separate_chimeric_bam_file, const bool is_rna_bam_file, const bool external_duplicate_marking, const unsigned int max_itd_length) {

	// open BAM file
	samFile* bam_file = sam_open(bam_file_path.c_str(), "rb");
	crash(bam_file == NULL, "failed to open SAM file");
	if (bam_file->is_cram)
		cram_set_option(bam_file->fp.cram, CRAM_OPT_REFERENCE, assembly_file_path.c_str());
	bam_hdr_t* bam_header = sam_hdr_read(bam_file);
	crash(bam_header == NULL, "failed to read SAM header");

	// add contigs which are not yet listed in <contigs>
	// and make a map tid -> contig, because the contig IDs in the BAM file need not necessarily match the contig IDs in the GTF file
	tid_to_contig_t tid_to_contig(bam_header->n_targets);
	vector<bool> interesting_tids(bam_header->n_targets);
	for (int target = 0; target < bam_header->n_targets; ++target) {
		string contig_name = removeChr(bam_header->target_name[target]);
		contigs.insert(pair<string,contig_t>(contig_name, contigs.size())); // this fails (i.e., nothing is inserted), if the contig already exists
		crash(contigs.size() == USHRT_MAX - 1, "too many contigs");
		if (contigs.size() > original_contig_names.size())
			original_contig_names.resize(contigs.size());
		original_contig_names[contigs[contig_name]] = bam_header->target_name[target];
		tid_to_contig[target] = contigs[contig_name];
		if (is_rna_bam_file) {
			if (tid_to_contig[target] >= (int) interesting_tids.size())
				interesting_tids.resize(tid_to_contig[target]+1);
			interesting_tids[tid_to_contig[target]] = is_interesting_contig(contig_name, interesting_contigs);
		}
	}
	coverage.resize(contigs, assembly);

	// make sure we have the sequence of all interesting contigs, otherwise later steps will crash
	for (contigs_t::iterator contig = contigs.begin(); contig != contigs.end(); ++contig)
		crash(assembly.find(contig->second) == assembly.end() && is_interesting_contig(contig->first, interesting_contigs), "could not find sequence of contig '" + contig->first + "'");

	// convert viral contigs to vector of booleans for faster lookup
	vector<bool> viral_contigs_bool(contigs.size());
	for (contigs_t::iterator contig = contigs.begin(); contig != contigs.end(); ++contig)
		viral_contigs_bool[contig->second] = is_interesting_contig(contig->first, viral_contigs);
	mapped_viral_reads_by_contig.resize(contigs.size());

	// read BAM records
	bam1_t* bam_record = bam_init1();
	crash(bam_record == NULL, "failed to allocate memory.");
	collated_bam_records_t collated_bam_records; // holds the first mate until we have found the second
	bool no_chimeric_reads = true;
	unsigned int missing_hi_tag = 0;
	unsigned int malformed_count = 0;
	string read_name;
	int sam_read1_status;
	while ((sam_read1_status = sam_read1(bam_file, bam_header, bam_record)) >= 0) {

		if (is_rna_bam_file)
			if ((bam_record->core.flag & BAM_FUNMAP) || (bam_record->core.flag & BAM_FPAIRED) && (bam_record->core.flag & BAM_FMUNMAP))
				continue; // ignore unmapped reads

		int64_t hit_index = 1;
		if (!separate_chimeric_bam_file) { // ignore HI tag in Chimeric.out.sam, because it only contains unique hits anyway
			uint8_t* hi_tag = bam_aux_get(bam_record, "HI");
			if (hi_tag != NULL) {
				hit_index = bam_aux2i(hi_tag);
			} else if (bam_record->core.flag & BAM_FSECONDARY) {
				missing_hi_tag++;
				continue; // ignore secondary alignments when HI tag is missing, because multi-mapping alignments could not be segregated
			}
		}
		read_name = (char*) bam_get_qname(bam_record);
		read_name += "," + to_string(hit_index); // append HI tag to name to enable segregation of multi-mapping reads

		// fix contig number to match ours
		bam_record->core.tid = tid_to_contig[bam_record->core.tid];

		// add supplementary alignments directly to the chimeric alignments without collating
		if (separate_chimeric_bam_file && !is_rna_bam_file && (bam_record->core.flag & BAM_FSECONDARY)) { // extract supplementary reads from Chimeric.out.sam
			add_chimeric_alignment(chimeric_alignments[read_name], bam_record, true/*supplementary*/);
			no_chimeric_reads = false;
			continue;
		}

		// add supplementary alignments directly to the chimeric alignments without collating
		if (is_rna_bam_file && (bam_record->core.flag & BAM_FSUPPLEMENTARY)) { // extract supplementary reads from Aligned.out.bam
			if (!separate_chimeric_bam_file) { // don't load alignments twice (from Chimeric.out.sam and from Aligned.out.bam)
				if (is_clipped_at_correct_end(bam_record))
					add_chimeric_alignment(chimeric_alignments[read_name], bam_record, true/*supplementary*/);
				else
					malformed_count++;
				no_chimeric_reads = false;
			}
			continue;
		}

		// count mapped reads on interesting contigs
		if (interesting_tids[bam_record->core.tid])
			mapped_reads++;

		// add discordant mates directly to the chimeric alignments without collating
		if (is_rna_bam_file && (bam_record->core.flag & BAM_FPAIRED) && !(bam_record->core.flag & BAM_FPROPER_PAIR)) { // extract discordant mates from Aligned.out.bam
			if (!separate_chimeric_bam_file) { // don't load alignments twice (from Chimeric.out.sam and from Aligned.out.bam)
				add_chimeric_alignment(chimeric_alignments[read_name], bam_record);
				no_chimeric_reads = false;
			}
			// compute coverage of discordant mates individually as if they were single-end reads
			if (!external_duplicate_marking || !(bam_record->core.flag & BAM_FDUP)) {
				bam_record->core.flag &= !BAM_FPAIRED;
				coverage.add_fragment(bam_record, NULL, true);
			}
			continue;
		}

		// for paired-end data we need to wait until we have read both mates
		bam1_t* previously_seen_mate = NULL;
		if (bam_record->core.flag & BAM_FPAIRED) {

			// try to insert the mate into the collated BAM records
			// if there was already a record with the same read name, insertion will fail (->second set to false) and
			// previously_seen_mate->first will point to the mate which was already in the collated BAM records
			pair<collated_bam_records_t::iterator,bool> find_previously_seen_mate = collated_bam_records.insert(read_name.c_str(), bam_record);
			if (!find_previously_seen_mate.second) { // this is the second mate we have seen
				previously_seen_mate = *find_previously_seen_mate.first;
				collated_bam_records.erase(find_previously_seen_mate.first);
			}

		}

		if ((bam_record->core.flag & BAM_FPAIRED) && previously_seen_mate == NULL) { // this is the first mate with the given read name, which we encounter
			
			bam_record = bam_init1(); // allocate memory for the next record
			crash(bam_record == NULL, "failed to allocate memory");

		} else { // single-end data or we have already read the first mate previously

			if (separate_chimeric_bam_file && !is_rna_bam_file) { // this is Chimeric.out.sam => load everything

				mates_t& mates = chimeric_alignments[read_name];
				add_chimeric_alignment(mates, bam_record);
				if (previously_seen_mate != NULL)
					add_chimeric_alignment(mates, previously_seen_mate);
				no_chimeric_reads = false;

			} else { // this is Aligned.out.bam => load only discordant mates and split reads, and only when there is no Chimeric.out.sam

				bool is_read_through_alignment = false;
				alignment_t tandem_alignment;

				if (bam_aux_get(bam_record, "SA") != NULL && is_clipped_at_correct_end(bam_record) ||
				    previously_seen_mate != NULL && bam_aux_get(previously_seen_mate, "SA") != NULL && is_clipped_at_correct_end(previously_seen_mate)) { // split-read with SA tag
					if (!separate_chimeric_bam_file) {
						mates_t& mates = chimeric_alignments[read_name];
						add_chimeric_alignment(mates, bam_record);
						if (previously_seen_mate != NULL)
							add_chimeric_alignment(mates, previously_seen_mate);
						no_chimeric_reads = false;
					}
				} else if (!clipped_sequence_is_adapter(bam_record, previously_seen_mate) &&
				           (previously_seen_mate == NULL || get_strand(bam_record) != get_strand(previously_seen_mate)) && // strands must be different, so we can distinguish mate1 from mate2
				           (is_tandem_duplication(bam_record, assembly, max_itd_length, tandem_alignment) || // is it a tandem duplication that STAR failed to align?
				            is_tandem_duplication(previously_seen_mate, assembly, max_itd_length, tandem_alignment))) {
					if (!separate_chimeric_bam_file || is_rna_bam_file && chimeric_alignments.find(read_name) == chimeric_alignments.end()) {
						mates_t& mates = chimeric_alignments[read_name];
						add_chimeric_alignment(mates, bam_record, get_strand(bam_record) == tandem_alignment.strand && !tandem_alignment.supplementary);
						if (previously_seen_mate != NULL)
							add_chimeric_alignment(mates, previously_seen_mate, get_strand(previously_seen_mate) == tandem_alignment.strand && !tandem_alignment.supplementary);
						mates.push_back(tandem_alignment);
					}
				} else { // could be a read-through alignment
					is_read_through_alignment = extract_read_through_alignment(chimeric_alignments, read_name, bam_record, previously_seen_mate, gene_annotation_index, separate_chimeric_bam_file);

					// count mapped reads on viral contigs to detect viral infection
					if (viral_contigs_bool[bam_record->core.tid]) {
						// only count perfectly matching alignments to ignore alignment artifacts
						for (bam1_t* mate = bam_record; mate != NULL; mate = (mate == previously_seen_mate) ? NULL : previously_seen_mate) {
							bool pristine_alignment = true;
							for (unsigned int i = 0; i < mate->core.n_cigar && pristine_alignment; i++) {
								uint32_t cigar_op = bam_cigar_op(bam_get_cigar(mate)[i]);
								if (cigar_op != BAM_CREF_SKIP && cigar_op != BAM_CMATCH && cigar_op != BAM_CDIFF)
									pristine_alignment = false;
							}
							if (pristine_alignment)
								mapped_viral_reads_by_contig[mate->core.tid]++;
						}
					}
				}

				if (!external_duplicate_marking || !(bam_record->core.flag & BAM_FDUP))
					coverage.add_fragment(bam_record, previously_seen_mate, is_read_through_alignment);
			}

			if (previously_seen_mate != NULL)
				bam_destroy1(previously_seen_mate);
		}
	}

	crash(sam_read1_status < -1, "failed to load alignments");

	// close BAM file
	bam_destroy1(bam_record);
	bam_hdr_destroy(bam_header);
	sam_close(bam_file);

	// sanity check: input files should not be empty
	crash(is_rna_bam_file && mapped_reads == 0, "no normal reads found");
	// sanity check: remove malformed alignments
	malformed_count += remove_malformed_alignments(chimeric_alignments);
	if (malformed_count > 0)
		cerr << "WARNING: " << malformed_count << " SAM records were malformed and ignored" << endl;
	// sanity check: there should be at least 1 chimeric read, or else Arriba is probably not being used properly
	if (separate_chimeric_bam_file && !is_rna_bam_file || // this is Chimeric.out.sam
	    !separate_chimeric_bam_file) // this is Aligned.out.bam and STAR was run with --chimOutType WithinBAM
		crash(no_chimeric_reads, "no split reads or discordant mates found (STAR must either be run with '--chimOutType WithinBAM' or the file 'Chimeric.out.sam' must be passed to Arriba via the argument -c)");
	// sanity check: multi-mapping chimeric reads should have the HI tag
	if (missing_hi_tag > 0)
		cerr << "WARNING: " << missing_hi_tag << " secondary alignments lack the 'HI' tag and were ignored (STAR must be run with '--outSAMattributes HI' for Arriba to make use of multi-mapping reads for fusion detection)" << endl;

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

unsigned int mark_multimappers(chimeric_alignments_t& chimeric_alignments) {
	unsigned int count = 0;
	if (!chimeric_alignments.empty())
		for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); next(chimeric_alignment) != chimeric_alignments.end(); ++chimeric_alignment)
			if (strip_hi_tag_from_read_name(chimeric_alignment->first) == strip_hi_tag_from_read_name(next(chimeric_alignment)->first)) {
				chimeric_alignment->second.multimapper = true;
				next(chimeric_alignment)->second.multimapper = true;
				count++;
			}
	return count;
}

