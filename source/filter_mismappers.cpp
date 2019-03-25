#include <cmath>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
#include "common.hpp"
#include "annotation.hpp"
#include "assembly.hpp"
#include "filter_mismappers.hpp"

using namespace std;

typedef set<position_t> splice_sites_t;
typedef unordered_map<gene_t,splice_sites_t> splice_sites_by_gene_t;

void get_downstream_splice_sites(const gene_t gene, const exon_annotation_index_t& exon_annotation_index, splice_sites_t& splice_sites) {

	// nothing to do, if there are no exons on the given contig
	if (exon_annotation_index[gene->contig].empty()) {
		splice_sites.clear();
		return;
	}

	// find all downstream-oriented splice-sites of the given gene (because alignment is oriented downstream)
	exon_contig_annotation_index_t::const_iterator exons = exon_annotation_index[gene->contig].lower_bound(gene->start);
	while (exons != exon_annotation_index[gene->contig].end() && exons->first <= gene->end) {
		if (is_breakpoint_spliced(gene, DOWNSTREAM, exons->first, exon_annotation_index))
			splice_sites.insert(exons->first);
		++exons;
	}
}

kmer_as_int_t kmer_to_int(const string& kmer, const string::size_type position, const char kmer_length) {
	kmer_as_int_t result = 0;
	for (char base = 0; base < kmer_length; ++base) {
		result = result<<2;
		switch (kmer.c_str()[position + base]) {
			case 'T': result += 0; break;
			case 'G': result += 1; break;
			case 'C': result += 2; break;
			default:  result += 3; break;
		}
	}
	return result;
}

void make_kmer_index(const fusions_t& fusions, const assembly_t& assembly, const char kmer_length, kmer_indices_t& kmer_indices) {

	// find genes which are involved in fusions which have not been discarded yet
	gene_set_t genes_to_filter;
	for (fusions_t::const_iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {
		if (fusion->second.filter != NULL)
			continue;
		if (fusion->second.gene1 == fusion->second.gene2)
			continue; // comparing sequence similarity only makes sense between different genes
		genes_to_filter.insert(fusion->second.gene1);
		genes_to_filter.insert(fusion->second.gene2);
	}

	// store positions of kmers in hash
	for (gene_set_t::iterator gene = genes_to_filter.begin(); gene != genes_to_filter.end(); ++gene) {
		const string& contig_sequence = assembly.at((**gene).contig);
		if ((int) kmer_indices.size() <= (**gene).contig)
			kmer_indices.resize((**gene).contig+1);
		for (position_t pos = (**gene).start; pos + kmer_length < (**gene).end; pos++)
			if (contig_sequence[pos] != 'N') // don't index masked regions, as long stretches of N's inflate the number of hits
				kmer_indices[(**gene).contig][kmer_to_int(contig_sequence, pos, kmer_length)].push_back(pos);
	}

	// sort kmer hits by increasing position, so that we can go through the list sequentially
	for (kmer_indices_t::iterator kmer_index = kmer_indices.begin(); kmer_index != kmer_indices.end(); ++kmer_index)
		for (kmer_index_t::iterator kmer_hits = kmer_index->begin(); kmer_hits != kmer_index->end(); ++kmer_hits) {
			sort(kmer_hits->second.begin(), kmer_hits->second.end());
			// when genes overlap the same kmer hit might be added multiple times => only keep unique entries
			auto last = unique(kmer_hits->second.begin(), kmer_hits->second.end());
			kmer_hits->second.erase(last, kmer_hits->second.end());
		}
}

bool align(int score, const string& read_sequence, int read_pos, const string& contig_sequence, const int gene_pos, const position_t gene_start, const position_t gene_end, const kmer_index_t& kmer_index, const char kmer_length, const splice_sites_t& splice_sites, const int min_score, int max_deletions) {

	int skipped_bases = 0;

	for (/* read_pos comes from parameters */;
	     read_pos + kmer_length < (int) read_sequence.length() && // don't run over end of read
	     read_pos + min_score <= (int) read_sequence.length() + score + 2*kmer_length; // give up, when we can impossibly get above min_score, because we are near the end of the read
	                                                                             // 2*kmer_length takes into account that the score can improve, if we can extend to the left (up to kmer_length)
	     read_pos++, score--, skipped_bases++) { // if a base cannot be aligned, go to the next, but give -1 penalty and increase the number of skipped bases

		auto kmer_hits = kmer_index.find(kmer_to_int(read_sequence, read_pos, kmer_length));
		if (kmer_hits == kmer_index.end())
			continue; // kmer not found on given contig

		for (auto kmer_hit = lower_bound(kmer_hits->second.begin(), kmer_hits->second.end(), gene_pos); kmer_hit != kmer_hits->second.end() && *kmer_hit < gene_end; ++kmer_hit) {

			int extended_score = score + kmer_length;
			if (read_pos == skipped_bases) // so far, all bases at the beginning of the read have been skipped
				extended_score += skipped_bases; // this effectively removes any penalties on leading mismatches (as in local alignment)
			if (extended_score >= min_score)
				return true;

			// extend match locally to the left
			{
				int extended_read_pos = read_pos - 1;
				int extended_gene_pos = *kmer_hit - 1;
				unsigned int mismatch_count = 0;
				while (extended_read_pos >= read_pos - skipped_bases && // only align yet unaligned bases
				       extended_gene_pos >= gene_start) {

					if (read_sequence[extended_read_pos] == contig_sequence[extended_gene_pos]) {

						// read and gene sequences match => increase score and go the next base
						if (read_pos == skipped_bases)
							extended_score += 1; // increase by 1, because leading skipped bases are not penalized in the main for-loop
						else
							extended_score += 2; // increase by 2, because intermittent skipped bases are already penalized in the main for-loop
						if (extended_score >= min_score)
							return true;

					} else { // there is a mismatch

						// allow one mismatch
						mismatch_count++;
						if (mismatch_count > 1)
							break;
						//extended_score--; // no need to penalize the score, because this is already done in the main for-loop
					}
					extended_read_pos--;
					extended_gene_pos--;
				}
			}

			// extend match locally to the right
			{
				int extended_read_pos = read_pos + kmer_length;
				int extended_gene_pos = *kmer_hit + kmer_length;
				unsigned int mismatch_count = 0;
				unsigned int consecutive_mismatches = 0;
				splice_sites_t::const_iterator next_splice_site = splice_sites.lower_bound(extended_gene_pos - 1);
				while (extended_read_pos < (int) read_sequence.length() && extended_gene_pos <= gene_end) {

					// try a spliced alignment, if we run over a splice-site
					if (next_splice_site != splice_sites.end()) {
						if (extended_gene_pos - 1 > *next_splice_site)
							++next_splice_site;
						if (next_splice_site != splice_sites.end() && extended_gene_pos - 1 == *next_splice_site)
							if (align(extended_score, read_sequence, extended_read_pos, contig_sequence, extended_gene_pos, gene_start, gene_end, kmer_index, kmer_length, splice_sites, min_score, max_deletions))
								return true;
					}

					if (read_sequence[extended_read_pos] == contig_sequence[extended_gene_pos]) {

						// read and gene sequences match => increase score and go the next base
						extended_score++;
						if (extended_score >= min_score)
							return true;
						consecutive_mismatches = 0;

					} else { // there is a mismatch

						mismatch_count++;
						if (mismatch_count == 1) // when there is more than one mismatch, do another k-mer lookup
							if (max_deletions > 0 && read_sequence.length() >= 30 && // do not allow too many deletions/introns and only if the read is reasonably long
							    align(extended_score, read_sequence, extended_read_pos, contig_sequence, extended_gene_pos, gene_start, gene_end, kmer_index, kmer_length, splice_sites, min_score, max_deletions-1))
								return true;
						extended_score--; // penalize mismatch
						consecutive_mismatches++;
						if (consecutive_mismatches >= 4)
							break;
					}

					extended_read_pos++;
					extended_gene_pos++;
				}
			}
		}
	}

	// we only get here, if the read could not be aligned
	return false;
}

bool align_both_strands(const string& read_sequence, const int read_length, const int max_mate_gap, const bool breakpoints_on_same_contig, const position_t alignment_start, const position_t alignment_end, const kmer_indices_t& kmer_indices, const assembly_t& assembly, const exon_annotation_index_t& exon_annotation_index, splice_sites_by_gene_t& splice_sites_by_gene, gene_set_t& genes, const char kmer_length, const float min_align_percent, int min_score) {
	min_score = min(min_score, (int) (min_align_percent * read_sequence.size() + 0.5));
	for (gene_set_t::iterator gene = genes.begin(); gene != genes.end(); ++gene) {

		// find all splice sites in the genes
		if (splice_sites_by_gene.find(*gene) == splice_sites_by_gene.end())
			get_downstream_splice_sites(*gene, exon_annotation_index, splice_sites_by_gene[*gene]);

		// align against gene and some buffer before and after the gene (but not beyond contig boundaries)
		position_t gene_start = max((**gene).start - max_mate_gap - read_length, 0);
		position_t gene_end = min((**gene).end + max_mate_gap + read_length, (int) assembly.at((**gene).contig).size() - 1);

		// in the case of intragenic events or overlapping genes,
		// both, the donor AND the acceptor gene overlap the breakpoint
		// => we do no alignment, because this would always discard the read
		if (breakpoints_on_same_contig &&
		    (alignment_start >= gene_start && alignment_start <= gene_end ||
		     alignment_end   >= gene_start && alignment_end   <= gene_end))
			continue;

		if (align(0, read_sequence, 0, assembly.at((**gene).contig), gene_start, gene_start, gene_end, kmer_indices[(**gene).contig], kmer_length, splice_sites_by_gene.at(*gene), min_score, 1)) { // align on forward strand
			return true;
		} else { // align on reverse strand
			string reverse_complement;
			string original = read_sequence;
			dna_to_reverse_complement(original, reverse_complement);
			if (align(0, reverse_complement, 0, assembly.at((**gene).contig), gene_start, gene_start, gene_end, kmer_indices[(**gene).contig], kmer_length, splice_sites_by_gene.at(*gene), min_score, 1))
				return true;
		}
	}
	return false;
}

void count_mismappers(vector<chimeric_alignments_t::iterator>& chimeric_alignments_list, short unsigned int& mismappers, short unsigned int& total_reads, short unsigned int& supporting_reads) {
	for (auto chimeric_alignment = chimeric_alignments_list.begin(); chimeric_alignment != chimeric_alignments_list.end(); ++chimeric_alignment) {
		if ((**chimeric_alignment).second.filter == NULL) {
			total_reads++;
		} else if ((**chimeric_alignment).second.filter == FILTERS.at("mismappers")) {
			total_reads++;
			mismappers++;
			if (supporting_reads > 0)
				supporting_reads--;
		}
	}
}

// extend split read and compare against reference to check if STAR clipped prematurely (mostly due to accumulation of SNPs)
bool extend_split_read(const alignment_t& split_read, const assembly_t& assembly, const float min_align_percent) {

	// get clipped segment and the reference sequence at the position of the clipped segment
	string clipped_sequence;
	string reference_sequence;
	int clipped_count;
	if (split_read.strand == FORWARD) {
		clipped_count = min((int) split_read.preclipping(), split_read.start); // don't run over contig boundary
		clipped_sequence = split_read.sequence.substr(split_read.preclipping() - clipped_count, clipped_count);
		reference_sequence = assembly.at(split_read.contig).substr(split_read.start - clipped_count, clipped_count);
	} else {
		clipped_count = min((int) split_read.postclipping(), (int) assembly.at(split_read.contig).size() - split_read.end); // don't run over contig boundary
		clipped_sequence = split_read.sequence.substr(split_read.sequence.size() - split_read.postclipping(), clipped_count);
		reference_sequence = assembly.at(split_read.contig).substr(split_read.end, clipped_count);
	}

	// count number of matching bases between clipped segment and reference
	unsigned int matching_bases = 0;
	for (unsigned int i = 0; i < clipped_sequence.size(); ++i)
		if (clipped_sequence[i] == reference_sequence[i])
			++matching_bases;

	return matching_bases >= floor(clipped_sequence.size() * min_align_percent);
}

unsigned int filter_mismappers(fusions_t& fusions, const kmer_indices_t& kmer_indices, const char kmer_length, const assembly_t& assembly, const exon_annotation_index_t& exon_annotation_index, const float max_mismapper_fraction, const int max_mate_gap) {

	const float min_align_percent = 0.8; // allow ~1 mismatch for every 10 matches
	const int min_score = 40; // consider this score or higher a match (even if less than min_align_percent match)

	splice_sites_by_gene_t splice_sites_by_gene;

	// align discordnat mate / clipped segment in gene of origin
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.gene1 == fusion->second.gene2)
			continue; // re-aligning the read only makes sense between different genes

		if (fusion->second.filter != NULL)
			continue;

		//TODO hotfix to prevent MTAP:CDKN2B-AS1 from being removed
		if (fusion->second.gene1->name == "MTAP" && fusion->second.gene2->name == "CDKN2B-AS1")
			continue;

		// re-align split reads
		vector<chimeric_alignments_t::iterator> all_split_reads;
		all_split_reads.insert(all_split_reads.end(), fusion->second.split_read1_list.begin(), fusion->second.split_read1_list.end());
		all_split_reads.insert(all_split_reads.end(), fusion->second.split_read2_list.begin(), fusion->second.split_read2_list.end());
		for (auto chimeric_alignment = all_split_reads.begin(); chimeric_alignment != all_split_reads.end(); ++chimeric_alignment) {

			if ((**chimeric_alignment).second.filter != NULL)
				continue; // read has already been filtered

			// introduce aliases for cleaner code
			alignment_t& split_read = (**chimeric_alignment).second[SPLIT_READ];
			alignment_t& supplementary = (**chimeric_alignment).second[SUPPLEMENTARY];
			alignment_t& mate1 = (**chimeric_alignment).second[MATE1];

			if (split_read.strand == FORWARD) {
				if (extend_split_read(split_read, assembly, min_align_percent) ||
				    align_both_strands(split_read.sequence.substr(0, split_read.preclipping()), split_read.sequence.size(), max_mate_gap, fusion->second.contig1 == fusion->second.contig2, supplementary.start, supplementary.end, kmer_indices, assembly, exon_annotation_index, splice_sites_by_gene, split_read.genes, kmer_length, min_align_percent, min_score) || // clipped segment aligns to donor
				    align_both_strands(mate1.sequence.substr(mate1.preclipping()), mate1.sequence.size(), max_mate_gap, fusion->second.contig1 == fusion->second.contig2, mate1.start, mate1.end, kmer_indices, assembly, exon_annotation_index, splice_sites_by_gene, supplementary.genes, kmer_length, min_align_percent, min_score)) { // non-spliced mate aligns to acceptor
					(**chimeric_alignment).second.filter = FILTERS.at("mismappers");
				}
			} else { // split_read.strand == REVERSE
				if (extend_split_read(split_read, assembly, min_align_percent) ||
				    align_both_strands(split_read.sequence.substr(split_read.sequence.length() - split_read.postclipping()), split_read.sequence.size(), max_mate_gap, fusion->second.contig1 == fusion->second.contig2, supplementary.start, supplementary.end, kmer_indices, assembly, exon_annotation_index, splice_sites_by_gene, split_read.genes, kmer_length, min_align_percent, min_score) || // clipped segment aligns to donor
				    align_both_strands(mate1.sequence.substr(0, mate1.sequence.length() - mate1.postclipping()), mate1.sequence.size(), max_mate_gap, fusion->second.contig1 == fusion->second.contig2, mate1.start, mate1.end, kmer_indices, assembly, exon_annotation_index, splice_sites_by_gene, supplementary.genes, kmer_length, min_align_percent, min_score)) { // non-spliced mate aligns to acceptor
					(**chimeric_alignment).second.filter = FILTERS.at("mismappers");
				}
			}
		}

		// re-align discordant mates
		for (auto chimeric_alignment = fusion->second.discordant_mate_list.begin(); chimeric_alignment != fusion->second.discordant_mate_list.end(); ++chimeric_alignment) {
			if ((**chimeric_alignment).second.filter != NULL)
				continue; // read has already been filtered

			if ((**chimeric_alignment).second.size() == 2) { // discordant mates

				// introduce aliases for cleaner code
				alignment_t& mate1 = (**chimeric_alignment).second[MATE1];
				alignment_t& mate2 = (**chimeric_alignment).second[MATE2];

				if (align_both_strands(mate1.sequence, mate1.sequence.size(), max_mate_gap, fusion->second.contig1 == fusion->second.contig2, mate1.start, mate1.end, kmer_indices, assembly, exon_annotation_index, splice_sites_by_gene, mate2.genes, kmer_length, min_align_percent, min_score) ||
				    align_both_strands(mate2.sequence, mate2.sequence.size(), max_mate_gap, fusion->second.contig1 == fusion->second.contig2, mate2.start, mate2.end, kmer_indices, assembly, exon_annotation_index, splice_sites_by_gene, mate1.genes, kmer_length, min_align_percent, min_score)) {
					(**chimeric_alignment).second.filter = FILTERS.at("mismappers");
				}
			}
		}

	}

	// discard all fusions with more than XX% mismappers
	unsigned int remaining = 0;
	for (fusions_t::iterator fusion = fusions.begin(); fusion != fusions.end(); ++fusion) {

		if (fusion->second.filter != NULL)
			continue; // fusion has already been filtered

		short unsigned int total_reads = 0;
		short unsigned int mismappers = 0;
		count_mismappers(fusion->second.split_read1_list, mismappers, total_reads, fusion->second.split_reads1);
		count_mismappers(fusion->second.split_read2_list, mismappers, total_reads, fusion->second.split_reads2);
		count_mismappers(fusion->second.discordant_mate_list, mismappers, total_reads, fusion->second.discordant_mates);

		// remove fusions with mostly mismappers
		if (mismappers > 0 && mismappers >= floor(max_mismapper_fraction * total_reads))
			fusion->second.filter = FILTERS.at("mismappers");
		else
			remaining++;

	}

	return remaining;
}

