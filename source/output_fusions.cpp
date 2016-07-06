#include <iostream>
#include <fstream>
#include <functional>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include "common.hpp"
#include "annotation.hpp"
#include "output_fusions.hpp"

using namespace std;

// try to determine which gene makes the 5' end of the transcript by
// looking at which promoter likely drives transcription
enum transcript_start_t { TRANSCRIPT_START_GENE1, TRANSCRIPT_START_GENE2, TRANSCRIPT_START_AMBIGUOUS };
transcript_start_t get_start_of_transcript(const fusion_t& fusion, const annotation_t& gene_annotation) {
//TODO inspect splice motifs in unclear cases
	if (gene_annotation[fusion.gene1].is_dummy) { // if gene1 is a dummy gene, then gene2 has priority; otherwise gene1 has priority
		// TODO this branch needs a make-over when we support strand-specific libraries
		if (gene_annotation[fusion.gene2].strand == FORWARD && fusion.direction2 == DOWNSTREAM) { // transcript = gene2(+) -> gene1(+/-)
			return TRANSCRIPT_START_GENE2;
		} else if (gene_annotation[fusion.gene2].strand == REVERSE && fusion.direction2 == UPSTREAM) { // transcript = gene2(-) -> gene1(+/-)
			return TRANSCRIPT_START_GENE2;
		} else if (gene_annotation[fusion.gene2].strand == FORWARD && gene_annotation[fusion.gene1].contig == gene_annotation[fusion.gene2].contig && gene_annotation[fusion.gene1].end < gene_annotation[fusion.gene2].start) { // potential read-through fusion with upstream unannotated exon
			return TRANSCRIPT_START_GENE1;
		} else if (gene_annotation[fusion.gene2].strand == REVERSE && gene_annotation[fusion.gene1].contig == gene_annotation[fusion.gene2].contig && gene_annotation[fusion.gene1].start > gene_annotation[fusion.gene2].end) { // potential read-through fusion with downstream unannotated exon
			return TRANSCRIPT_START_GENE1;
		} else { // ambiguous, since orientation of dummy gene is unclear
			return TRANSCRIPT_START_AMBIGUOUS;
		}
	} else if (!fusion.exonic1 && fusion.exonic2 || // if breakpoint1 is intronic and breakpoint2 isn't, then gene2 has priority
	           !fusion.spliced1 && fusion.spliced2) { // if both breakpoints are exonic, but only breakpoint2 is at a splice site, then gene2 has priority
		if (gene_annotation[fusion.gene2].strand == FORWARD && fusion.direction2 == DOWNSTREAM) { // transcript = gene2(+) -> gene1(+/-)
			return TRANSCRIPT_START_GENE2;
		} else if (gene_annotation[fusion.gene2].strand == REVERSE && fusion.direction2 == UPSTREAM) { // transcript = gene2(-) -> gene1(+/-)
			return TRANSCRIPT_START_GENE2;
		} else if (gene_annotation[fusion.gene1].strand == FORWARD && fusion.direction1 == DOWNSTREAM) { // transcript = gene1(+) -> gene2(+/-)
			return TRANSCRIPT_START_GENE1;
		} else if (gene_annotation[fusion.gene1].strand == REVERSE && fusion.direction1 == UPSTREAM) { // transcript = gene1(-) -> gene2(+/-)
			return TRANSCRIPT_START_GENE1;
		} else { // end-to-end-fused genes
			return TRANSCRIPT_START_AMBIGUOUS;
		}
	// in all other cases gene1 has priority
	} else if (gene_annotation[fusion.gene1].strand == FORWARD && fusion.direction1 == DOWNSTREAM) { // transcript = gene1(+) -> gene2(+/-)
		return TRANSCRIPT_START_GENE1;
	} else if (gene_annotation[fusion.gene1].strand == REVERSE && fusion.direction1 == UPSTREAM) { // transcript = gene1(-) -> gene2(+/-)
		return TRANSCRIPT_START_GENE1;
	} else if (gene_annotation[fusion.gene2].strand == FORWARD && fusion.direction2 == DOWNSTREAM) { // transcript = gene2(+) -> gene1(+/-)
		return TRANSCRIPT_START_GENE2;
	} else if (gene_annotation[fusion.gene2].strand == REVERSE && fusion.direction2 == UPSTREAM) { // transcript = gene2(-) -> gene1(+/-)
		return TRANSCRIPT_START_GENE2;
	} else { // end-to-end-fused genes
		return TRANSCRIPT_START_AMBIGUOUS;
	}
}

bool get_dna_sequence_by_region(const position_t breakpoint, const direction_t direction, const gene_t gene, const annotation_t& gene_annotation, annotation_index_t& exon_annotation_index, const unsigned int length, string& sequence) {

	position_t start;
	if (direction == UPSTREAM) {
		if (breakpoint + 2 < gene_annotation[gene].start || breakpoint > gene_annotation[gene].end || breakpoint + length > gene_annotation[gene].end)
			return false; // we cannot get the sequence outside the gene
		contig_annotation_index_t::iterator next_exon_boundary = exon_annotation_index[gene_annotation[gene].contig].upper_bound(breakpoint);
		if (next_exon_boundary != exon_annotation_index[gene_annotation[gene].contig].end() && next_exon_boundary->first < breakpoint + length) { // breakpoint is close to an exon end
			if (next_exon_boundary->first <= breakpoint + 2 && next_exon_boundary->second.find(gene) == next_exon_boundary->second.end()) { // breakpoint is close to the 5' end of an exon
				next_exon_boundary = exon_annotation_index[gene_annotation[gene].contig].upper_bound(breakpoint + 3);
				if (next_exon_boundary != exon_annotation_index[gene_annotation[gene].contig].end() && next_exon_boundary->first < breakpoint + length) // breakpoint is also close to the 3' end of an exon
					return false; // breakpoint is too close to exon boundary
			} else {
				return false; // breakpoint is too close to exon boundary
			} 
		}
		start = breakpoint;
	} else { // direction == DOWNSTREAM
		if (breakpoint + 2 < gene_annotation[gene].start || breakpoint - 2 > gene_annotation[gene].end || gene_annotation[gene].start + length > breakpoint)
			return false; // we cannot get the sequence before the start
		contig_annotation_index_t::iterator next_exon_boundary = exon_annotation_index[gene_annotation[gene].contig].upper_bound(breakpoint - length);
		if (next_exon_boundary != exon_annotation_index[gene_annotation[gene].contig].end() && next_exon_boundary->first < breakpoint) {
			if (next_exon_boundary->first < breakpoint - 2 || next_exon_boundary->second.find(gene) == next_exon_boundary->second.end())
				return false; // breakpoint is too close to exon boundary
		}
		start = breakpoint - length;
	}

	start = max(0, start - gene_annotation[gene].start);
	if (start + length > gene_annotation[gene].sequence.length())
		return false;

	sequence = gene_annotation[gene].sequence.substr(start, length);
	std::transform(sequence.begin(), sequence.end(), sequence.begin(), (int (*)(int))std::toupper);

	return true;
}

string get_fusion_transcript_sequence(const fusion_t& fusion, annotation_t& gene_annotation, annotation_index_t& exon_annotation_index, const unsigned int length, const transcript_start_t transcript_start) {

	// we can only construct a fusion transcript, if we have split reads
	if (fusion.split_reads1 + fusion.split_reads2 == 0) {
		// look for discarded split-reads
		bool has_discarded_split_reads = false;
		for (auto i = fusion.chimeric_alignments.begin(); i != fusion.chimeric_alignments.end() && !has_discarded_split_reads; ++i)
			if ((**i).size() == 3 && (**i).filters.find(FILTERS.at("same_gene")) == (**i).filters.end())
				has_discarded_split_reads = true;
		if (!has_discarded_split_reads)
			return ".";
	}

	// get the sequence next to breakpoint1
	string sequence1;
	if (!get_dna_sequence_by_region(fusion.breakpoint1, fusion.direction1, fusion.gene1, gene_annotation, exon_annotation_index, length, sequence1))
		return ".";

	// get the sequence next to breakpoint2
	string sequence2;
	if (!get_dna_sequence_by_region(fusion.breakpoint2, fusion.direction2, fusion.gene2, gene_annotation, exon_annotation_index, length, sequence2))
		return ".";

	// try to determine the strand from the orientation of the gene
	// and make the reverse complement, if necessary
	if (transcript_start == TRANSCRIPT_START_GENE1) {
		return ((gene_annotation[fusion.gene1].strand == FORWARD) ? sequence1 : dna_to_reverse_complement(sequence1)) + "|" + ((fusion.direction2 == UPSTREAM) ? sequence2 : dna_to_reverse_complement(sequence2));
	} else if (transcript_start == TRANSCRIPT_START_GENE2) {
		return ((gene_annotation[fusion.gene2].strand == FORWARD) ? sequence2 : dna_to_reverse_complement(sequence2)) + "|" + ((fusion.direction1 == UPSTREAM) ? sequence1 : dna_to_reverse_complement(sequence1));
	} else { // start of gene is ambiguous
		return ".";
	}
}

bool sort_fusions_by_support(const fusion_t* x, const fusion_t* y) {
	if (x->supporting_reads() != y->supporting_reads())
		return x->supporting_reads() > y->supporting_reads();
	else if (x->evalue != y->evalue)
		return x->evalue < y->evalue;
	else
		return x->gene1 + x->gene2 < y->gene1 + y->gene2; // this does not really sort, it only ensures that fusions between the
		                                                  // same pair of genes are grouped together, if e-value and supporting reads are equal
}

typedef tuple<gene_t /*gene1*/, gene_t /*gene2*/> gene_pair_t;
class sort_fusions_by_rank_of_best_t {
public:
	unordered_map< gene_pair_t, fusion_t* > best;
	bool operator()(const fusion_t* x, const fusion_t* y) {
		fusion_t* best_x = best.at(make_tuple(x->gene1, x->gene2));
		fusion_t* best_y = best.at(make_tuple(y->gene1, y->gene2));
		if (best_x != best_y)
			return sort_fusions_by_support(best_x, best_y);
		else
			return sort_fusions_by_support(x, y);
	}
};

string gene_to_name(const gene_t gene, const contig_t contig, const position_t breakpoint, const annotation_t& gene_annotation, annotation_index_t& gene_annotation_index) {
	// if the gene is not a dummy gene (intergenic region), simply return the name of the gene
	if (!gene_annotation[gene].is_dummy) {
		return gene_annotation[gene].name;
	} else { // is a dummy gene => find flanking genes and report distances to them

		string result;

		// lookup position in gene annotation index
		contig_annotation_index_t::iterator index_hit2 = gene_annotation_index[contig].lower_bound(breakpoint);
		contig_annotation_index_t::reverse_iterator index_hit1(index_hit2);

		// go upstream until we find a non-dummy gene
		while (index_hit1 != gene_annotation_index[contig].rend() && (index_hit1->second.empty() || gene_annotation[*index_hit1->second.begin()].is_dummy))
			++index_hit1;

		// append upstream flanking genes with distances to gene name
		if (index_hit1 != gene_annotation_index[contig].rend()) {
			for (gene_set_t::iterator gene = index_hit1->second.begin(); gene != index_hit1->second.end(); ++gene) {
				if (!gene_annotation[*gene].is_dummy) {
					if (!result.empty())
						result += ",";
					result += gene_annotation[*gene].name + "(" + to_string(breakpoint - gene_annotation[*gene].end) + ")";
				}
			}
		}

		// go downstream until we find a non-dummy gene
		while (index_hit2 != gene_annotation_index[contig].end() && (index_hit2->second.empty() || gene_annotation[*index_hit2->second.begin()].is_dummy))
			++index_hit2;

		// append downstream flanking genes with distances to gene name
		if (index_hit2 != gene_annotation_index[contig].end()) {
			for (gene_set_t::iterator gene = index_hit2->second.begin(); gene != index_hit2->second.end(); ++gene) {
				if (!gene_annotation[*gene].is_dummy) {
					if (!result.empty())
						result += ",";
					result += gene_annotation[*gene].name + "(" + to_string(gene_annotation[*gene].start - breakpoint) + ")";
				}
			}
		}

		return result;
	}
}

void write_fusions_to_file(fusions_t& fusions, const string& output_file, annotation_t& gene_annotation, annotation_index_t& gene_annotation_index, annotation_index_t& exon_annotation_index, vector<string> contigs_by_id, const bool print_supporting_reads, const bool write_discarded_fusions) {
//TODO add "chr", if necessary
//TODO add fusion_strand
//TODO update help after changing columns

	// make a vector of pointers to all fusions
	// the vector will hold the fusions in sorted order
	vector<fusion_t*> sorted_fusions;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
		if (write_discarded_fusions || i->second.filters.empty())
			sorted_fusions.push_back(&(i->second));
	}

	// don't sort the discarded fusions
	if (!write_discarded_fusions) {

		// sometimes we report multiple breakpoints for the same pair of genes
		// in such cases, we want to make sure that all breakpoints of the same pair of genes
		// are grouped together in the output file, even if some of them have very low support
		// all such breakpoints will get the rank of the best ranking breakpoints
		// => find out what the best ranking breakpoints are for each pair of genes
		//    we store these breakpoints in a sorting object (sort_fusions_by_rank_of_best),
		//    which we use as a parameter to the sort function later
		sort_fusions_by_rank_of_best_t sort_fusions_by_rank_of_best;
		for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
			if (!i->second.filters.empty())
				continue; // skip discarded fusions

			gene_pair_t gene_pair = make_tuple(i->second.gene1, i->second.gene2);
			if (sort_fusions_by_rank_of_best.best.find(gene_pair) == sort_fusions_by_rank_of_best.best.end()) {
				sort_fusions_by_rank_of_best.best[gene_pair] = &(i->second);
			} else {
				fusion_t*& current_best = sort_fusions_by_rank_of_best.best[gene_pair];
				if (sort_fusions_by_support(&(i->second), current_best))
					current_best = &(i->second);
			}
		}

		// sort all gene pairs by the rank of the best scoring breakpoints of a given gene pair
		sort(sorted_fusions.begin(), sorted_fusions.end(), ref(sort_fusions_by_rank_of_best));
	}

	// write sorted list to file
	ofstream out(output_file);
	if (!out.is_open()) {
		cerr << "ERROR: Failed to open output file '" << output_file << "'." << endl;
		exit(1);
	}
	out << "#gene1\tgene2\tstrand1\tstrand2\tbreakpoint1\tbreakpoint2\tsite1\tsite2\tdirection1\tdirection2\tsplit_reads1\tsplit_reads2\tdiscordant_mates\te-value\tfilter\tfusion_transcript\tread_identifiers" << endl;
	for (vector<fusion_t*>::iterator i = sorted_fusions.begin(); i != sorted_fusions.end(); ++i) {

		if ((**i).filters.empty() == write_discarded_fusions) // write either filtered or unfiltered fusions, but not both
			continue;

		// write the gene which likely makes the 5' end of the transcript first
		transcript_start_t transcript_start = get_start_of_transcript(**i, gene_annotation);

		gene_t gene1 = (**i).gene1; gene_t gene2 = (**i).gene2;
		contig_t contig1 = (**i).contig1; contig_t contig2 = (**i).contig2;
		position_t breakpoint1 = (**i).breakpoint1; position_t breakpoint2 = (**i).breakpoint2;
		direction_t direction1 = (**i).direction1; direction_t direction2 = (**i).direction2;
		unsigned int split_reads1 = (**i).split_reads1; unsigned int split_reads2 = (**i).split_reads2;
		string site1, site2;
		if (gene_annotation[(**i).gene1].is_dummy) {
			site1 = "intergenic";
		} else if ((**i).spliced1) {
			site1 = "splice-site";
		} else if ((**i).exonic1) {
			site1 = "exonic";
		} else {
			site1 = "intronic";
		}
		if (gene_annotation[(**i).gene2].is_dummy) {
			site2 = "intergenic";
		} else if ((**i).spliced2) {
			site2 = "splice-site";
		} else if ((**i).exonic2) {
			site2 = "exonic";
		} else {
			site2 = "intronic";
		}

		if (transcript_start == TRANSCRIPT_START_GENE2) {
			swap(gene1, gene2);
			swap(direction1, direction2);
			swap(contig1, contig2);
			swap(breakpoint1, breakpoint2);
			swap(site1, site2);
			swap(split_reads1, split_reads2);
		}

		out << gene_to_name(gene1, contig1, breakpoint1, gene_annotation, gene_annotation_index) << "\t" << gene_to_name(gene2, contig2, breakpoint2, gene_annotation, gene_annotation_index) << "\t"
		    << ((gene_annotation[gene1].is_dummy) ? "." : ((gene_annotation[gene1].strand == FORWARD) ? "+" : "-")) << "\t" << ((gene_annotation[gene2].is_dummy) ? "." : ((gene_annotation[gene2].strand == FORWARD) ? "+" : "-")) << "\t"
		    << contigs_by_id[contig1] << ":" << (breakpoint1+1) << "\t" << contigs_by_id[contig2] << ":" << (breakpoint2+1) << "\t"
		    << site1 << "\t" << site2 << "\t"
		    << ((direction1 == UPSTREAM) ? "upstream" : "downstream") << "\t" << ((direction2 == UPSTREAM) ? "upstream" : "downstream") << "\t"
		    << split_reads1 << "\t" << split_reads2 << "\t" << (**i).discordant_mates << "\t" << (**i).evalue;

		// count the number of reads discarded by a given filter
		map<string,unsigned int> filters;
		for (auto filter = (**i).filters.begin(); filter != (**i).filters.end(); ++filter)
			filters[**filter] = 0;
		for (auto chimeric_alignment = (**i).chimeric_alignments.begin(); chimeric_alignment != (**i).chimeric_alignments.end(); ++chimeric_alignment) {
			if (!(**chimeric_alignment).filters.empty())
				filters[**(**chimeric_alignment).filters.begin()]++; // adjust this when we support more than one filter per read
		}

		// output filters
		out << "\t";
		if (filters.empty()) {
			out << ".";
		} else {
			for (auto filter = filters.begin(); filter != filters.end(); ++filter) {
				if (filter != filters.begin())
					out << ",";
				out << filter->first;
				if (filter->second != 0)
					out << "(" << filter->second << ")";
			}
		}

		// print a fusion-spanning sequence
		out << "\t";
		if (!write_discarded_fusions) {
			out << get_fusion_transcript_sequence(**i, gene_annotation, exon_annotation_index, 30, transcript_start);
		} else {
			out << ".";
		}

		// if requested, print identifiers of supporting reads
		out << "\t";
		if (print_supporting_reads) {
			for (auto j = (**i).chimeric_alignments.begin(); j != (**i).chimeric_alignments.end(); ++j) {
				if (j != (**i).chimeric_alignments.begin())
					out << ",";
				out << (**j).name;
			}
		} else {
			out << ".";
		}

		out << endl;
	}
	out.close();
}

