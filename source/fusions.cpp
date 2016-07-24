#include <unordered_map>
#include <tuple>
#include <vector>
#include <list>
#include <string>
#include <cmath>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "fusions.hpp"

using namespace std;

unsigned int find_fusions(chimeric_alignments_t& chimeric_alignments, fusions_t& fusions, annotation_index_t& exon_annotation_index) {

	typedef unordered_map< tuple<gene_t/*gene1*/,gene_t/*gene2*/>, vector<chimeric_alignments_t::iterator> > discordant_mates_by_gene_pair_t;
	discordant_mates_by_gene_pair_t discordant_mates_by_gene_pair; // contains the discordant mates for each pair of genes

	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {

		contig_t contig1, contig2;
		position_t breakpoint1, breakpoint2;
		direction_t direction1, direction2;
		gene_set_t genes1, genes2;
		bool exonic1, exonic2;
		position_t anchor_start1, anchor_start2;

		if (i->second.size() == 3) { // split read

			// extract info about fusion from alignments
			contig1 = i->second[SPLIT_READ].contig;
			contig2 = i->second[SUPPLEMENTARY].contig;
			breakpoint1 = (i->second[SPLIT_READ].strand == FORWARD) ? i->second[SPLIT_READ].start : i->second[SPLIT_READ].end;
			breakpoint2 = (i->second[SUPPLEMENTARY].strand == FORWARD) ? i->second[SUPPLEMENTARY].end : i->second[SUPPLEMENTARY].start;
			genes1 = i->second[SPLIT_READ].genes;
			genes2 = i->second[SUPPLEMENTARY].genes;
			direction1 = (i->second[SPLIT_READ].strand == FORWARD) ? UPSTREAM : DOWNSTREAM;
			direction2 = (i->second[SUPPLEMENTARY].strand == FORWARD) ? DOWNSTREAM : UPSTREAM;
			exonic1 = i->second[SPLIT_READ].exonic;
			exonic2 = i->second[SUPPLEMENTARY].exonic;
			anchor_start1 = (i->second[MATE1].strand == FORWARD) ? i->second[MATE1].start : i->second[MATE1].end;
			anchor_start2 = (i->second[SUPPLEMENTARY].strand == FORWARD) ? i->second[SUPPLEMENTARY].start : i->second[SUPPLEMENTARY].end;

			// make sure the breakpoint with the lower coordinate is always first
			// otherwise the same fusion could generate two entries in the fusions hashmap
			bool swapped = false;
			if (contig1 > contig2 || (contig1 == contig2 && breakpoint1 > breakpoint2)) {
				swap(contig1, contig2);
				swap(breakpoint1, breakpoint2);
				swap(genes1, genes2);
				swap(direction1, direction2);
				swap(exonic1, exonic2);
				swap(anchor_start1, anchor_start2);
				swapped = true;
			}

			// make a fusion from the given breakpoints
			bool overlap_duplicate1 = false;
			for (gene_set_t::iterator gene1 = genes1.begin(); gene1 != genes1.end(); ++gene1) {
				bool overlap_duplicate2 = false;
				for (gene_set_t::iterator gene2 = genes2.begin(); gene2 != genes2.end(); ++gene2) {

					// copy properties of supporting read to fusion
					fusion_t& fusion = fusions[make_tuple(*gene1, *gene2, contig1, contig2, breakpoint1, breakpoint2, direction1, direction2)];
					fusion.gene1 = *gene1; fusion.gene2 = *gene2;
					fusion.direction1 = direction1; fusion.direction2 = direction2;
					fusion.contig1 = contig1; fusion.contig2 = contig2;
					fusion.breakpoint1 = breakpoint1; fusion.breakpoint2 = breakpoint2;
					fusion.exonic1 = exonic1; fusion.exonic2 = exonic2;

					if (!exonic1 && !exonic2) {
						fusion.spliced1 = false;
						fusion.spliced2 = false;
					} else {
						fusion.spliced1 = is_breakpoint_spliced(*gene1, direction1, contig1, breakpoint1, exon_annotation_index);
						fusion.spliced2 = is_breakpoint_spliced(*gene2, direction2, contig2, breakpoint2, exon_annotation_index);
					}

					if (fusion.split_reads1 + fusion.split_reads2 + fusion.discordant_mates == 0) {
						if (i->second.filters.empty())
							fusion.filters.clear();
						else
							fusion.filters.insert(*i->second.filters.begin()); // TODO adjust this when we support more than one filter per read
					}
					fusion.chimeric_alignments.push_back(&(i->second));

					// expand the size of the anchor
					if (fusion.direction1 == DOWNSTREAM && (anchor_start1 < fusion.anchor_start1 || fusion.anchor_start1 == 0)) {
						fusion.anchor_start1 = anchor_start1;
					} else if (fusion.direction1 == UPSTREAM && (anchor_start1 > fusion.anchor_start1 || fusion.anchor_start1 == 0)) {
						fusion.anchor_start1 = anchor_start1;
					}
					if (fusion.direction2 == DOWNSTREAM && (anchor_start2 < fusion.anchor_start2 || fusion.anchor_start2 == 0)) {
						fusion.anchor_start2 = anchor_start2;
					} else if (fusion.direction2 == UPSTREAM && (anchor_start2 > fusion.anchor_start2 || fusion.anchor_start2 == 0)) {
						fusion.anchor_start2 = anchor_start2;
					}

					// when the breakpoint falls into a region where genes overlap,
					// mark all genes except the first as "overlap_duplicate"
					fusion.overlap_duplicate1 = overlap_duplicate1;
					fusion.overlap_duplicate2 = overlap_duplicate2;
					overlap_duplicate2 = true;

					// increase split read counters for the given fusion
					if (i->second.filters.empty())
						if (swapped)
							fusion.split_reads2++;
						else
							fusion.split_reads1++;
				}
				overlap_duplicate1 = true;
			}

		} else if (i->second.size() == 2) { // discordant mates

			// extract info about fusion from alignments
			contig1 = i->second[MATE1].contig;
			contig2 = i->second[MATE2].contig;
			breakpoint1 = (i->second[MATE1].strand == FORWARD) ? i->second[MATE1].end : i->second[MATE1].start;
			breakpoint2 = (i->second[MATE2].strand == FORWARD) ? i->second[MATE2].end : i->second[MATE2].start;
			genes1 = i->second[MATE1].genes;
			genes2 = i->second[MATE2].genes;
			direction1 = (i->second[MATE1].strand == FORWARD) ? DOWNSTREAM : UPSTREAM;
			direction2 = (i->second[MATE2].strand == FORWARD) ? DOWNSTREAM : UPSTREAM;
			exonic1 = i->second[MATE1].exonic;
			exonic2 = i->second[MATE2].exonic;
			anchor_start1 = (i->second[MATE1].strand == FORWARD) ? i->second[MATE1].start : i->second[MATE1].end;
			anchor_start2 = (i->second[MATE2].strand == FORWARD) ? i->second[MATE2].start : i->second[MATE2].end;

			// make sure the breakpoint with the lower coordinate is always first
			// otherwise the same fusion could generate two entries in the fusions hashmap
			if (contig1 > contig2 || (contig1 == contig2 && breakpoint1 > breakpoint2)) {
				swap(contig1, contig2);
				swap(breakpoint1, breakpoint2);
				swap(genes1, genes2);
				swap(direction1, direction2);
				swap(exonic1, exonic2);
				swap(anchor_start1, anchor_start2);
			}

			// make a fusion from the given breakpoints
			bool overlap_duplicate1 = false;
			for (gene_set_t::iterator gene1 = genes1.begin(); gene1 != genes1.end(); ++gene1) {
				bool overlap_duplicate2 = false;
				for (gene_set_t::iterator gene2 = genes2.begin(); gene2 != genes2.end(); ++gene2) {

					// copy properties of supporting read to fusion
					bool is_new_fusion = fusions.find(make_tuple(*gene1, *gene2, contig1, contig2, breakpoint1, breakpoint2, direction1, direction2)) == fusions.end();
					fusion_t& fusion = fusions[make_tuple(*gene1, *gene2, contig1, contig2, breakpoint1, breakpoint2, direction1, direction2)];
					fusion.gene1 = *gene1; fusion.gene2 = *gene2;
					fusion.direction1 = direction1; fusion.direction2 = direction2;
					fusion.contig1 = contig1; fusion.contig2 = contig2;
					fusion.breakpoint1 = breakpoint1; fusion.breakpoint2 = breakpoint2;
					fusion.exonic1 = exonic1; fusion.exonic2 = exonic2;
					fusion.spliced1 = false; fusion.spliced2 = false;

					if (i->second.filters.empty()) {
						fusion.filters.clear();
					} else if (is_new_fusion || !fusion.filters.empty()) {
						fusion.filters.insert(*i->second.filters.begin()); // TODO adjust this when we support more than one filter per read
					}

					// expand the size of the anchor
					if (fusion.direction1 == DOWNSTREAM && (anchor_start1 < fusion.anchor_start1 || fusion.anchor_start1 == 0)) {
						fusion.anchor_start1 = anchor_start1;
					} else if (fusion.direction1 == UPSTREAM && (anchor_start1 > fusion.anchor_start1 || fusion.anchor_start1 == 0)) {
						fusion.anchor_start1 = anchor_start1;
					}
					if (fusion.direction2 == DOWNSTREAM && (anchor_start2 < fusion.anchor_start2 || fusion.anchor_start2 == 0)) {
						fusion.anchor_start2 = anchor_start2;
					} else if (fusion.direction2 == UPSTREAM && (anchor_start2 > fusion.anchor_start2 || fusion.anchor_start2 == 0)) {
						fusion.anchor_start2 = anchor_start2;
					}

					// when the breakpoint falls into a region where genes overlap,
					// mark all genes except the first as "overlap_duplicate"
					fusion.overlap_duplicate1 = overlap_duplicate1;
					fusion.overlap_duplicate2 = overlap_duplicate2;
					overlap_duplicate2 = true;

					// store the discordant mates in a hashmap for fast lookup
					// we will need this later to find all the discordant mates supporting a given fusion
					discordant_mates_by_gene_pair[make_tuple(*gene1, *gene2)].push_back(i);
				}
				overlap_duplicate1 = true;
			}
		}
	}

	// for each fusion, count the supporting discordant mates


	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (!i->second.filters.empty())
			continue; // don't look for discordant mates, if the fusion has been filtered

		// get list of discordant mates supporting a fusion between the given gene pair
		discordant_mates_by_gene_pair_t::iterator discordant_mates = discordant_mates_by_gene_pair.find(make_tuple(i->second.gene1, i->second.gene2));
		if (discordant_mates != discordant_mates_by_gene_pair.end()) {

			// discard those discordant mates which point in the wrong direction (away from the breakpoint)
			// or which are too far away from the breakpoint (3 standard deviations of the average mate gap)
			for (auto discordant_mate = discordant_mates->second.begin(); discordant_mate != discordant_mates->second.end(); ++discordant_mate) {

				alignment_t* mate1 = &((*discordant_mate)->second[MATE1]); // introduce some aliases for cleaner code
				alignment_t* mate2 = &((*discordant_mate)->second[MATE2]);

				// make sure mate1 points to the mate with the lower coordinate
				// this ensures that the coordinate of the correct mate is compared against the coordinate of the breakpoint
				if (mate1->contig > mate2->contig || (mate1->contig == mate2->contig && mate1->start > mate2->start))
					swap(mate1, mate2); //TODO this generates fusions with supporting reads all set to 0, when the fusion is intragenic

				if (((i->second.direction1 == DOWNSTREAM && mate1->strand == FORWARD && (mate1->end-2 <= i->second.breakpoint1 && i->second.split_reads1 + i->second.split_reads2 > 0 || mate1->end-200 <= i->second.breakpoint1 && i->second.split_reads1 + i->second.split_reads2 == 0)) || (i->second.direction1 == UPSTREAM && mate1->strand == REVERSE && (mate1->start+2 >= i->second.breakpoint1 && i->second.split_reads1+i->second.split_reads2 > 0 || mate1->start+200 >= i->second.breakpoint1 && i->second.split_reads1+i->second.split_reads2 == 0))) &&
				    ((i->second.direction2 == DOWNSTREAM && mate2->strand == FORWARD && (mate2->end-2 <= i->second.breakpoint2 && i->second.split_reads1 + i->second.split_reads2 > 0 || mate2->end-200 <= i->second.breakpoint2 && i->second.split_reads1 + i->second.split_reads2 == 0)) || (i->second.direction2 == UPSTREAM && mate2->strand == REVERSE && (mate2->start+2 >= i->second.breakpoint2 && i->second.split_reads1+i->second.split_reads2 > 0 || mate2->start+200 >= i->second.breakpoint2 && i->second.split_reads1+i->second.split_reads2 == 0)))) {

					i->second.chimeric_alignments.push_back(&(**discordant_mate).second);

					if ((*discordant_mate)->second.filters.empty())
						i->second.discordant_mates++;
					else if (!i->second.filters.empty())
						i->second.filters.insert(*(*discordant_mate)->second.filters.begin()); // TODO adjust this when we support more than one filter per read

					// expand the size of the anchor
					if (i->second.direction1 == DOWNSTREAM && (mate1->start < i->second.anchor_start1 || i->second.anchor_start1 == 0)) {
						i->second.anchor_start1 = mate1->start;
					} else if (i->second.direction1 == UPSTREAM && (mate1->end > i->second.anchor_start1 || i->second.anchor_start1 == 0)) {
						i->second.anchor_start1 = mate1->end;
					}
					if (i->second.direction2 == DOWNSTREAM && (mate2->start < i->second.anchor_start2 || i->second.anchor_start2 == 0)) {
						i->second.anchor_start2 = mate2->start;
					} else if (i->second.direction2 == UPSTREAM && (mate2->end > i->second.anchor_start2 || i->second.anchor_start2 == 0)) {
						i->second.anchor_start2 = mate2->end;
					}
				}
			}
		}
	}

	// count fusions
	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i)
		if (i->second.filters.empty())
			remaining++;

	return remaining;
}

