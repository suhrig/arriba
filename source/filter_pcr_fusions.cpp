#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include "common.hpp"
#include "annotation.hpp"
#include "filter_pcr_fusions.hpp"

using namespace std;














unsigned int filter_pcr_fusions(fusions_t& fusions, const annotation_t& gene_annotation, const float max_pcr_fusion_score, const unsigned int max_exonic_breakpoints, const unsigned int max_partners_with_many_exonic_breakpoints, const unsigned int min_split_reads) {

	vector<unsigned int> exonic_breakpoint_count(gene_annotation.size()); // count the number of fusions with exonic (non-spliced) breakpoints for each gene
	vector<unsigned int> partners_with_many_exonic_breakpoints(gene_annotation.size()); // count the number of gene partners which have many exonic breakpoints for each gene
	unordered_map< tuple<gene_t/*gene1*/, gene_t/*gene2*/>, unsigned int > exonic_breakpoints_by_gene_pair; // count the number of exonic breakpoints for each gene pair
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
		if (i->second.gene1 != i->second.gene2 &&
		    !i->second.spliced1 && !i->second.spliced2 &&
		    (!i->second.overlap_duplicate1 || !i->second.overlap_duplicate2) &&
		    i->second.exonic1 && i->second.exonic2 &&
		    i->second.filters.find(FILTERS.at("merge_adjacent")) == i->second.filters.end()) {

			unsigned int split_reads = i->second.split_reads1 + i->second.split_reads2;
			if (split_reads == 0) {
				// look for split reads in the filtered reads
				for (auto chimeric_alignment = i->second.chimeric_alignments.begin(); chimeric_alignment != i->second.chimeric_alignments.end(); ++chimeric_alignment) {
					if ((**chimeric_alignment).size() == 3) {
						split_reads++;
						break; // one read is enough
					}
				}
			}

			if (split_reads > 0) {
				if (!i->second.overlap_duplicate2)
					exonic_breakpoint_count[i->second.gene1]++;
				if (!i->second.overlap_duplicate1)
					exonic_breakpoint_count[i->second.gene2]++;
				if (++exonic_breakpoints_by_gene_pair[make_tuple(i->second.gene1, i->second.gene2)] == max_exonic_breakpoints) {
					partners_with_many_exonic_breakpoints[i->second.gene1]++;
					partners_with_many_exonic_breakpoints[i->second.gene2]++;
				}
			}
		}
	}

	// calculate score which reflects the likelihood of PCR fusions for each gene
	vector<float> pcr_fusion_scores(gene_annotation.size());
	for (unsigned int i = 0; i < gene_annotation.size(); ++i) {
		pcr_fusion_scores[i] = partners_with_many_exonic_breakpoints[i];
		if (exonic_breakpoint_count[i] > 0)
			pcr_fusion_scores[i] *= log10(exonic_breakpoint_count[i]); // we take the logarithm so that the factors have about equal order of magnitude / weight
		else
			pcr_fusion_scores[i] = 0;
	}

	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (!i->second.filters.empty())
			continue; // fusion has already been filtered

//TODO shouldn't this say !spliced1 OR !spliced2
//TODO how about we don't care if it is spliced in the PCR amplified gene and only care about if it is spliced in the other gene?
		if ((!i->second.spliced1 && !i->second.spliced2 || i->second.split_reads1 + i->second.split_reads2 <= min_split_reads || pcr_fusion_scores[i->second.gene1] >= max_pcr_fusion_score && pcr_fusion_scores[i->second.gene2] >= max_pcr_fusion_score) && // discard fusions, which are not spliced / have few reads / are between genes with high PCR fusion scores
		    (i->second.exonic1 || i->second.exonic2) && // one of the breakpoints must be exonic (in theory, both should be, but there are too many exceptions)
		    pcr_fusion_scores[i->second.gene1] + pcr_fusion_scores[i->second.gene2] >= max_pcr_fusion_score && // PCR fusion score must be above threshold
		    (partners_with_many_exonic_breakpoints[i->second.gene1] >= max_partners_with_many_exonic_breakpoints || partners_with_many_exonic_breakpoints[i->second.gene2] >= max_partners_with_many_exonic_breakpoints)) { // one of the partners must have many other partners with many exonic breakpoints
			i->second.filters.insert(FILTERS.at("pcr_fusions"));
		} else {
			++remaining;
		}
	}

	return remaining;
}

