#include <algorithm>
#include <string>
#include <vector>
#include "common.hpp"
#include "filter_top_expressed_viral_contigs.hpp"

using namespace std;

struct sort_contigs_by_expression_t {
	vector<float>* expression;
	bool operator()(const contig_t x, const contig_t y) const {
		float expression_x = expression->at(x);
		float expression_y = expression->at(y);
		if (expression_x != expression_y) {
			return expression_x < expression_y;
		} else {
			return x < y; // ensures deterministic behavior in case of ties
		}
	}
};

unsigned int filter_top_expressed_viral_contigs(chimeric_alignments_t& chimeric_alignments, unsigned int top_count, const contigs_t& contigs, const string& viral_contigs, const vector<unsigned long int>& mapped_viral_reads_by_contig, const assembly_t& assembly) {

	// calculate expression of viral contigs, i.e., normalize mapped reads to contig length
	vector<float> expression_by_contig;
	expression_by_contig.reserve(mapped_viral_reads_by_contig.size());
	for (size_t contig = 0; contig < mapped_viral_reads_by_contig.size(); ++contig)
		if (assembly.find(contig) != assembly.end())
			expression_by_contig.push_back(1.0 * mapped_viral_reads_by_contig.at(contig)/assembly.at(contig).size());
		else
			expression_by_contig.push_back(0);

	// put viral contigs in vector for sorting by expression
	vector<contig_t> contigs_sorted_by_expression;
	contigs_sorted_by_expression.reserve(expression_by_contig.size());
	for (size_t contig = 0; contig < expression_by_contig.size(); ++contig)
		contigs_sorted_by_expression.push_back(contig);

	// partially sort using nth_element
	sort_contigs_by_expression_t sort_contigs_by_expression;
	sort_contigs_by_expression.expression = &expression_by_contig;
	if (top_count > mapped_viral_reads_by_contig.size())
		top_count = mapped_viral_reads_by_contig.size();
	top_count = mapped_viral_reads_by_contig.size() - top_count;
	nth_element(contigs_sorted_by_expression.begin(), contigs_sorted_by_expression.begin()+top_count, contigs_sorted_by_expression.end(), sort_contigs_by_expression);
	float min_expression_threshold = expression_by_contig[contigs_sorted_by_expression[top_count]];

	// convert viral_contigs to vector of booleans for faster lookup
	vector<bool> viral_contigs_bool(contigs.size());
	for (contigs_t::const_iterator contig = contigs.begin(); contig != contigs.end(); ++contig)
		viral_contigs_bool[contig->second] = is_interesting_contig(contig->first, viral_contigs);

	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != FILTER_none)
			continue; // the read has already been filtered

		// at least one mate must map to host genome
		for (mates_t::iterator mate = chimeric_alignment->second.begin(); mate != chimeric_alignment->second.end(); ++mate) {
			if (viral_contigs_bool[mate->contig]) {
				if (expression_by_contig[mate->contig] < min_expression_threshold || expression_by_contig[mate->contig] == 0) {
					chimeric_alignment->second.filter = FILTER_top_expressed_viral_contigs;
					goto next_read;
				}
			}
		}

		remaining++;

		next_read: continue;
	}
	return remaining;
}

