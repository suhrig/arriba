#include <algorithm>
#include <vector>
#include "common.hpp"
#include "filter_top_expressed_viral_contigs.hpp"
#include "filter_mismappers.hpp"

using namespace std;

struct sort_contigs_by_expression_t {
	vector<float>* expression;
	bool operator()(const contig_t x, const contig_t y) const {
		float expression_x = expression->at(x);
		float expression_y = expression->at(y);
		if (expression_x != expression_y) {
			return expression_x > expression_y;
		} else {
			return x > y; // ensures deterministic behavior in case of ties
		}
	}
};

// determine if viruses are related based on fraction of shared kmers in their genome
bool related_viral_strains(const string& virus1, const string& virus2) {

	// choose smaller of the two viruses for making a list of its kmers
	const string* small_virus = &virus1;
	const string* big_virus = &virus2;
	if (small_virus->size() > big_virus->size())
		swap(small_virus, big_virus);

	// get all kmers of smaller virus
	const char kmer_length = 12;
	map<kmer_as_int_t, unsigned int/*non-zero if shared*/> small_virus_kmers;
	for (size_t i = 0; i + kmer_length <= small_virus->size(); i++)
		small_virus_kmers[kmer_to_int(*small_virus, i, kmer_length)] = 0;

	// calculate fraction of kmers of small virus also found in big virus
	unsigned int shared_kmers = 0;
	const unsigned int min_shared_kmers = small_virus_kmers.size() / 10; // consider viruses related if at least this fraction of kmers is shared
	for (size_t i = 0; i + kmer_length <= big_virus->size(); i++) {
		auto small_virus_kmer_count = small_virus_kmers.find(kmer_to_int(*big_virus, i, kmer_length));
		if (small_virus_kmer_count != small_virus_kmers.end() && small_virus_kmer_count->second++ == 0) // don't count repetitive kmers more than once
			if (++shared_kmers >= min_shared_kmers)
				return true;
	}

	// we only get here if the viruses don't share many kmers
	return false;
}

unsigned int filter_top_expressed_viral_contigs(chimeric_alignments_t& chimeric_alignments, unsigned int top_count, const vector<bool>& viral_contigs, const vector<bool>& interesting_contigs, const vector<unsigned long int>& mapped_viral_reads_by_contig, const assembly_t& assembly) {

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

	// sort viral contigs by expression
	sort_contigs_by_expression_t sort_contigs_by_expression;
	sort_contigs_by_expression.expression = &expression_by_contig;
	sort(contigs_sorted_by_expression.begin(), contigs_sorted_by_expression.end(), sort_contigs_by_expression);

	// merge related viral strains and count them as one when selecting the top N expressed viruses
	unsigned int corrected_top_count = 0;
	for (unsigned int i = 1; i < contigs_sorted_by_expression.size() && expression_by_contig[contigs_sorted_by_expression[i]] > 0 && top_count > 0; ++i) {
		corrected_top_count++;
		if (assembly.find(contigs_sorted_by_expression[i]) == assembly.end() ||
		    assembly.find(contigs_sorted_by_expression[i-1]) == assembly.end() ||
		    !related_viral_strains(assembly.at(contigs_sorted_by_expression[i]), assembly.at(contigs_sorted_by_expression[i-1])))
			top_count--;
	}
	if (corrected_top_count != 0)
		corrected_top_count--;
	float min_expression_threshold = expression_by_contig[contigs_sorted_by_expression[corrected_top_count]];

	// viral integration into the host genome is mostly random and many integration sites are intergenic
	// therefore, the ratio of intergenic to genic integration sites should be high for infiltrating viruses (in contrast to alignment artifacts)
	// => if there are viruses with a high intergenic-to-genic ratio, we do not apply this filter to them
	float min_fraction_of_intergenic_integration_sites = 0.33;
	unsigned int top_count_for_viruses_with_high_fraction_of_intergenic_integration_sites = 50; // don't apply this filter to the top 50 highest expressed viruses
	if (top_count_for_viruses_with_high_fraction_of_intergenic_integration_sites > mapped_viral_reads_by_contig.size())
		top_count_for_viruses_with_high_fraction_of_intergenic_integration_sites = mapped_viral_reads_by_contig.size();
	top_count_for_viruses_with_high_fraction_of_intergenic_integration_sites = mapped_viral_reads_by_contig.size() - top_count_for_viruses_with_high_fraction_of_intergenic_integration_sites;
	float min_expression_threshold_for_viruses_with_high_fraction_of_intergenic_integration_sites = expression_by_contig[contigs_sorted_by_expression[top_count_for_viruses_with_high_fraction_of_intergenic_integration_sites]];
	vector<gene_set_t> integration_sites_by_virus(viral_contigs.size());
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {
		alignment_t* viral_mapped_read = NULL;
		alignment_t* host_mapped_read = NULL;
		if (viral_contigs[chimeric_alignment->second[MATE1].contig]) {
			viral_mapped_read = &chimeric_alignment->second[MATE1];
		} else if (interesting_contigs[chimeric_alignment->second[MATE1].contig]) {
			host_mapped_read = &chimeric_alignment->second[MATE1];
		}
		unsigned int mate2 = (chimeric_alignment->second.size() == 3) ? SUPPLEMENTARY : MATE2;
		if (viral_contigs[chimeric_alignment->second[mate2].contig]) {
			viral_mapped_read = &chimeric_alignment->second[mate2];
		} else if (interesting_contigs[chimeric_alignment->second[mate2].contig]) {
			host_mapped_read = &chimeric_alignment->second[mate2];
		}
		if (viral_mapped_read != NULL && host_mapped_read != NULL)
			integration_sites_by_virus[viral_mapped_read->contig].insert(host_mapped_read->genes.begin(), host_mapped_read->genes.end());
	}
	vector<float> fraction_of_intergenic_integration_sites_by_virus(integration_sites_by_virus.size());
	for (contig_t contig = 0; contig < integration_sites_by_virus.size(); ++contig) {
		unsigned int intergenic = 0;
		unsigned int genic = 0;
		for (auto gene = integration_sites_by_virus[contig].begin(); gene != integration_sites_by_virus[contig].end(); ++gene) {
			if ((**gene).is_dummy) {
				intergenic++;
			} else {
				genic++;
			}
		}
		if (intergenic > 0)
			fraction_of_intergenic_integration_sites_by_virus[contig] = 1.0 * intergenic / (genic + intergenic);
	}

	unsigned int remaining = 0;
	for (chimeric_alignments_t::iterator chimeric_alignment = chimeric_alignments.begin(); chimeric_alignment != chimeric_alignments.end(); ++chimeric_alignment) {

		if (chimeric_alignment->second.filter != FILTER_none)
			continue; // the read has already been filtered

		// at least one mate must map to host genome
		for (mates_t::iterator mate = chimeric_alignment->second.begin(); mate != chimeric_alignment->second.end(); ++mate) {
			if (viral_contigs[mate->contig]) {
				if (expression_by_contig[mate->contig] == 0 || expression_by_contig[mate->contig] < min_expression_threshold) {
					if (fraction_of_intergenic_integration_sites_by_virus[mate->contig] < min_fraction_of_intergenic_integration_sites ||
					    expression_by_contig[mate->contig] == 0 ||
					    expression_by_contig[mate->contig] < min_expression_threshold_for_viruses_with_high_fraction_of_intergenic_integration_sites) {
						chimeric_alignment->second.filter = FILTER_top_expressed_viral_contigs;
						goto next_read;
					}
				}
			}
		}

		remaining++;

		next_read: continue;
	}
	return remaining;
}

