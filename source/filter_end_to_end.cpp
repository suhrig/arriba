#include "common.hpp"
#include "annotation.hpp"
#include "filter_end_to_end.hpp"

using namespace std;

unsigned int filter_end_to_end_fusions(fusions_t& fusions, const annotation_t& gene_annotation) {

	unsigned int remaining = 0;
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {
		if (!i->second.filters.empty())
			continue; // fusion has already been filtered

		if (i->second.discordant_mates + i->second.split_reads1 == 0 || i->second.discordant_mates + i->second.split_reads2 == 0 || i->second.split_reads1 + i->second.split_reads2 == 0) {
			if ((gene_annotation[i->second.gene1].is_dummy || (gene_annotation[i->second.gene1].strand == FORWARD && i->second.direction1 == UPSTREAM) || (gene_annotation[i->second.gene1].strand == REVERSE && i->second.direction1 == DOWNSTREAM)) &&
			    (gene_annotation[i->second.gene2].is_dummy || (gene_annotation[i->second.gene2].strand == FORWARD && i->second.direction2 == UPSTREAM) || (gene_annotation[i->second.gene2].strand == REVERSE && i->second.direction2 == DOWNSTREAM))) {
				i->second.filters.insert(FILTERS.at("end_to_end"));
				continue;
			}
		}

		++remaining;
	}

	return remaining;
}

