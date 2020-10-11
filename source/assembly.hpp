#ifndef _ASSEMBLY_H
#define _ASSEMBLY_H 1

#include <string>
#include "common.hpp"

using namespace std;

inline char dna_to_complement(const char dna) {
	if (dna == 'a') { return 't'; }
	else if (dna == 't') { return 'a'; }
	else if (dna == 'c') { return 'g'; }
	else if (dna == 'g') { return 'c'; }
	else if (dna == 'A') { return 'T'; }
	else if (dna == 'T') { return 'A'; }
	else if (dna == 'C') { return 'G'; }
	else if (dna == 'G') { return 'C'; }
	else if (dna == '[') { return ']'; }
	else if (dna == ']') { return '['; }
	else return dna;
}

void dna_to_reverse_complement(const string& dna, string& reverse_complement);

string dna_to_reverse_complement(const string& dna);

void load_assembly(assembly_t& assembly, const string& fasta_file_path, contigs_t& contigs, vector<string>& original_contig_names, const string& interesting_contigs);

#endif /* _ASSEMBLY_H */
