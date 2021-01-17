#include <climits>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "read_compressed_file.hpp"
#include "assembly.hpp"

using namespace std;

void dna_to_reverse_complement(const string& dna, string& reverse_complement) {
	if (!reverse_complement.empty())
		reverse_complement.clear();
	reverse_complement.reserve(dna.length());
	for (string::const_reverse_iterator i = dna.rbegin(); i != dna.rend(); ++i)
		reverse_complement += dna_to_complement(*i);
}

string dna_to_reverse_complement(const string& dna) {
	string reverse_complement;
	dna_to_reverse_complement(dna, reverse_complement);
	return reverse_complement;
}

void load_assembly(assembly_t& assembly, const string& fasta_file_path, contigs_t& contigs, vector<string>& original_contig_names, const string& interesting_contigs) {

	// read FastA file line by line
	autodecompress_file_t fasta_file(fasta_file_path);
	string line;
	contig_t current_contig = USHRT_MAX;
	while (fasta_file.getline(line)) {
		if (!line.empty()) {

			// get contig name
			if (line[0] == '>') {
				crash(contigs.size() == USHRT_MAX - 1, "too many contigs");
				istringstream iss(line.substr(1));
				string contig_name;
				iss >> contig_name;
				pair<contigs_t::iterator,bool> new_contig = contigs.insert(pair<string,contig_t>(removeChr(contig_name), contigs.size()));
				current_contig = new_contig.first->second;
				if (original_contig_names.size() < contigs.size())
					original_contig_names.resize(contigs.size());
				original_contig_names[current_contig] = contig_name;
				if (!is_interesting_contig(contig_name, interesting_contigs))
					current_contig = USHRT_MAX; // skip uninteresting contigs

			// get sequence
			} else if (current_contig != USHRT_MAX) { // skip line if contig is undefined or not interesting
				std::transform(line.begin(), line.end(), line.begin(), (int (*)(int))std::toupper); // convert sequence to uppercase
				assembly[current_contig] += line;
			}
		}
	}
}

