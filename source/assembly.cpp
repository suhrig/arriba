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

void load_assembly(assembly_t& assembly, const string& fasta_file_path, const contigs_t& contigs, const vector<bool>& interesting_contigs) {

	// open FastA file
	stringstream fasta_file;
	autodecompress_file(fasta_file_path, fasta_file);

	// read line by line
	string line;
	contig_t current_contig = -1;
	while (getline(fasta_file, line)) {
		if (!line.empty()) {

			// get contig name
			if (line[0] == '>') {
				istringstream iss(line.substr(1));
				string contig_name;
				iss >> contig_name;
				if (contigs.find(removeChr(contig_name)) != contigs.end()) {
					current_contig = contigs.at(removeChr(contig_name));
				} else {
					cerr << "WARNING: unknown contig: " << contig_name << endl;
				}

			// get sequence
			} else if (current_contig != -1 && interesting_contigs[current_contig]) { // skip line if contig is undefined or not interesting
				std::transform(line.begin(), line.end(), line.begin(), (int (*)(int))std::toupper); // convert sequence to uppercase
				assembly[current_contig] += line;
			}
		}
	}

	// check if we found the sequence for all interesting contigs
	for (contigs_t::const_iterator contig = contigs.begin(); contig != contigs.end(); ++contig)
		if (interesting_contigs[contig->second])
			if (assembly.find(contig->second) == assembly.end()) {
				cerr << "ERROR: could not find sequence of contig '" << contig->first << "'" << endl;
				exit(1);
			}
}

