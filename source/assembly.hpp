#ifndef _ASSEMBLY_H
#define _ASSEMBLY_H 1

#include <string>
#include <vector>
#include "common.hpp"

using namespace std;

void load_assembly(assembly_t& assembly, const string& fasta_file_path, const contigs_t& contigs, const vector<bool>& interesting_contigs);

#endif /* _ASSEMBLY_H */
