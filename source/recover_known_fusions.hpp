#ifndef _RECOVER_KNOWN_FUSIONS_H
#define _RECOVER_KNOWN_FUSIONS_H 1

#include <string>
#include <unordered_map>
#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned int recover_known_fusions(fusions_t& fusions, const string& known_fusions_file_path, const unordered_map<string,gene_t>& genes, const bool low_tumor_content);

#endif /* _RECOVER_KNOWN_FUSIONS_H */
