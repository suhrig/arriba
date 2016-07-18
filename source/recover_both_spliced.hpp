#ifndef _RECOVER_BOTH_SPLICED_H
#define _RECOVER_BOTH_SPLICED_H 1

#include "common.hpp"
#include "annotation.hpp"

using namespace std;

unsigned int recover_both_spliced(fusions_t& fusions, const annotation_t& gene_annotation, const bool low_tumor_content);

#endif /* _RECOVER_BOTH_SPLICED_H */
