#ifndef _EXTRACT_READS_FASTQ_H
#define _EXTRACT_READS_FASTQ_H 1

#include <fstream>
#include "sam.h"

using namespace std;

void extract_fastq(ofstream& fastq_file1, ofstream& fastq_file2, bam1_t* read1, bam1_t* read2, const uint32_t min_clipped_length);

#endif /* _EXTRACT_READS_FASTQ_H */
