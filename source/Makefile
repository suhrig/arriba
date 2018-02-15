CXX := g++
CXXFLAGS := -pthread -std=c++0x -O2 -I. -I../samtools-1.3/htslib-1.3/htslib
LDFLAGS := -lz

all: arriba extract_reads

ARRIBA_OBJECTS := annotation.o assembly.o options_arriba.o options.o read_chimeric_alignments.o filter_multi_mappers.o filter_uninteresting_contigs.o filter_inconsistently_clipped.o filter_homopolymer.o filter_duplicates.o read_stats.o fusions.o filter_proximal_read_through.o filter_same_gene.o filter_small_insert_size.o filter_long_gap.o filter_hairpin.o filter_mismatches.o filter_low_entropy.o filter_relative_support.o filter_both_intronic.o filter_both_novel.o filter_intragenic_both_exonic.o filter_min_support.o recover_known_fusions.o recover_both_spliced.o filter_blacklisted_ranges.o filter_end_to_end.o filter_pcr_fusions.o merge_adjacent_fusions.o select_best.o filter_short_anchor.o filter_nonexpressed.o filter_homologs.o filter_mismappers.o recover_many_spliced.o filter_genomic_support.o recover_isoforms.o output_fusions.o read_compressed_file.o ../samtools-1.3/libbam.a ../samtools-1.3/htslib-1.3/libhts.a
arriba: $(ARRIBA_OBJECTS) arriba.cpp
	$(CXX) $(CXXFLAGS) -o arriba arriba.cpp $(ARRIBA_OBJECTS) $(LDFLAGS)

EXTRACT_READS_OBJECTS := annotation.o assembly.o options.o options_extract_reads.o extract_reads_chimeric.o extract_reads_fastq.o extract_reads_read_through_fusion.o read_compressed_file.o ../samtools-1.3/libbam.a ../samtools-1.3/htslib-1.3/libhts.a
extract_reads: $(EXTRACT_READS_OBJECTS) extract_reads.cpp
	$(CXX) $(CXXFLAGS) -o extract_reads extract_reads.cpp $(EXTRACT_READS_OBJECTS) $(LDFLAGS)

%.o: %.cpp %.hpp common.hpp annotation.hpp annotation.t.hpp assembly.hpp read_compressed_file.hpp
	$(CXX) -c $(CXXFLAGS) $< $(LDFLAGS)

../samtools-1.3/htslib-1.3/libhts.a:
	$(MAKE) -C ../samtools-1.3 htslib-1.3/libhts.a

../samtools-1.3/libbam.a:
	$(MAKE) -C ../samtools-1.3 libbam.a

clean:
	rm -f *.o ../arriba ../extract_reads

