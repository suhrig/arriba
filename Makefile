HTSLIB := htslib-1.8
SOURCE := source
CXX := g++
CXXFLAGS := -pthread -std=c++0x -O2 -I$(SOURCE) -I$(HTSLIB)/htslib
LDFLAGS := -lz -llzma -lbz2

all: arriba

arriba: $(SOURCE)/arriba.cpp $(SOURCE)/annotation.o $(SOURCE)/assembly.o $(SOURCE)/options.o $(SOURCE)/read_chimeric_alignments.o $(SOURCE)/filter_multi_mappers.o $(SOURCE)/filter_uninteresting_contigs.o $(SOURCE)/filter_inconsistently_clipped.o $(SOURCE)/filter_homopolymer.o $(SOURCE)/filter_duplicates.o $(SOURCE)/read_stats.o $(SOURCE)/fusions.o $(SOURCE)/filter_proximal_read_through.o $(SOURCE)/filter_same_gene.o $(SOURCE)/filter_small_insert_size.o $(SOURCE)/filter_long_gap.o $(SOURCE)/filter_hairpin.o $(SOURCE)/filter_mismatches.o $(SOURCE)/filter_low_entropy.o $(SOURCE)/filter_relative_support.o $(SOURCE)/filter_both_intronic.o $(SOURCE)/filter_non_coding_neighbors.o $(SOURCE)/filter_intragenic_both_exonic.o $(SOURCE)/filter_min_support.o $(SOURCE)/recover_known_fusions.o $(SOURCE)/recover_both_spliced.o $(SOURCE)/filter_blacklisted_ranges.o $(SOURCE)/filter_end_to_end.o $(SOURCE)/filter_pcr_fusions.o $(SOURCE)/merge_adjacent_fusions.o $(SOURCE)/select_best.o $(SOURCE)/filter_short_anchor.o $(SOURCE)/filter_nonexpressed.o $(SOURCE)/filter_homologs.o $(SOURCE)/filter_mismappers.o $(SOURCE)/recover_many_spliced.o $(SOURCE)/filter_genomic_support.o $(SOURCE)/recover_isoforms.o $(SOURCE)/output_fusions.o $(SOURCE)/read_compressed_file.o $(HTSLIB)/libhts.a
	$(CXX) $(CXXFLAGS) -o arriba $^ $(LDFLAGS)

%.o: %.cpp %.hpp $(SOURCE)/common.hpp $(SOURCE)/annotation.hpp $(SOURCE)/annotation.t.hpp $(SOURCE)/assembly.hpp $(SOURCE)/read_compressed_file.hpp
	$(CXX) -c $(CXXFLAGS) -o $@ $< $(LDFLAGS)

$(HTSLIB)/libhts.a:
	$(MAKE) -C $(HTSLIB) libhts.a

clean:
	rm -f $(SOURCE)/*.o arriba
	$(MAKE) -C $(HTSLIB) clean

