# input directories
HTSLIB := htslib
SOURCE := source
STATIC_LIBS := static_libs_centos6.10

# compiler flags
CXX := g++
CXXFLAGS := -Wall -Wno-parentheses -pthread -std=c++0x -O2

# the following two variables are needed to pass the location of headers/libraries in a bioconda build environment
CPPFLAGS := -I$(HTSLIB)/htslib
LDFLAGS := 

# the LIBS* variables define which libraries should be linked statically/dynamically
LIBS_SO := -lz -lm -lbz2 -llzma
LIBS_A := $(HTSLIB)/libhts.a

all: arriba

arriba: $(SOURCE)/arriba.cpp $(SOURCE)/annotation.o $(SOURCE)/assembly.o $(SOURCE)/options.o $(SOURCE)/read_chimeric_alignments.o $(SOURCE)/filter_uninteresting_contigs.o $(SOURCE)/filter_inconsistently_clipped.o $(SOURCE)/filter_homopolymer.o $(SOURCE)/filter_duplicates.o $(SOURCE)/read_stats.o $(SOURCE)/fusions.o $(SOURCE)/filter_proximal_read_through.o $(SOURCE)/filter_same_gene.o $(SOURCE)/filter_small_insert_size.o $(SOURCE)/filter_long_gap.o $(SOURCE)/filter_hairpin.o $(SOURCE)/filter_mismatches.o $(SOURCE)/filter_low_entropy.o $(SOURCE)/filter_relative_support.o $(SOURCE)/filter_both_intronic.o $(SOURCE)/filter_non_coding_neighbors.o $(SOURCE)/filter_intragenic_both_exonic.o $(SOURCE)/filter_min_support.o $(SOURCE)/recover_known_fusions.o $(SOURCE)/recover_both_spliced.o $(SOURCE)/filter_blacklisted_ranges.o $(SOURCE)/filter_end_to_end.o $(SOURCE)/filter_pcr_fusions.o $(SOURCE)/merge_adjacent_fusions.o $(SOURCE)/select_best.o $(SOURCE)/filter_short_anchor.o $(SOURCE)/filter_no_coverage.o $(SOURCE)/filter_homologs.o $(SOURCE)/filter_mismappers.o $(SOURCE)/recover_many_spliced.o $(SOURCE)/filter_genomic_support.o $(SOURCE)/recover_isoforms.o $(SOURCE)/output_fusions.o $(SOURCE)/read_compressed_file.o $(LIBS_A)
	$(CXX) $(CXXFLAGS) -I$(SOURCE) $(CPPFLAGS) -o arriba $^ $(LDFLAGS) $(LIBS_SO)

%.o: %.cpp $(wildcard $(SOURCE)/*.hpp)
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

$(HTSLIB)/libhts.a:
	$(MAKE) -C $(HTSLIB) CPPFLAGS="$(CPPFLAGS)" LDFLAGS="$(LDFLAGS)" libhts.a

clean:
	rm -f $(SOURCE)/*.o arriba
	$(MAKE) -C $(HTSLIB) clean

release:
	$(MAKE) LIBS_SO="" LIBS_A="$(LIBS_A) $(wildcard $(STATIC_LIBS)/*.a)" CPPFLAGS="-DHAVE_LIBDEFLATE $(CPPFLAGS) -I$(STATIC_LIBS) -I../$(STATIC_LIBS)"

bioconda:
	$(MAKE) LIBS_SO="-ldl -lhts -ldeflate $(LIBS_SO)" LIBS_A="" CPPFLAGS="-DHAVE_LIBDEFLATE $(CPPFLAGS)" LDFLAGS="$(LDFLAGS)"

