# input directories
SOURCE := source
STATIC_LIBS := $(shell mkdir -p libraries && echo libraries)

# compiler flags
CXX := g++
CXXFLAGS := -Wall -Wno-parentheses -pthread -std=c++0x -O2

# make a statically linked binary by default and a dynamically linked one for bioconda
all:
	$(MAKE) LIBS_A="$(STATIC_LIBS)/libhts.a $(STATIC_LIBS)/libdeflate.a $(STATIC_LIBS)/libz.a $(STATIC_LIBS)/libbz2.a $(STATIC_LIBS)/liblzma.a" arriba
bioconda:
	$(MAKE) LIBS_SO="-ldl -lhts -ldeflate -lz -lbz2 -llzma -lm" arriba

# make arriba executable
arriba: $(SOURCE)/arriba.cpp $(SOURCE)/annotation.o $(SOURCE)/assembly.o $(SOURCE)/options.o $(SOURCE)/read_chimeric_alignments.o $(SOURCE)/filter_duplicates.o $(SOURCE)/filter_uninteresting_contigs.o $(SOURCE)/filter_viral_contigs.o $(SOURCE)/filter_top_expressed_viral_contigs.o $(SOURCE)/filter_low_coverage_viral_contigs.o $(SOURCE)/filter_inconsistently_clipped.o $(SOURCE)/filter_homopolymer.o $(SOURCE)/read_stats.o $(SOURCE)/fusions.o $(SOURCE)/filter_proximal_read_through.o $(SOURCE)/filter_same_gene.o $(SOURCE)/filter_small_insert_size.o $(SOURCE)/filter_long_gap.o $(SOURCE)/filter_hairpin.o $(SOURCE)/filter_multimappers.o $(SOURCE)/filter_mismatches.o $(SOURCE)/filter_low_entropy.o $(SOURCE)/filter_relative_support.o $(SOURCE)/filter_both_intronic.o $(SOURCE)/filter_non_coding_neighbors.o $(SOURCE)/filter_intragenic_both_exonic.o $(SOURCE)/recover_internal_tandem_duplication.o $(SOURCE)/filter_min_support.o $(SOURCE)/recover_known_fusions.o $(SOURCE)/recover_both_spliced.o $(SOURCE)/filter_blacklisted_ranges.o $(SOURCE)/filter_end_to_end.o $(SOURCE)/filter_in_vitro.o $(SOURCE)/merge_adjacent_fusions.o $(SOURCE)/select_best.o $(SOURCE)/filter_marginal_read_through.o $(SOURCE)/filter_short_anchor.o $(SOURCE)/filter_no_coverage.o $(SOURCE)/filter_homologs.o $(SOURCE)/filter_mismappers.o $(SOURCE)/recover_many_spliced.o $(SOURCE)/filter_genomic_support.o $(SOURCE)/recover_isoforms.o $(SOURCE)/annotate_tags.o $(SOURCE)/annotate_protein_domains.o $(SOURCE)/output_fusions.o $(SOURCE)/read_compressed_file.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(SOURCE) -I$(STATIC_LIBS)/htslib -I$(STATIC_LIBS)/tsl -o arriba $^ $(LDFLAGS) $(LIBS_A) $(LIBS_SO)
%.o: %.cpp $(wildcard $(SOURCE)/*.hpp) $(LIBS_A) $(STATIC_LIBS)/tsl/htrie_map.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -I$(STATIC_LIBS)/htslib -I$(STATIC_LIBS)/tsl -o $@ $<

# download and compile dependencies for a static build
WGET := $(shell (which wget && echo "--no-check-certificate -O -") || echo "curl -k -L")
$(STATIC_LIBS)/tsl/htrie_map.h:
	$(WGET) 'https://github.com/Tessil/hat-trie/archive/v0.6.0.tar.gz' | tar -xzf - -C $(STATIC_LIBS) && \
	cp -r $(STATIC_LIBS)/hat-trie-*/include/tsl $(STATIC_LIBS)
$(STATIC_LIBS)/libdeflate.a:
	$(WGET) 'https://github.com/ebiggers/libdeflate/archive/v1.14.tar.gz' | tar -xzf - -C $(STATIC_LIBS) && \
	cd $(STATIC_LIBS)/libdeflate-*/ && $(MAKE) libdeflate.a && cp libdeflate.a libdeflate.h ..
$(STATIC_LIBS)/libz.a:
	$(WGET) 'https://zlib.net/fossils/zlib-1.3.1.tar.gz' | tar -xzf - -C $(STATIC_LIBS) && \
	cd $(STATIC_LIBS)/zlib-*/ && ./configure && $(MAKE) libz.a && cp zlib.h zconf.h libz.a ..
$(STATIC_LIBS)/libbz2.a:
	$(WGET) 'https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz' | tar -xzf - -C $(STATIC_LIBS) && \
	cd $(STATIC_LIBS)/bzip2-*/ && $(MAKE) libbz2.a bzip2 && cp libbz2.a bzlib.h ..
$(STATIC_LIBS)/liblzma.a:
	$(WGET) 'https://sourceforge.net/projects/lzmautils/files/xz-5.6.4.tar.gz' | tar -xzf - -C $(STATIC_LIBS) && \
	cd $(STATIC_LIBS)/xz-*/ && ./configure && $(MAKE) && cp -r src/liblzma/.libs/liblzma.a src/liblzma/api/lzma src/liblzma/api/lzma.h ..
$(STATIC_LIBS)/libhts.a: $(STATIC_LIBS)/libdeflate.a $(STATIC_LIBS)/libz.a $(STATIC_LIBS)/libbz2.a $(STATIC_LIBS)/liblzma.a
	$(WGET) 'https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2' | $(STATIC_LIBS)/bzip2-*/bzip2 -d -c - | tar -xf - -C $(STATIC_LIBS) && \
	cd $(STATIC_LIBS)/htslib-*/ && $(MAKE) config.h && sed -i -e 's/CURL/DEFLATE/' config.h && $(MAKE) NONCONFIGURE_OBJS="" CPPFLAGS="$(CPPFLAGS) -I.." libhts.a && cp -r libhts.a htslib ..

# cleanup routine
clean:
	rm -rf $(SOURCE)/*.o arriba $(STATIC_LIBS)

