bin_PROGRAMS = phylonium
phylonium_SOURCES = phylonium.cxx io.cxx io.h sequence.cxx sequence.h esa.cxx esa.h process.cxx process.h evo_model.cxx evo_model.h global.h
phylonium_CPPFLAGS = $(OPENMP_CFLAGS) -Wall -Wextra -I$(top_srcdir)/libs
phylonium_CXXFLAGS = $(OPENMP_CXXFLAGS) -O2 -ggdb
phylonium_LDADD = \
	$(top_builddir)/libs/libcompat.a \
	$(top_builddir)/libs/libpfasta.a \
	$(top_builddir)/libs/librevseqcmp.a \
	$(top_builddir)/libs/libseqcmp.a

if ENABLE_X86_SIMD

phylonium_LDADD+=\
	$(top_builddir)/libs/librevseqcmp_avx2.a \
	$(top_builddir)/libs/librevseqcmp_ssse3.a \
	$(top_builddir)/libs/libseqcmp_avx2.a \
	$(top_builddir)/libs/libseqcmp_sse2.a

if ENABLE_AVX512
phylonium_LDADD+=$(top_builddir)/libs/libseqcmp_avx512.a
endif
endif

.PHONY: format debug

debug: CXXFLAGS=-O0 -fno-omit-frame-pointer -ggdb
debug: phylonium

format:
	clang-format -i *.h *.cxx
