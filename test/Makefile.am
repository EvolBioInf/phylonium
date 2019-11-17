check_PROGRAMS= simf unittests

simf_SOURCES= simf.cxx
simf_CPPFLAGS= -std=c++14 -Wall -Wextra
simf_CXXFLAGS= -ggdb

unittests_SOURCES= \
	$(top_srcdir)/src/esa.cxx \
	$(top_srcdir)/src/esa.h \
	$(top_srcdir)/src/evo_model.cxx \
	$(top_srcdir)/src/evo_model.h \
	$(top_srcdir)/src/process.cxx \
	$(top_srcdir)/src/process.h \
	$(top_srcdir)/src/sequence.cxx \
	$(top_srcdir)/src/sequence.h \
	catch.hpp \
	Tmain.cxx \
	Tprocess.cxx \
	Tsequence.cxx
unittests_CPPFLAGS= -I$(top_srcdir)/src -I$(top_srcdir)/libs -std=c++14 -Wall -Wextra
unittests_CXXFLAGS= -ggdb
unittests_LDADD = \
	$(top_builddir)/libs/libcompat.a \
	$(top_builddir)/libs/libpfasta.a \
	$(top_builddir)/libs/librevseqcmp.a \
	$(top_builddir)/libs/libseqcmp.a

if ENABLE_X86_SIMD

unittests_LDADD+=\
	$(top_builddir)/libs/librevseqcmp_avx2.a \
	$(top_builddir)/libs/librevseqcmp_ssse3.a \
	$(top_builddir)/libs/libseqcmp_avx2.a \
	$(top_builddir)/libs/libseqcmp_sse2.a

if ENABLE_AVX512
unittests_LDADD+=$(top_builddir)/libs/libseqcmp_avx512.a
endif
endif

dist_noinst_DATA= simple.sh

.PHONY: format

format:
	clang-format -i *.h *.cxx
