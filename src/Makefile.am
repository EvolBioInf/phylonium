bin_PROGRAMS = phylonium
phylonium_SOURCES = phylonium.cxx io.cxx io.h sequence.cxx sequence.h esa.cxx esa.h process.cxx process.h global.h
phylonium_CPPFLAGS = $(OPENMP_CFLAGS) -std=c++14 -Wall -Wextra -Weffc++ -I$(top_srcdir)/libs
phylonium_CXXFLAGS = $(OPENMP_CXXFLAGS) -O3 -ggdb
phylonium_LDADD = $(top_builddir)/libs/libpfasta.a $(top_builddir)/libs/libcompat.a

.PHONY: format debug

debug: CXXFLAGS=-O0 -fno-omit-frame-pointer -ggdb
debug: phylonium

format:
	clang-format -i *.h *.cxx
