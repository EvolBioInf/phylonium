check_PROGRAMS= simf unittests

simf_SOURCES= simf.cxx
simf_CPPFLAGS= -std=c++14 -Wall -Wextra
simf_CXXFLAGS= -ggdb

unittests_SOURCES= \
	Tmain.cxx \
	Tsequence.cxx \
	$(top_srcdir)/src/sequence.cxx \
	$(top_srcdir)/src/sequence.h
unittests_CPPFLAGS= -I$(top_srcdir)/src -std=c++14 -Wall -Wextra
unittests_CXXFLAGS= -ggdb

dist_noinst_DATA= simple.sh

.PHONY: format

format:
	clang-format -i *.h *.cxx
