check_PROGRAMS= simf
simf_SOURCES= simf.cxx
simf_CPPFLAGS= -std=c++14 -Wall -Wextra
simf_CXXFLAGS= -ggdb

dist_noinst_DATA= simple.sh

.PHONY: format

format:
	clang-format -i *.h *.cxx
