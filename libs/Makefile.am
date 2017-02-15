noinst_LIBRARIES= libpfasta.a libcompat.a

libpfasta_a_SOURCES= pfasta.c pfasta.h
libpfasta_a_CPPFLAGS= -I$(top_srcdir)/opt

libcompat_a_SOURCES = compat-stdlib.h
if !HAVE_REALLOCARRAY
libcompat_a_SOURCES+= reallocarray.c
endif
