noinst_LIBRARIES= libpfasta.a libcompat.a libseqcmp.a libseqcmp_sse2.a libseqcmp_avx2.a libseqcmp_avx512.a

libpfasta_a_SOURCES= pfasta.c pfasta.h
libpfasta_a_CPPFLAGS= -I$(top_srcdir)/opt

libcompat_a_SOURCES = compat-stdlib.h
if !HAVE_REALLOCARRAY
libcompat_a_SOURCES+= reallocarray.c
endif

libseqcmp_a_SOURCES = seqcmp.c
libseqcmp_a_LIBADD = libseqcmp_sse2.a libseqcmp_avx2.a libseqcmp_avx512.a
libseqcmp_sse2_a_SOURCES = seqcmp_sse2.c
libseqcmp_sse2_a_CFLAGS = -mpopcnt -msse2
libseqcmp_avx2_a_SOURCES = seqcmp_avx2.c
libseqcmp_avx2_a_CFLAGS = -mpopcnt -mavx2
libseqcmp_avx512_a_SOURCES = seqcmp_avx512.c
libseqcmp_avx512_a_CFLAGS = -mpopcnt -mavx512bw -mavx512vl

format:
	clang-format -i *.h *.c

