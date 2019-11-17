noinst_LIBRARIES= \
	libcompat.a \
	libpfasta.a \
	librevseqcmp.a \
	libseqcmp.a

libpfasta_a_SOURCES= pfasta.c pfasta.h
libpfasta_a_CPPFLAGS= -I$(top_srcdir)/opt

libcompat_a_SOURCES = compat-stdlib.h
if !HAVE_REALLOCARRAY
libcompat_a_SOURCES+= reallocarray.c
endif

libseqcmp_a_SOURCES = seqcmp.c seqcmp.h

if ENABLE_X86_SIMD
noinst_LIBRARIES+= libseqcmp_sse2.a libseqcmp_avx2.a
libseqcmp_a_LIBADD = libseqcmp_sse2.a libseqcmp_avx2.a
libseqcmp_sse2_a_SOURCES = seqcmp_sse2.c
libseqcmp_sse2_a_CFLAGS = -mpopcnt -msse2
libseqcmp_avx2_a_SOURCES = seqcmp_avx2.c
libseqcmp_avx2_a_CFLAGS = -mpopcnt -mavx2

if ENABLE_AVX512
noinst_LIBRARIES+=libseqcmp_avx512.a
libseqcmp_a_LIBADD+=libseqcmp_avx512.a
libseqcmp_avx512_a_SOURCES = seqcmp_avx512.c
libseqcmp_avx512_a_CFLAGS = -mpopcnt -mavx512bw -mavx512vl
endif
endif

librevseqcmp_a_SOURCES = revseqcmp.c revseqcmp.h
librevseqcmp_a_LIBADD =

if ENABLE_X86_SIMD
noinst_LIBRARIES+= librevseqcmp_ssse3.a librevseqcmp_avx2.a
librevseqcmp_a_LIBADD+= librevseqcmp_ssse3.a librevseqcmp_avx2.a
librevseqcmp_ssse3_a_SOURCES = revseqcmp_ssse3.c
librevseqcmp_ssse3_a_CFLAGS = -mpopcnt -mssse3
librevseqcmp_avx2_a_SOURCES = revseqcmp_avx2.c
librevseqcmp_avx2_a_CFLAGS = -mpopcnt -mavx2
endif

format:
	clang-format -i *.h *.c
