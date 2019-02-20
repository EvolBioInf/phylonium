
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

size_t seqcmp(const char *begin, const char *other, size_t length);
size_t seqcmp_sse2(const char *begin, const char *other, size_t length);
size_t seqcmp_avx2(const char *begin, const char *other, size_t length);
size_t seqcmp_avx512(const char *begin, const char *other, size_t length);

typedef size_t(seqcmp_fn)(const char *begin, const char *other, size_t length);

#ifdef __cplusplus
}
#endif
