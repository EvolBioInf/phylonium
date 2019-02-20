
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

int is_complement(char c, char d);

size_t revseqcmp(const char *begin, const char *other, size_t length);
size_t revseqcmp_ssse3(const char *begin, const char *other, size_t length);

typedef size_t(revseqcmp_fn)(const char *begin, const char *other,
							 size_t length);

#ifdef __cplusplus
}
#endif
