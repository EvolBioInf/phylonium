/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2019 © Fabian Klötzl
 */
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include "config.h"

size_t seqcmp(const char *begin, const char *other, size_t length);

#ifdef ENABLE_X86_SIMD
size_t seqcmp_sse2(const char *begin, const char *other, size_t length);
size_t seqcmp_avx2(const char *begin, const char *other, size_t length);

#ifdef ENABLE_AVX512
size_t seqcmp_avx512(const char *begin, const char *other, size_t length);
#endif
#endif

typedef size_t(seqcmp_fn)(const char *begin, const char *other, size_t length);

#ifdef __cplusplus
}
#endif
