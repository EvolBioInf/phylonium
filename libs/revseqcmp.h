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

int is_complement(char c, char d);

size_t revseqcmp(const char *begin, const char *other, size_t length);

#ifdef ENABLE_X86_SIMD
size_t revseqcmp_ssse3(const char *begin, const char *other, size_t length);
#endif

typedef size_t(revseqcmp_fn)(const char *begin, const char *other,
							 size_t length);

#ifdef __cplusplus
}
#endif
