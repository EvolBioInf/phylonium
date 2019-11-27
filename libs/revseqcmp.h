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

/** @brief Check whether two characters are complementary.
 * @param c - One nucleotide.
 * @param d - A nucleotide from the other sequence.
 * @returns true iff the two nucleotides are complements.
 */
extern inline int is_complement(char c, char d)
{
	int xorr = c ^ d;
	return (xorr & 6) == 4;
}

size_t revseqcmp(const char *begin, const char *other, size_t length);

#ifdef ENABLE_X86_SIMD
size_t revseqcmp_ssse3(const char *begin, const char *other, size_t length);
size_t revseqcmp_avx2(const char *begin, const char *other, size_t length);
#endif

typedef size_t(revseqcmp_fn)(const char *begin, const char *other,
							 size_t length);

#ifdef __cplusplus
}
#endif
