/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2019 © Fabian Klötzl
 */
#include <emmintrin.h>
#include <stddef.h>
#include <string.h>
#include <tmmintrin.h>
#include "revseqcmp.h"

/** @brief Check whether two characters are complementary.
 * @param c - One nucleotide.
 * @param d - A nucleotide from the other sequence.
 * @returns true iff the two nucleotides are complements.
 */
int is_complement(char c, char d)
{
	int xorr = c ^ d;
	return (xorr & 6) == 4;
}

size_t revseqcmp_ssse3(const char *self, const char *other, size_t length)
{
	size_t substitutions = 0;
	size_t offset = 0;

	size_t vec_offset = 0;
	size_t vec_length = length / sizeof(__m128i);

	substitutions += sizeof(__m128i) * vec_length;
	for (; vec_offset < vec_length; vec_offset++) {
		__m128i b;
		memcpy(&b, self + vec_offset * sizeof(__m128i), sizeof(b));
		__m128i o;
		size_t pos = length - (vec_offset + 1) * sizeof(__m128i);
		memcpy(&o, other + pos, sizeof(o));

		__m128i mask =
			_mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
		__m128i reversed = _mm_shuffle_epi8(o, mask);

		__m128i v1 = _mm_xor_si128(b, reversed);
		__m128i mask6 = _mm_set1_epi8(6);
		__m128i v2 = _mm_and_si128(v1, mask6);
		__m128i mask4 = _mm_set1_epi8(4);
		__m128i v3 = _mm_cmpeq_epi8(v2, mask4);

		unsigned int vmask = _mm_movemask_epi8(v3);
		// substitutions += sizeof(__m128i) - __builtin_popcount(vmask);
		substitutions -= __builtin_popcount(vmask);
	}

	offset += vec_offset * sizeof(__m128i);

	for (; offset < length; offset++) {
		if (!is_complement(self[offset], other[length - 1 - offset])) {
			substitutions++;
		}
	}

	return substitutions;
}
