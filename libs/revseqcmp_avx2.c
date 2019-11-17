/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2019 © Fabian Klötzl
 */
#include <emmintrin.h>
#include <immintrin.h>
#include <stddef.h>
#include <string.h>
#include <tmmintrin.h>
#include "revseqcmp.h"

size_t revseqcmp_avx2(const char *self, const char *other, size_t length)
{
	size_t substitutions = 0;
	size_t offset = 0;

	typedef __m256i vec_type;
	size_t vec_size = sizeof(vec_type);
	size_t vec_offset = 0;
	size_t vec_length = length / vec_size;

	substitutions += vec_size * vec_length;

	vec_type mask =
		_mm256_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0,
						1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
	vec_type mask6 = _mm256_set1_epi8(6);
	vec_type mask4 = _mm256_set1_epi8(4);

	for (; vec_offset < vec_length; vec_offset++) {
		vec_type b;
		memcpy(&b, self + vec_offset * vec_size, vec_size);
		vec_type o;
		size_t pos = length - (vec_offset + 1) * vec_size;
		memcpy(&o, other + pos, vec_size);

		vec_type reversed = _mm256_shuffle_epi8(o, mask);
		vec_type swapped = _mm256_permute2x128_si256(b, b, 1);

		vec_type v1 = _mm256_xor_si256(swapped, reversed);
		vec_type v2 = _mm256_and_si256(v1, mask6);
		vec_type v3 = _mm256_cmpeq_epi8(v2, mask4);

		unsigned int vmask = _mm256_movemask_epi8(v3);
		substitutions -= __builtin_popcount(vmask);
	}

	offset += vec_offset * vec_size;

	for (; offset < length; offset++) {
		if (!is_complement(self[offset], other[length - 1 - offset])) {
			substitutions++;
		}
	}

	return substitutions;
}
