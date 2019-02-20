/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 (C) Fabian Klötzl
 */

#include <assert.h>
#include <immintrin.h>
#include <string.h>
#include "seqcmp.h"

typedef __m256i vec_type;

size_t seqcmp_avx2(const char *begin, const char *other, size_t length)
{
	assert(begin != NULL);
	assert(other != NULL);

	size_t substitutions = 0;
	size_t offset = 0;

	const size_t vec_bytes = sizeof(vec_type);
	size_t vec_offset = 0;
	size_t vec_length = (length / vec_bytes) & ~(size_t)1; // round down

	size_t equal = 0;
	for (; vec_offset < vec_length; vec_offset++) {
		vec_type begin_chunk;
		vec_type other_chunk;
		memcpy(&begin_chunk, begin + vec_offset * vec_bytes, vec_bytes);
		memcpy(&other_chunk, other + vec_offset * vec_bytes, vec_bytes);

		vec_type comp = _mm256_cmpeq_epi8(begin_chunk, other_chunk);

		unsigned int vmask = _mm256_movemask_epi8(comp);
		equal += __builtin_popcount(vmask);

		vec_offset++;
		// second pass
		memcpy(&begin_chunk, begin + vec_offset * vec_bytes, vec_bytes);
		memcpy(&other_chunk, other + vec_offset * vec_bytes, vec_bytes);

		comp = _mm256_cmpeq_epi8(begin_chunk, other_chunk);

		vmask = _mm256_movemask_epi8(comp);
		equal += __builtin_popcount(vmask);
	}

	substitutions = vec_offset * vec_bytes - equal;

	offset += vec_offset * vec_bytes;

	for (; offset < length; offset++) {
		if (begin[offset] != other[offset]) {
			substitutions++;
		}
	}

	return substitutions;
}
