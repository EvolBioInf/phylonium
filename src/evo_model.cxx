/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 © Fabian Klötzl
 */
#include "evo_model.h"
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iterator>
#include <numeric>

constexpr bool is_complement(char c, char d);

#ifdef __AVX512BW__
#ifdef __AVX512VL__
#include <immintrin.h>

typedef __m256i vec_type;

size_t count_subst_avx512(const char *self, const char *other, size_t length)
{
	size_t substitutions = 0;

	const size_t vec_bytes = sizeof(vec_type);
	size_t vec_offset = 0;
	size_t vec_length = length / vec_bytes;

	size_t equal = 0;
	for (; vec_offset < vec_length; vec_offset++) {
		vec_type self_chunk;
		vec_type other_chunk;
		memcpy(&self_chunk, self + vec_offset * vec_bytes, vec_bytes);
		memcpy(&other_chunk, other + vec_offset * vec_bytes, vec_bytes);

		unsigned int vmask = _mm256_cmpeq_epi8_mask(self_chunk, other_chunk);
		equal += __builtin_popcount(vmask);
	}

	substitutions = vec_length * vec_bytes - equal;

	size_t offset = vec_offset * vec_bytes;
	for (; offset < length; offset++) {
		if (self[offset] != other[offset]) {
			substitutions++;
		}
	}

	return substitutions;
}

#endif
#endif

#ifdef __AVX2__
#include <immintrin.h>

typedef __m256i vec_type;

size_t count_subst_avx2(const char *self, const char *other, size_t length)
{
	size_t substitutions = 0;
	size_t offset = 0;

	const size_t vec_bytes = sizeof(vec_type);
	size_t vec_offset = 0;
	size_t vec_length = (length / vec_bytes) & ~(size_t)1; // round down

	size_t equal = 0;
	for (; vec_offset < vec_length; vec_offset++) {
		vec_type self_chunk;
		vec_type other_chunk;
		memcpy(&self_chunk, self + vec_offset * vec_bytes, vec_bytes);
		memcpy(&other_chunk, other + vec_offset * vec_bytes, vec_bytes);

		vec_type comp = _mm256_cmpeq_epi8(self_chunk, other_chunk);

		unsigned int vmask = _mm256_movemask_epi8(comp);
		equal += __builtin_popcount(vmask);

		vec_offset++;
		// second pass
		memcpy(&self_chunk, self + vec_offset * vec_bytes, vec_bytes);
		memcpy(&other_chunk, other + vec_offset * vec_bytes, vec_bytes);

		comp = _mm256_cmpeq_epi8(self_chunk, other_chunk);

		vmask = _mm256_movemask_epi8(comp);
		equal += __builtin_popcount(vmask);
	}

	substitutions = vec_offset * vec_bytes - equal;

	offset += vec_offset * vec_bytes;

	for (; offset < length; offset++) {
		if (self[offset] != other[offset]) {
			substitutions++;
		}
	}

	return substitutions;
}

#endif

#ifdef __SSE2__
#ifdef __SSSE3__
#include <emmintrin.h>
#include <tmmintrin.h>

size_t intr(const char *self, const char *other, size_t length)
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

#endif
#endif

gsl_rng *RNG;

/** @brief Turn a nucleotide into a 2bit representation.
 * This function uses a clever way to transform nucleotides (ACGT) into a 2bit
 * representation. Maybe even too clever.
 * @param nucl - The nucleotide.
 * @returns returns a 2bit representation.
 */
int evo_model::hash(char nucl) noexcept
{
	if (nucl < 'A') return -1;
	nucl &= 6;
	nucl ^= nucl >> 1;
	return nucl >> 1;
}

/** @brief Compare two characters and update the counts.
 * @param a - A nucleotide from one sequence.
 * @param b - The nucleotide from the other sequence.
 */
void evo_model::account(char a, char b) noexcept
{
	homologs++;
	if (a != b) {
		substitutions++;
	}
}

/** Tell the compiler a branch is unlikely to be taken. */
#define UNLIKELY(X) __builtin_expect(X, 0)

/** @brief Compare two sequences.
 * @param sa - A pointer to a nucleotide sequence.
 * @param sb - A pointer to another nucleotide sequence.
 * @param length - The length of the homologous section.
 */
void evo_model::account(const char *sa, const char *sb, size_t length) noexcept
{
	size_t mutations = 0;
	size_t offset = 0;

#ifdef __AVX512BW__
#ifdef __AVX512VL__

	mutations = count_subst_avx512(sa, sb, length);
	offset += length;

#endif
#endif

#ifdef __AVX2__

	mutations += count_subst_avx2(sa + offset, sb + offset, length - offset);
	offset += length;

#endif

	for (; offset < length; offset++) {
		if (UNLIKELY(sa[offset] != sb[offset])) {
			mutations++;
		}
	}

	homologs += length;
	substitutions += mutations;
}

/** @brief Check whether two characters are complementary.
 * @param c - One nucleotide.
 * @param d - A nucleotide from the other sequence.
 * @returns true iff the two nucleotides are complements.
 */
constexpr bool is_complement(char c, char d)
{
	auto xorr = c ^ d;
	return (xorr & 6) == 4;
}

/** @brief Compare one sequence with the reverse complement of another.
 * @param sa - The forward sequence.
 * @param sb - The sequence of which the reverse complement is of interest.
 * @param b_offset - The offset of the reverse complement (TODO: get rid of
 * this).
 * @param length - The length of the homologous region.
 */
void evo_model::account_rev(const char *sa, const char *sb, size_t b_offset,
							size_t length) noexcept
{
	size_t mutations = 0;
	size_t offset = 0;

#ifdef __SSE2__
	mutations += intr(sa, sb + b_offset - length, length);
	offset += length;
#endif

	for (; offset < length; offset++) {
		if (UNLIKELY(!is_complement(sa[offset], sb[b_offset - 1 - offset]))) {
			mutations++;
		}
	}

	homologs += length;
	substitutions += mutations;
}

/** @brief Integrate the other count into this one.
 * @param other - The other substitution counts.
 * @returns a reference to the updated counts.
 */
evo_model &evo_model::operator+=(const evo_model &other) noexcept
{
	homologs += other.homologs;
	substitutions += other.substitutions;

	return *this;
}

/** @brief Return the number of homologous nucleotides.
 * @return the number of homologous nucleotides.
 */
size_t evo_model::total() const noexcept
{
	return homologs;
}

/** @brief Estimate the substitution rate.
 * @returns the rate of substitutions.
 */
double evo_model::estimate_raw(bool zero_on_error) const noexcept
{
	size_t nucl = total();
	if (nucl == 0) return zero_on_error ? 0.0 : NAN;

	size_t SNPs = substitutions;
	return SNPs / (double)nucl;
}

/** @brief Estimate the evolutionary distance via the Jukes-Cantor correction.
 * @returns the evolutionary distance.
 */
double evo_model::estimate_JC(bool zero_on_error) const noexcept
{
	auto dist = estimate_raw(zero_on_error);
	dist = -0.75 * log(1.0 - (4.0 / 3.0) * dist); // jukes cantor

	// fix negative zero
	return dist <= 0.0 ? 0.0 : dist;
}

/** @brief Bootstrap a new distance using a method by Haubold & Klötzl (2016).
 * @returns a new instance of homologous nucleotides.
 */
evo_model evo_model::bootstrap() const
{
	auto ret = *this; // copy

	auto subst_rate = substitutions / (double)homologs;

	std::array<double, 2> p = {subst_rate, 1 - subst_rate};
	std::array<unsigned int, 2> neu = {};

	gsl_ran_multinomial(RNG, 2, homologs, p.data(),
						reinterpret_cast<unsigned int *>(neu.data()));

	ret.substitutions = neu[0];

	return ret;
}

/** @brief Computes the 'coverage'.
 * @returns the number of nucleotides covered by this homology model.
 */
double evo_model::coverage(size_t length) const noexcept
{
	return (double)total() / length;
}
