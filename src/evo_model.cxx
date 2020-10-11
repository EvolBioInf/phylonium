/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 - 2020 © Fabian Klötzl
 */
#include "evo_model.h"
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include "revseqcmp.h"
#include "seqcmp.h"

extern std::mt19937 prng;

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
	size_t mutations = seqcmp(sa, sb, length);

	homologs += length;
	substitutions += mutations;
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
	size_t mutations = revseqcmp(sa, sb + b_offset - length, length);

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

	std::binomial_distribution<> d(homologs, subst_rate);

	ret.substitutions = d(prng);

	return ret;
}

/** @brief Computes the 'coverage'.
 * @returns the number of nucleotides covered by this homology model.
 */
double evo_model::coverage(size_t length) const noexcept
{
	return (double)total() / length;
}
