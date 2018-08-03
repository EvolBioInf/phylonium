#include "evo_model.h"
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>

gsl_rng *RNG;

int evo_model::hash(char nucl) noexcept
{
	if (nucl < 'A') return -1;
	nucl &= 6;
	nucl ^= nucl >> 1;
	return nucl >> 1;
}

void evo_model::account(char a, char b) noexcept
{
	homologs++;
	if (a != b) {
		substitutions++;
	}
}

#define UNLIKELY(X) __builtin_expect(X, 0)

void evo_model::account(const char *sa, const char *sb, size_t length) noexcept
{
	size_t mutations = 0;
	for (size_t k = 0; k < length; k++) {
		if (UNLIKELY(sa[k] != sb[k])) {
			mutations++;
		}
	}

	homologs += length;
	substitutions += mutations;
}

constexpr bool is_complement(char c, char d)
{
	auto xorr = c ^ d;
	return (xorr & 6) == 4;
}

void evo_model::account_rev(const char *sa, const char *sb, size_t b_offset,
							size_t length) noexcept
{
	size_t mutations = 0;
	for (size_t k = 0; k < length; k++) {
		if (UNLIKELY(!is_complement(sa[k], sb[b_offset - 1 - k]))) {
			mutations++;
		}
	}

	homologs += length;
	substitutions += mutations;
}

evo_model &evo_model::operator+=(const evo_model &other) noexcept
{
	homologs += other.homologs;
	substitutions += other.substitutions;

	return *this;
}

size_t evo_model::total() const noexcept
{
	return homologs;
}

double evo_model::estimate_raw() const noexcept
{
	size_t nucl = total();
	if (nucl == 0) return 0;

	size_t SNPs = substitutions;
	return SNPs / (double)nucl;
}

double evo_model::estimate_JC() const noexcept
{
	auto dist = estimate_raw();
	dist = -0.75 * log(1.0 - (4.0 / 3.0) * dist); // jukes cantor

	// fix negative zero
	return dist <= 0.0 ? 0.0 : dist;
}

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
