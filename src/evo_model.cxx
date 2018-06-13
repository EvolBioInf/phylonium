#include "evo_model.h"
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iterator>
#include <numeric>

int evo_model::hash(char nucl) noexcept
{
	if (nucl < 'A') return -1;
	nucl &= 6;
	nucl ^= nucl >> 1;
	return nucl >> 1;
}

void evo_model::account(char a, char b) noexcept
{
	auto hash_a = hash(a);
	auto hash_b = hash(b);
	if (hash_a < 0 || hash_b < 0) return;
	if (b < a) {
		std::swap(a, b);
		std::swap(hash_a, hash_b);
	}

	// ensure a <= b

	int index;
	switch (hash_a) {
		case 0: index = AtoA; break;
		case 1: index = CtoC; break;
		case 2: index = GtoG; break;
		case 3: index = TtoT; break;
	}

	index += hash_b - hash_a;
	counts[index]++;
}

void evo_model::account(const char *sa, const char *sb, size_t length) noexcept
{
	size_t mutations = 0;
	for (size_t k = 0; k < length; k++) {
		if (sa[k] != sb[k]) {
			mutations++;
		}
	}

	counts[AtoA] += length - mutations;
	counts[AtoC] += mutations;
}

constexpr bool is_complement(char c, char d)
{
	if (c == 'A') return d == 'T';
	if (c == 'C') return d == 'G';
	if (c == 'G') return d == 'C';
	if (c == 'T') return d == 'A';
	return false;
	// return (~(c ^ d) & 2) && (c != d);
}

void evo_model::account_rev(const char *sa, const char *sb, size_t b_offset,
							size_t length) noexcept
{
	size_t mutations = 0;
	for (size_t k = 0; k < length; k++) {
		if (!is_complement(sa[k], sb[b_offset - 1 - k])) {
			mutations++;
		}
	}

	counts[AtoA] += length - mutations;
	counts[AtoC] += mutations;
}

evo_model &evo_model::operator+=(const evo_model &other) noexcept
{
	for (int i = 0; i < MUTCOUNTS; i++) {
		counts[i] += other.counts[i];
	}

	for (int i = MUTCOUNTS; i < sizeof(counts) / sizeof(counts[0]); i++)
		assert(counts[i] == 0);

	return *this;
}

size_t evo_model::total() const noexcept
{
	return std::accumulate(std::begin(counts), std::end(counts), 0);
}

double evo_model::estimate_raw() const noexcept
{
	size_t nucl = total();
	if (nucl == 0) return 0;
	size_t SNPs = 0;

	SNPs += counts[AtoC];
	SNPs += counts[AtoG];
	SNPs += counts[AtoT];
	SNPs += counts[CtoG];
	SNPs += counts[CtoT];
	SNPs += counts[GtoT];

	return SNPs / (double)nucl;
}

double evo_model::estimate_JC() const noexcept
{
	auto dist = estimate_raw();
	dist = -0.75 * log(1.0 - (4.0 / 3.0) * dist); // jukes cantor

	// fix negative zero
	return dist <= 0.0 ? 0.0 : dist;
}
