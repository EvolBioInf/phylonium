
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "esa.h"
#include "global.h"
#include "sequence.h"

char complement(char c);

class evo_model
{
  protected:
	enum {
		AtoA,
		AtoC,
		AtoG,
		AtoT,
		CtoC,
		CtoG,
		CtoT,
		GtoG,
		GtoT,
		TtoT,
		MUTCOUNTS,
		CtoA = AtoC,
		GtoA = AtoG,
		GtoC = CtoG,
		TtoA = AtoT,
		TtoC = CtoT,
		TtoG = GtoT
	};
	int counts[MUTCOUNTS] = {0};

  public:
	evo_model() = default;

	static int hash(char nucl)
	{
		if (nucl < 'A') return -1;
		nucl &= 6;
		nucl ^= nucl >> 1;
		return nucl >> 1;
	}

	void account(char a, char b) noexcept
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

	void account(const char *sa, const char *sb, size_t length) noexcept
	{
		size_t mutations = 0;
		for (size_t k = 0; k < length; k++) {
			if (sa[k] != sb[k]) {
				mutations++;
			}
		}

		counts[AtoA] = length - mutations;
		counts[AtoC] = mutations;
	}

	constexpr bool is_complement(char c, char d) const noexcept
	{
		return (~(c ^ d) & 2) && (c != d);
	}

	void account_rev(const char *sa, const char *sb, size_t b_offset, size_t length) noexcept
	{
		size_t mutations = 0;
		for (size_t k = 0; k < length; k++) {
			if (!is_complement(sa[k], sb[b_offset - 1 - k])) {
				mutations++;
			}
		}

		counts[AtoA] = length - mutations;
		counts[AtoC] = mutations;
	}

	evo_model &operator+=(const evo_model &other) noexcept
	{
		for (int i = 0; i < MUTCOUNTS; i++) {
			counts[i] += other.counts[i];
		}

		return *this;
	}

	size_t total() const noexcept
	{
		return std::accumulate(std::begin(counts), std::end(counts), 0);
	}

	double estimate_raw() const noexcept
	{
		size_t nucl = total();
		if (nucl == 0) return NAN;
		size_t SNPs = 0;

		SNPs += counts[AtoC];
		SNPs += counts[AtoG];
		SNPs += counts[AtoT];
		SNPs += counts[CtoG];
		SNPs += counts[CtoT];
		SNPs += counts[GtoT];

		return SNPs / (double)nucl;
	}

	double estimate_JC() const noexcept
	{
		auto dist = estimate_raw();
		dist = -0.75 * log(1.0 - (4.0 / 3.0) * dist); // jukes cantor

		// fix negative zero
		return dist <= 0.0 ? 0.0 : dist;
	}
};

double shuprop(size_t, double, size_t);

/**
 * @brief Calculates the minimum anchor length.
 *
 * Given some parameters calculate the minimum length for anchors according
 * to the distribution from Haubold et al. (2009).
 *
 * @param p - The probability with which an anchor is allowed to be random.
 * @param g - The the relative amount of GC in the subject.
 * @param l - The length of the subject.
 * @returns The minimum length of an anchor.
 */
size_t min_anchor_length(double p, double g, size_t l)
{
	size_t x = 1;

	double prop = 0.0;
	for (; prop < 1 - p; x++) {
		prop = shuprop(x, g / 2, l);
	}

	return x;
}

/**
 * @brief Calculates the binomial coefficient of n and k.
 *
 * We used to use gsl_sf_lnchoose(xx,kk) for this functionality.
 * After all, why implement something that has already been done?
 * Well, the reason is simplicity: GSL is used for only this one
 * function and the input (n<=20) is not even considered big.
 * Hence its much easier to have our own implementation and ditch
 * the GSL dependency even if that means our code is a tiny bit
 * less optimized and slower.
 *
 * @param n - The n part of the binomial coefficient.
 * @param k - analog.
 * @returns (n choose k)
 */
size_t binomial_coefficient(size_t n, size_t k)
{
	if (n <= 0 || k > n) {
		return 0;
	}

	if (k == 0 || k == n) {
		return 1;
	}

	if (k > n - k) {
		k = n - k;
	}

	size_t res = 1;

	for (size_t i = 1; i <= k; i++) {
		res *= n - k + i;
		res /= i;
	}

	return res;
}

/**
 * @brief Given `x` this function calculates the probability of a shustring
 * with a length less than `x`.
 *
 * Let X be the longest shortest unique substring (shustring) at any position.
 * Then this function computes P{X <= x} with respect to the given parameter
 * set. See Haubold et al. (2009).
 *
 * @param x - The maximum length of a shustring.
 * @param g - The the half of the relative amount of GC in the DNA.
 * @param l - The length of the subject.
 * @returns The probability of a certain shustring length.
 */
double shuprop(size_t x, double p, size_t l)
{
	double xx = (double)x;
	double ll = (double)l;
	size_t k;

	double s = 0.0;

	for (k = 0; k <= x; k++) {
		double kk = (double)k;
		double t = pow(p, kk) * pow(0.5 - p, xx - kk);

		s += pow(2, xx) * (t * pow(1 - t, ll)) *
			 (double)binomial_coefficient(x, k);
		if (s >= 1.0) {
			s = 1.0;
			break;
		}
	}

	return s;
}

class homology
{
  public:
	enum class dir { forward, reverse } direction = dir::forward;
	size_t index_reference = 0;
	size_t index_reference_projected = 0;
	size_t index_query = 0;
	size_t length = 0;

	homology() = default;
	homology(size_t ir, size_t iq, size_t l = 0) noexcept
		: direction{dir::forward},
		  index_reference{ir},
		  index_reference_projected{ir},
		  index_query{iq},
		  length{l}
	{
	}

	auto extend(size_t stride) noexcept
	{
		return length += stride;
	}

	/** @brief Transform coordinates, if necessary.
	 *
	 */
	void reverseEh(size_t reference_length) noexcept
	{
		if (index_reference < reference_length) return;

		// hack, see esa::esa()
		index_reference_projected =
			2 * reference_length + 1 - length - index_reference;
		direction = dir::reverse;
	}

	bool overlaps(const homology &other) const noexcept
	{
		if (index_reference_projected == other.index_reference_projected) {
			return true;
		}
		if (index_reference_projected < other.index_reference_projected) {
			return index_reference_projected + length >
				   other.index_reference_projected;
		}
		// else: index_reference_projected > other.index_reference_projected
		return index_reference_projected <
			   other.index_reference_projected + other.length;
	}

	bool starts_left_of(const homology &other) const noexcept
	{
		return index_reference_projected < other.index_reference_projected;
	}
};

/**
 * @brief Compute homologies between two sequences.
 */
auto anchor_homologies(const esa &ref, double gc, const sequence &seq)
{
	auto hv = std::vector<homology>();

	size_t border = ref.size() / 2;

	size_t query_length = seq.size();

	size_t last_pos_Q = 0;
	size_t last_pos_S = 0;
	size_t last_length = 0;
	// This variable indicates that the last anchor was the right anchor of a
	// pair.
	bool last_was_right_anchor = false;

	size_t this_pos_Q = 0;
	size_t this_pos_S;
	size_t this_length;

	size_t threshold = min_anchor_length(RANDOM_ANCHOR_PROP, gc, ref.size());

	auto current = homology(0, 0);

	// Iterate over the complete query.
	while (this_pos_Q < query_length) {
		auto inter = ref.get_match_cached(seq.c_str() + this_pos_Q,
										  query_length - this_pos_Q);

		this_length = std::max(inter.l, 0);

		if (inter.i == inter.j && this_length >= threshold) {
			// anchor
			this_pos_S = ref.SA[inter.i];

			if (this_pos_S > last_pos_S &&
				this_pos_Q - last_pos_Q == this_pos_S - last_pos_S &&
				(this_pos_S < border) == (last_pos_S < border)) {
				// right anchor

				current.extend(this_pos_Q - last_pos_Q - last_length +
							   this_length);

				last_was_right_anchor = true;
			} else {
				// no anchor pair, left anchor!
				if (last_was_right_anchor || last_length / 2 >= threshold) {
					// push last
					current.reverseEh(border);
					hv.push_back(std::move(current));
				}

				current = homology(this_pos_S, this_pos_Q, this_length);

				// this is not right anchor
				last_was_right_anchor = false;
			}

			// Cache values for later
			last_pos_Q = this_pos_Q;
			last_pos_S = this_pos_S;
			last_length = this_length;
		}

		// Advance
		this_pos_Q += this_length + 1;
	}

	// Very special case: The sequences are identical
	if (last_length >= query_length) {
		current = homology(last_pos_S, 0, query_length);
	}

	if (last_was_right_anchor || last_length / 2 >= threshold) {
		current.reverseEh(border);
		hv.push_back(std::move(current));
	}

	return hv;
}

evo_model compare(const sequence &sa, const homology &ha, const sequence &sb,
				  const homology &hb);

/**
 *
 */
void process(const genome &reference, const std::vector<genome> &genomes)
{
	auto N = genomes.size();
	auto start_time = std::chrono::steady_clock::now();
	auto subject = join(reference);
	auto ref = esa(subject); // seq + # + reverse
	auto homologies = std::vector<std::vector<homology>>(N);

	auto gc = gc_content(subject.get_nucl());

	if (FLAGS & flags::verbose) {
		std::cerr << "ref: " << reference.get_name() << std::endl;
		std::cerr << "length: " << subject.get_nucl().size() << std::endl;
	}

	auto queries = std::vector<sequence>();
	auto inserter = std::back_inserter(queries);
	queries.reserve(N);
	std::transform(begin(genomes), end(genomes), inserter, join);

// now compare every sequence to the subject
#pragma omp parallel for num_threads(THREADS)
	for (size_t j = 0; j < N; j++) {
		// std::cerr << "query: " << j << std::endl;
		auto query = queries[j];
		auto hvlocal = anchor_homologies(ref, gc, query);
		std::sort(begin(hvlocal), end(hvlocal),
				  [](const homology &self, const homology &other) {
					  return self.index_reference_projected <
							 other.index_reference_projected;
				  });

		homologies[j] = std::move(hvlocal);
	}

	if (FLAGS & flags::verbose) {
		auto end_time = std::chrono::steady_clock::now();
		std::cerr << "aligning took "
				  << std::chrono::duration_cast<std::chrono::milliseconds>(
						 end_time - start_time)
						 .count()
				  << "ms.\n";
	}

	auto matrix = std::vector<evo_model>(N * N);
	auto M = [&matrix, N = N ](size_t i, size_t j)->evo_model &
	{
		return matrix[i * N + j];
	};

#pragma omp parallel for num_threads(THREADS)
	for (size_t i = 0; i < N; i++) {
		for (size_t j = i + 1; j < N; j++) {
			auto mutations = evo_model();
			// compare i and j
			// compare homologies[i] and homologies[j]
			const auto &other = homologies[j];
			auto right_ptr = begin(other);

			auto pile = std::vector<homology>();

			for (const auto &homo : homologies[i]) {
				auto ends_left_of_homo = [&homo](const auto &other_homo) {
					return other_homo.starts_left_of(homo) &&
						   !other_homo.overlaps(homo);
				};

				auto overlaps_homo = [&homo](const auto &other_homo) {
					return other_homo.overlaps(homo);
				};

				// erase homologies which are done
				pile.erase(
					std::remove_if(begin(pile), end(pile), ends_left_of_homo),
					end(pile));

				// skip elements left of homo!
				right_ptr =
					std::find_if_not(right_ptr, end(other), ends_left_of_homo);

				// add new homologies
				auto far_right_ptr =
					find_if_not(right_ptr, end(other), overlaps_homo);

				/* Note that this cannot be merged with the find above
				 * into a single copy_if. The list of homologies is sorted!
				 * Thus this combination is linear, whereas copy_if would
				 * have a O(n^2) complexity.
				 */
				std::copy(right_ptr, far_right_ptr, std::back_inserter(pile));
				right_ptr = far_right_ptr;

				// compare homo against pile
				for (const auto &other : pile) {
					if (!homo.overlaps(other)) continue;
					mutations += compare(queries[i], homo, queries[j], other);
				}
			}

			M(j, i) = M(i, j) += mutations;
		}
	}

	auto dist_matrix = std::vector<double>(N * N, NAN);
	std::transform(begin(matrix), end(matrix), begin(dist_matrix),
				   [](const evo_model &em) { return em.estimate_JC(); });

	std::cout << N << std::endl;
	for (size_t i = 0; i < N; i++) {
		std::cout << genomes[i].get_name();
		for (size_t j = 0; j < N; j++) {
			std::cout << "  " << std::setw(8) << std::setprecision(4)
					  << (i == j ? 0.0 : dist_matrix[i * N + j]);
		}
		std::cout << std::endl;
	}

	if (FLAGS & flags::verbose) {
		std::cerr << "\nCoverages:\n";
		for (size_t i = 0; i < N; i++) {
			for (size_t j = 0; j < N; j++) {
				std::cerr << std::setw(8) << std::setprecision(4)
						  << matrix[i * N + j].total() /
								 ((double)genomes[j].get_length())
						  << "  ";
			}
			std::cerr << std::endl;
		}
	}
}

char complement(char c)
{
	if (evo_model::hash(c) < 0) return 0;
	static const char inv[] = "TGCA";
	return inv[evo_model::hash(c)];
}

evo_model compare(const sequence &sa, const homology &ha, const sequence &sb,
				  const homology &hb)
{
	auto count = evo_model{};
	size_t common_start =
		std::max(ha.index_reference_projected, hb.index_reference_projected);

	size_t common_end = std::min(ha.index_reference_projected + ha.length,
								 hb.index_reference_projected + hb.length);

	size_t length = common_end - common_start;

	if (ha.direction == hb.direction &&
		ha.direction == homology::dir::forward) {
		// simple comparison
		auto a_offset =
			common_start - ha.index_reference_projected + ha.index_query;
		auto b_offset =
			common_start - hb.index_reference_projected + hb.index_query;

		count.account(sa.c_str() + a_offset, sb.c_str() + b_offset, length);
	} else if (ha.direction == hb.direction &&
			   ha.direction == homology::dir::reverse) {
		auto a_offset = ha.index_query + ha.length -
						(common_start - ha.index_reference_projected);
		auto b_offset = hb.index_query + hb.length -
						(common_start - hb.index_reference_projected);

		// no need for double complement.
		for (size_t i = 0; i < length; i++) {
			count.account(sb.c_str()[b_offset - 1 - i],
						  sa.c_str()[a_offset - 1 - i]);
		}
	} else if (hb.direction == homology::dir::reverse) {
		auto a_offset =
			common_start - ha.index_reference_projected + ha.index_query;

		// reverse b
		auto b_offset = hb.index_query + hb.length -
						(common_start - hb.index_reference_projected);

		count.account_rev(sa.c_str() + a_offset, sb.c_str(), b_offset, length);
	} else if (ha.direction == homology::dir::reverse) {
		auto b_offset =
			common_start - hb.index_reference_projected + hb.index_query;

		// reverse a
		auto a_offset = ha.index_query + ha.length -
						(common_start - ha.index_reference_projected);

		count.account_rev(sb.c_str() + b_offset, sa.c_str(), a_offset, length);
	}

	return count;
}
