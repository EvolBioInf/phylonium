
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "esa.h"
#include "global.h"
#include "sequence.h"

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
size_t minAnchorLength(double p, double g, size_t l)
{
	size_t x = 1;

	double prop = 0.0;
	while (prop < 1 - p) {
		prop = shuprop(x, g / 2, l);
		x++;
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

auto find_non_matches(const esa &ref, double gc, sequence &query_seq)
{
	auto nm = std::vector<sequence>{};

	size_t border = ref.size() / 2;

	const char *query = query_seq.c_str();
	size_t query_length = query_seq.size();

	size_t last_pos_Q = 0;
	size_t last_pos_S = 0;
	size_t last_length = 0;
	// This variable indicates that the last anchor was the right anchor of a
	// pair.
	bool last_was_right_anchor = true;

	size_t this_pos_Q = 0;
	size_t this_pos_S;
	size_t this_length;

	size_t threshold =
		minAnchorLength(1 - std::sqrt(1 - RANDOM_ANCHOR_PROP), gc, ref.size());

	// Iterate over the complete query.
	while (this_pos_Q < query_length) {
		auto inter =
			ref.get_match_cached(query + this_pos_Q, query_length - this_pos_Q);

		this_length = inter.l <= 0 ? 0 : inter.l;
		if (this_pos_Q + this_length > query_length) {
			this_length = query_length - this_pos_Q;
		}

		if (inter.i == inter.j && this_length >= threshold) {
			// anchor
			this_pos_S = ref.SA[inter.i];

			if (this_pos_S > last_pos_S &&
				this_pos_Q - last_pos_Q == this_pos_S - last_pos_S &&
				(this_pos_S < border) == (last_pos_S < border)) {
				// right anchor

				last_was_right_anchor = true;
			} else {
				// no anchor pair, left anchor!
				if (last_was_right_anchor || last_length / 2 >= threshold) {
					// push gap
					auto last_end = last_pos_Q + last_length;
					if (this_pos_Q - last_end >= MIN_SPLIT_LENGTH) {
						auto subsequence =
							query_seq.sub(last_end, this_pos_Q - last_end);
						nm.push_back(subsequence);
					}
				}

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

	// no anchor?
	if (last_was_right_anchor && last_pos_Q == 0) {
		// push whole query
		nm.push_back(query_seq);
	}

	return nm;
}

std::vector<sequence> filter(sequence &reference, std::vector<sequence> &set)
{
	std::vector<sequence> non_matches;

	auto ref_esa = esa{reference};
	auto gc = gc_content(reference.get_nucl());

	// omp parallel for
	for (ssize_t i = 0; i < set.size(); ++i) {
		auto self = set[i];
		auto local_non_matches = find_non_matches(ref_esa, gc, self);

		{ // critical path
			auto inserter = std::back_inserter(non_matches);
			std::move(local_non_matches.begin(), local_non_matches.end(),
					  inserter);
		}
	}

	return non_matches;
}
