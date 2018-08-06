
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include "esa.h"
#include "evo_model.h"
#include "global.h"
#include "sequence.h"

char complement(char c);

double shuprop(size_t, double, size_t);

/** @brief Find an element based on its index.
 * @param first - Iterator tho the beginning of a range.
 * @param last - Iterator to the end of a range.
 * @param p - A predicate function return true or false depending on the index.
 * @returns an iterator to the first element for which the predicate returns
 * true, or last if none.
 */
template <class InputIt, class UnaryPredicate>
InputIt find_if_i(InputIt first, InputIt last, UnaryPredicate p)
{
	auto begin = first;
	for (; first != last; ++first) {
		if (p(first - begin)) {
			return first;
		}
	}
	return last;
}

/** @brief Remove elements from a range using their index.
 * @param first - Iterator to the beginning of a range.
 * @param last - Iterator to the end of a range.
 * @param p - A predicate function returning true for the indices to be removed.
 * @returns The border at the end of the remaining elements.
 */
template <class ForwardIt, class UnaryPredicate>
ForwardIt remove_if_i(ForwardIt first, ForwardIt last, UnaryPredicate p)
{
	auto begin = first;
	first = find_if_i(first, last, p);
	if (first != last)
		for (ForwardIt i = first; ++i != last;)
			if (!p(i - begin)) *first++ = std::move(*i);
	return first;
}

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

	while (shuprop(x, g / 2, l) < 1 - p) {
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
 * @param p - The the half of the relative amount of GC in the DNA.
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
	/** Store the direction of a match. */
	enum class dir { forward, reverse } direction = dir::forward;
	/** The beginning on the reference genome. */
	size_t index_reference = 0;
	/** The beginning on the forward strand of the reference genome. */
	size_t index_reference_projected = 0;
	/** The beginning on the query genome. */
	size_t index_query = 0;
	/** The length of the homologous region. */
	size_t length = 0;

	/** @brief Reasonable default constructor. */
	homology() = default;

	/** @brief Create a new homology from coordinates. */
	homology(size_t ir, size_t iq, size_t l = 0) noexcept
		: direction{dir::forward}, index_reference{ir},
		  index_reference_projected{ir}, index_query{iq}, length{l}
	{
	}

	/** @brief Extend the current homologous region to the right.
	 * @param stride - amount of elongation.
	 * @returns the new length.
	 */
	auto extend(size_t stride) noexcept
	{
		return length += stride;
	}

	/** @brief Transform coordinates, if necessary.
	 * If the homology is within the reverse region of the reference project the
	 * coordinates onto the forward strand.
	 * @param reference_length - The length of the reference.
	 */
	void reverseEh(size_t reference_length) noexcept
	{
		if (index_reference < reference_length) return;

		// hack, see esa::esa()
		index_reference_projected =
			2 * reference_length + 1 - length - index_reference;
		direction = dir::reverse;
	}

	/** @brief Check for overlap.
	 * @param other - another homology.
	 * @returns true iff this homology overlaps the argument.
	 */
	bool overlaps(const homology &other) const noexcept
	{
		if (index_reference_projected == other.index_reference_projected) {
			return true;
		}
		if (starts_left_of(other)) return !ends_left_of(other);
		if (other.starts_left_of(*this)) return !other.ends_left_of(*this);

		// no more case
		// should not be executed
		return false;
	}

	/** @brief Determine the order of two homologies.
	 * @param other - a homology to compare with.
	 * @returns true if this homology starts left of the argument wrt. the
	 * reference genome.
	 */
	bool starts_left_of(const homology &other) const noexcept
	{
		return index_reference_projected < other.index_reference_projected;
	}

	/** @brief Determine the order of two homologies.
	 * @param other - a homology to compare with.
	 * @returns true if this homology ends left of the argument wrt. the
	 * reference genome.
	 */
	bool ends_left_of(const homology &other) const noexcept
	{
		return index_reference_projected + length <=
			   other.index_reference_projected;
	}
};

/** @brief Compute the length of the common prefix.
 * Given two strings return the number of characters being equal at the
 * beginning of both string.
 * @param S - One string.
 * @param Q - Another string.
 * @param remaining - The maximum number of characters to investigate.
 * @returns the LCP.
 */
size_t lcp(const char *S, const char *Q, size_t remaining)
{
	size_t length = 0;
	while (length < remaining && S[length] == Q[length]) {
		length++;
	}
	return length;
}

/**
 * @brief Compute homologies between two sequences.
 *
 * A homology is defined as a region on two sequences which two flanking
 * anchors. These anchors are long exact matches. For a detailed explanation see
 * Haubold et al 2014.
 *
 * @param ref - the reference sequence and ESA.
 * @param gc - the GC content of the reference.
 * @param seq - the query.
 * @returns a list of homologies.
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

	auto anchor = [&]() {
		auto inter = ref.get_match_cached(seq.c_str() + this_pos_Q,
										  query_length - this_pos_Q);
		this_length = std::max(inter.l, 0);
		this_pos_S = ref.SA[inter.i];
		return inter.i == inter.j && this_length >= threshold;
	};

	auto lucky_anchor = [&]() {
		size_t advance = this_pos_Q - last_pos_Q;
		size_t gap = this_pos_Q - last_pos_Q - last_length;

		size_t try_pos_S = last_pos_S + advance;

		if (try_pos_S >= static_cast<size_t>(ref.size()) || gap > threshold) {
			return false;
		}

		this_pos_S = try_pos_S;
		this_length = lcp(&seq.c_str()[this_pos_Q], &ref.S.c_str()[try_pos_S],
						  query_length - this_pos_Q);

		return this_length >= threshold;
	};

	// Iterate over the complete query.
	while (this_pos_Q < query_length) {
		if (lucky_anchor() || anchor()) {
			// anchor

			size_t end_S = last_pos_S + last_length;
			size_t end_Q = last_pos_Q + last_length;
			if (this_pos_S > end_S &&
				this_pos_Q - end_Q == this_pos_S - end_S &&
				(this_pos_S < border) == (last_pos_S < border)) {
				// right anchor

				current.extend(this_pos_Q - end_Q + this_length);

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
evo_model compare(const sequence &sa, const std::vector<homology> &ha,
				  const sequence &sb, const std::vector<homology> &hb);

/** @brief Remove all homologies that overlap any other.
 *
 * The result of `anchor_homologies` is a set of homologies. However, some of
 * these may overlap on the reference interfering with mutation counting later
 * on. Thus matches have to be filtered.
 *
 * This function removes all homologies which overlap with any other one.
 *
 * @param pile - in out parameter with homologies.
 */
void filter_overlaps_strict(std::vector<homology> &pile)
{
	if (pile.size() < 2) return;

	auto next = std::begin(pile) + 1;
	size_t border = 0;

	// capture mode has to be by-reference, otherwise bad things happen.
	// See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=81482
	auto filter = [&](const homology &homo) {
		bool overlaps_left = border > homo.index_reference_projected;
		border = std::max(border, homo.index_reference_projected + homo.length);
		bool overlaps_right = homo.overlaps(*next);
		next++;
		return overlaps_left || overlaps_right;
	};

	// the last homology can only overlap to the left
	auto remove_last = (pile.end() - 2)->overlaps(*(pile.end() - 1));
	auto split = std::remove_if(pile.begin(), pile.end() - 1, filter);

	if (!remove_last) {
		std::swap(*split, *(pile.end() - 1));
		split++;
	}

	pile.erase(split, pile.end());
}

/** @brief Maximizes the number of homologous nucleotides by carefully removing
 * all but one of the overlapping homologies.
 *
 * The result of `anchor_homologies` is a set of homologies. However, some of
 * these may overlap on the reference interfering with mutation counting later
 * on. Thus matches have to be filtered.
 *
 * This function selectively reduces overlapping stacks to only a single
 * homology. It does this by chaining neighboring homologies. The chain with
 * most nucleotides is then kept, all homologies not in it get removed.
 *
 * @param pile - in out parameter with homologies.
 */
void filter_overlaps_max(std::vector<homology> &pile)
{
	if (pile.size() < 2) return;

	size_t size = pile.size();

	auto predecessor = std::vector<ssize_t>(size, -1);
	auto score = std::vector<size_t>(size, 0);

	assert(predecessor.size() == size);
	assert(score.size() == size);

	score[0] = pile[0].length;

	for (ssize_t i = 1; i < size; i++) {
		// find maximum, so far
		auto max_value = (ssize_t)0;
		auto max_index = (ssize_t)0;

		for (ssize_t k = 0; k < i; k++) {
			if (!pile[k].ends_left_of(pile[i])) continue;

			if (score[k] > max_value) {
				max_value = score[k];
				max_index = k;
			}
		}

		assert(0 <= max_index && max_index < i);
		predecessor[i] = max_index;
		score[i] = score[max_index] + pile[i].length;
	}

	auto visited = std::vector<bool>(size, false);

	auto that = std::max_element(score.begin(), score.end());

	ssize_t index = that - score.begin();
	while (index >= 0) {
		visited[index] = true;
		index = predecessor[index];
	}

	auto split = remove_if_i(pile.begin(), pile.end(),
							 [&](size_t index) { return !visited[index]; });

	pile.erase(split, pile.end());
}

/** @brief Process the input and do a lot of stuff.
 * @param subject - The subject sequence.
 * @param queries - All sequences.
 * @returns a "matrix" of mutation counts.
 */
std::vector<evo_model> process(const sequence &subject,
							   const std::vector<sequence> &queries)
{

	auto N = queries.size();
	auto start_time = std::chrono::steady_clock::now();
	auto ref = esa(subject); // seq + # + reverse
	auto homologies = std::vector<std::vector<homology>>(N);

	auto gc = gc_content(subject.get_nucl());

	if (FLAGS & flags::verbose) {
		std::cerr << "ref: " << subject.get_name() << std::endl;
		std::cerr << "length: " << subject.get_nucl().size() << std::endl;
	}

// now compare every sequence to the subject
#pragma omp parallel for num_threads(THREADS)
	for (size_t j = 0; j < N; j++) {
		// std::cerr << "query: " << j << std::endl;
		auto query = queries[j];
		auto hvlocal = anchor_homologies(ref, gc, query);
		std::sort(begin(hvlocal), end(hvlocal),
				  [](const homology &self, const homology &other) {
					  return self.starts_left_of(other);
				  });

		filter_overlaps_max(hvlocal);

#pragma omp critical
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
	auto M = [&matrix, N = N](size_t i, size_t j) -> evo_model & {
		return matrix[i * N + j];
	};

#pragma omp parallel for num_threads(THREADS)
	for (size_t i = 0; i < N; i++) {
		for (size_t j = i + 1; j < N; j++) {

			M(j, i) = M(i, j) +=
				compare(queries[i], homologies[i], queries[j], homologies[j]);
		}
	}

	return matrix;
}

char complement(char c)
{
	if (evo_model::hash(c) < 0) return 0;
	static const char inv[] = "TGCA";
	return inv[evo_model::hash(c)];
}

/** @brief Compare two sequences based on their precomputed homologies wrt. the
 * reference.
 * @param sa - One sequence.
 * @param ha - The homologies of sequence A with the reference.
 * @param sb - Another sequence.
 * @param hb - The homologies of sequence B with the reference.
 * @returns mutation counts.
 */
evo_model compare(const sequence &sa, const std::vector<homology> &ha,
				  const sequence &sb, const std::vector<homology> &hb)
{
	auto mutations = evo_model();

	const auto &other = hb;
	auto right_ptr = begin(other);

	auto pile = std::vector<homology>();

	for (const auto &homo : ha) {
		auto ends_left_of_homo = [&homo](const auto &other_homo) {
			return other_homo.ends_left_of(homo);
		};

		auto overlaps_homo = [&homo](const auto &other_homo) {
			return other_homo.overlaps(homo);
		};

		// erase homologies which are done
		auto split =
			std::remove_if(pile.begin(), pile.end(), ends_left_of_homo);
		pile.erase(split, end(pile));

		// skip elements left of homo!
		right_ptr = std::find_if_not(right_ptr, other.end(), ends_left_of_homo);

		// add new homologies
		auto far_right_ptr = find_if_not(right_ptr, other.end(), overlaps_homo);

		/* Note that this cannot be merged with the find above
		 * into a single copy_if. The list of homologies is sorted!
		 * Thus this combination is linear, whereas copy_if would
		 * have a O(n^2) complexity.
		 */
		std::copy(right_ptr, far_right_ptr, std::back_inserter(pile));
		right_ptr = far_right_ptr;

		// compare homo against pile
		for (const auto &other_homo : pile) {
			mutations += compare(sa, homo, sb, other_homo);
		}
	}

	return mutations;
}

/** @brief Compare the homologous region of two sequences.
 * @param sa - One sequence.
 * @param ha - The homology of sequence A with the reference.
 * @param sb - Another sequence.
 * @param hb - The homology of sequence B with the reference.
 * @returns mutation counts.
 */
evo_model compare(const sequence &sa, const homology &ha, const sequence &sb,
				  const homology &hb)
{
	if (!ha.overlaps(hb)) {
		return evo_model{};
	}

	auto count = evo_model{};
	size_t common_start =
		std::max(ha.index_reference_projected, hb.index_reference_projected);

	size_t common_end = std::min(ha.index_reference_projected + ha.length,
								 hb.index_reference_projected + hb.length);

	assert(common_start < common_end);
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
		count.account(sa.c_str() + a_offset - length,
					  sb.c_str() + b_offset - length, length);
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
