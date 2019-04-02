/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 © Fabian Klötzl
 */
/**
 * @file
 * @brief ESA functions
 *
 * This file contains various functions that operate on an enhanced suffix
 * array. The basic algorithms originate from the book of Ohlebusch
 * "Bioinformatics Algorithms" (2013). Most of these were heavily modified
 * for improved performance. One example is the lcp-cache.
 *
 * The ESA structure defined in esa.h contains a `cache` field. This cache is
 * used to quickly look up lcp-intervals. Consider the queries "AAGT" and
 * "AACG". In both cases the interval for "AA" has to be looked up in the
 * ESA. If we simply store the interval for "AA" in the cache, once and use it
 * for each query we are significantly faster (up to 7 times).
 */
#include <memory>
#include <string>

#include "esa.h"
#include "global.h"

/** @brief The prefix length up to which LCP-intervals are cached. */
const size_t CACHE_LENGTH = 4;

/** @brief Map a code to the character. */
constexpr char code2char(ssize_t code) noexcept
{
	switch (code & 0x3) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
	}
	return '\0';
}

/** @brief Map a character to a two bit code. */
constexpr ssize_t char2code(char c) noexcept
{
	ssize_t result = -1;
	switch (c) {
		case 'A': result = 0; break;
		case 'C': result = 1; break;
		case 'G': result = 2; break;
		case 'T': result = 3; break;
	}
	return result;
}

/** @brief Initializes an ESA.
 *
 * This function initializes an ESA with respect to the provided sequence.
 * @param seq - The sequence
 */
esa::esa(const sequence &seq) : m_size{seq.size() * 2 + 1}
{
	m_master = seq;
	S = m_master.get_nucl() + '#' + reverse(m_master.get_nucl());
	SA = std::make_unique<saidx_t[]>(m_size);
	divsufsort(reinterpret_cast<const unsigned char *>(S.c_str()), SA.get(),
			   m_size);

	init_LCP();
	init_CLD();
	init_FVC();
	init_cache();
}

/** @brief Fills the LCP-Interval cache.
 *
 * Traversing the virtual suffix tree, created by SA, LCP and CLD is rather
 * slow. Hence we create a cache, holding the LCP-interval for a prefix of a
 * certain length ::CACHE_LENGTH. This function it the entry point for the
 * cache filling routine.
 */
void esa::init_cache()
{
	this->cache = std::make_unique<lcp_interval[]>(1 << (2 * CACHE_LENGTH));

	char str[CACHE_LENGTH + 1];
	str[CACHE_LENGTH] = '\0';

	saidx_t m = left_child(m_size);
	lcp_interval ij = {.l = LCP[m], .i = 0, .j = m_size - 1, .m = m};

	init_cache_dfs(str, 0, ij);
}

/** @brief Fills the cache — one char at a time.
 *
 * This function is a depth first search on the virtual suffix tree and fills
 * the cache. Or rather it calls it self until some value to cache is
 * calculated. This function is a recursive version of get_interval but with
 * more edge cases.
 *
 * @param str - The current prefix.
 * @param pos - The length of the prefix.
 * @param in - The LCP-interval of prefix[0..pos-1].
 */
void esa::init_cache_dfs(char *str, size_t pos, lcp_interval in)
{
	// we are not yet done, but the current strings do not exist in the subject.
	if (pos < CACHE_LENGTH && in.i == -1 && in.j == -1) {
		init_cache_fill(str, pos, in);
		return;
	}

	// we are past the caching length
	if (pos >= CACHE_LENGTH) {
		init_cache_fill(str, pos, in);
		return;
	}

	lcp_interval ij;

	// iterate over all nucleotides
	for (int code = 0; code < 4; ++code) {
		str[pos] = code2char(code);
		ij = get_interval(in, str[pos]);

		// fail early
		if (ij.i == -1 && ij.j == -1) {
			// if the current extension cannot be found, will with previous one
			init_cache_fill(str, pos + 1, in);
			continue;
		}

		// singleton
		if (ij.i == ij.j) {
			// fix length
			ij.l = pos + 1;
			init_cache_fill(str, pos + 1, ij);
			continue;
		}

		if (ij.l <= (ssize_t)(pos + 1)) {
			// Continue one level deeper
			// This is the usual case
			init_cache_dfs(str, pos + 1, ij);
			continue;
		}

		// The LCP-interval is deeper than expected
		// Check if it still fits into the cache
		if ((size_t)ij.l >= CACHE_LENGTH) {
			// If the lcp-interval exceeds the cache depth, stop here and fill
			init_cache_fill(str, pos + 1, in);
			continue;
		}

		/* At this point the prefix `str` of length `pos` has been found.
		 * However, the call to `getInterval` above found an interval with
		 * an LCP value bigger than `pos`. This means that not all elongations
		 * (more precise: just one) of `str` appear in the subject. Thus fill
		 * all values with the matched result to far and continue only with
		 * the one special substring.
		 */
		init_cache_fill(str, pos + 1, in);

		char non_acgt = 0;

		// fast forward
		size_t k = pos + 1;
		for (; k < (size_t)ij.l; k++) {
			// In some very edgy edge cases the lcp-interval `ij`
			// contains a `;` or another non-acgt character. Since we
			// cannot cache those, break.
			char c = S[SA[ij.i] + k];
			if (char2code(c) < 0) {
				non_acgt = 1;
				break;
			}

			str[k] = c;
		}

		// We are skipping intervals here. Maybe for each of them we should also
		// fill the cache. However, I haven't yet figured out how to do that
		// properly and whether it is worth it.

		if (non_acgt) {
			init_cache_fill(str, k, ij);
		} else {
			init_cache_dfs(str, k, ij);
		}
	}
}

/** @brief Fills the cache with a given value.
 *
 * Given a prefix and a value this function fills the cache beyond this point
 * the value.
 *
 * @param str - The current prefix.
 * @param pos - The length of the prefix.
 * @param in - The LCP-interval of prefix[0..pos-1].
 */
void esa::init_cache_fill(char *str, size_t pos, lcp_interval in)
{
	if (pos < CACHE_LENGTH) {
		for (int code = 0; code < 4; ++code) {
			str[pos] = code2char(code);
			init_cache_fill(str, pos + 1, in);
		}
	} else {
		ssize_t code = 0;
		for (size_t i = 0; i < CACHE_LENGTH; ++i) {
			code <<= 2;
			code |= char2code(str[i]);
		}

		cache[code] = in;
	}
}

/**
 * @brief Initializes the FVC (first variant character) array.
 *
 * The FVC is of my own invention and simply defined as
 * `FVC[i] = S[SA[i]+LCP[i]]`. This expression is constantly used in
 * get_interval. By precomputing the result, we have less memory
 * accesses, less cache misses, and thus improved runtimes of up to 15%
 * faster matching. This comes at a negligible cost of increased memory.
 */
void esa::init_FVC()
{
	FVC = std::make_unique<char[]>(m_size);
	FVC[0] = '\0';

	auto FVC_p = FVC.get();
	auto LCP_p = LCP.get();
	auto SA_p = SA.get();
	for (size_t i = m_size; i; i--, FVC_p++, SA_p++, LCP_p++) {
		*FVC_p = S[*SA_p + *LCP_p];
	}
}

/** @brief Initializes the CLD (child) array.
 *
 * See Ohlebusch.
 */
void esa::init_CLD()
{
	CLD = std::make_unique<saidx_t[]>(m_size + 1);

	typedef struct pair_s {
		saidx_t idx, lcp;
	} pair_t;

	auto stack = std::make_unique<pair_t[]>(m_size + 1);
	pair_t *top = stack.get(); // points at the topmost filled element
	pair_t last;

	right_child(0) = m_size;

	top->idx = 0;
	top->lcp = -1;

	// iterate over all elements
	for (ssize_t k = 1; k < m_size + 1; k++) {
		while (LCP[k] < top->lcp) {
			// top->lcp is a leaf
			last = *top--;

			// link all elements of same lcp value in a chain
			while (top->lcp == last.lcp) {
				right_child(top->idx) = last.idx;
				last = *top--;
			}

			// store the l-index of last
			if (LCP[k] < top->lcp) {
				right_child(top->idx) = last.idx;
			} else {
				left_child(k) = last.idx;
			}
		}

		// continue one level deeper
		top++;
		top->idx = k;
		top->lcp = LCP[k];
	}
}

/**
 * This function computed the LCP array, given the suffix array. Thereto it uses
 * a special `phi` array, which makes it slightly faster than the original
 * linear-time algorithm by Kasai et al.
 */
void esa::init_LCP()
{
	saidx_t len = m_size;

	// Allocate new memory
	// The LCP array is one element longer than S.
	LCP = std::make_unique<saidx_t[]>(len + 1);

	LCP[0] = -1;
	LCP[len] = -1;

	// Allocate temporary arrays
	auto PHI = std::make_unique<saidx_t[]>(len);
	auto &PLCP = PHI; //?

	PHI[SA[0]] = -1;
	saidx_t k;
	ssize_t i;

	for (i = 1; i < len; i++) {
		PHI[SA[i]] = SA[i - 1];
	}

	ssize_t l = 0;
	for (i = 0; i < len; i++) {
		k = PHI[i];
		if (k != -1) {
			while (S[k + l] == S[i + l]) {
				l++;
			}
			PLCP[i] = l;
			l--;
			if (l < 0) l = 0;
		} else {
			PLCP[i] = -1;
		}
	}

	// unpermutate the LCP array
	for (i = 1; i < len; i++) {
		LCP[i] = PLCP[SA[i]];
	}
}

/** @brief For the lcp-interval of string `w` compute the interval for `wa`
 *
 * Say, we already know the LCP-interval ij for a string `w`. Now we want to
 * check if `wa` may also be found in the ESA and thus in the subject. So we
 * look for the sub interval of `ij` in which all strings feature an `a` as
 * the next character. If such a sub interval is found, its boundaries are
 * returned.
 *
 * @param ij - The lcp-interval for `w`.
 * @param a - The next character.
 * @returns The lcp-interval one level deeper.
 */
lcp_interval esa::get_interval(lcp_interval ij, char a) const
{
	saidx_t i = ij.i;
	saidx_t j = ij.j;

	// check for singleton or empty interval
	if (i == j) {
		if (S[SA[i] + ij.l] != a) {
			ij.i = ij.j = -1;
		}
		return ij;
	}

	int m = ij.m;
	int l = ij.l;

	char c = S[SA[i] + l];
	goto SoSueMe;

	do {
		c = FVC[i];

	SoSueMe:
		if (c == a) {
			/* found ! */

			if (i != m - 1) {
				// found interval contains >1 element
				saidx_t n = left_child(m);

				ij = (lcp_interval){.l = LCP[n], .i = i, .j = m - 1, .m = n};
			} else {
				// empty or singleton
				// doing left_child(m) is not valid in this case!
				ij = (lcp_interval){.l = LCP[i], .i = i, .j = i, .m = -1};
			}

			return ij;
		}

		if (c > a) {
			break;
		}

		i = m;

		if (i == j) {
			break; // singleton interval, or `a` not found
		}

		m = right_child(m);
	} while (/*m != "bottom" && */ LCP[m] == l);

	// final sanity check
	if (i != ij.i ? FVC[i] == a : S[SA[i] + l] == a) {
		ij.i = i;
		ij.j = j;
		/* Also return the length of the LCP interval including `a` and
		 * possibly even more characters. Note: l + 1 <= LCP[m] */
		ij.l = LCP[m];
		ij.m = m;
	} else {
		ij.i = ij.j = -1;
	}

	return ij;
}

/** @brief Compute the longest match of a query with the subject.
 *
 * The *longest match* is the core concept of `andi`. Its simply defined as the
 * longest prefix of a query Q appearing anywhere in the subject S. Talking
 * about genetic sequences, a match is a homologous region, likely followed by a
 * SNP.
 *
 * This function returns the interval for where the longest match of the query
 * can be found in the ESA. Thereto it expects a starting interval for the
 * search.
 *
 * @param query - The query sequence.
 * @param qlen - The length of the query. Should correspond to `strlen(query)`.
 * @param k - The starting index into the query.
 * @param ij - The LCP interval for the string `query[0..k]`.
 * @returns The LCP interval for the longest prefix.
 */
lcp_interval esa::get_match_from(const char *query, size_t qlen, saidx_t k,
								 lcp_interval ij) const
{
	if (ij.i == -1 && ij.j == -1) {
		return ij;
	}

	// fail early on singleton intervals.
	if (ij.i == ij.j) {

		// try to extend the match. See line 513 below.
		saidx_t p = SA[ij.i];
		size_t k = ij.l;
		// const char *S = this->S;

		for (; k < qlen && S[p + k]; k++) {
			if (S[p + k] != query[k]) {
				ij.l = k;
				return ij;
			}
		}

		ij.l = k;
		return ij;
	}

	saidx_t l, i, j;

	lcp_interval res = ij;

	// Loop over the query until a mismatch is found
	do {
		// Get the subinterval for the next character.
		ij = get_interval(ij, query[k]);
		i = ij.i;
		j = ij.j;

		// If our match cannot be extended further, return.
		if (i == -1 && j == -1) {
			res.l = k;
			return res;
		}

		res.i = ij.i;
		res.j = ij.j;

		l = qlen;
		if (i < j && ij.l < l) {
			/* Instead of making another look up we can use the LCP interval
			 * calculated in get_interval */
			l = ij.l;
		}

		// By definition, the kth letter of the query was matched.
		k++;

		// Extend the match
		for (int p = SA[i]; k < l; k++) {
			if (S[p + k] != query[k]) {
				res.l = k;
				return res;
			}
		}
	} while (k < (ssize_t)qlen);

	res.l = qlen;
	return res;
}

/** @brief Get a match.
 *
 * Given an ESA and a string Q find the longest prefix of Q that matches
 * somewhere in C. This search is done entirely via jumping around in the ESA,
 * and thus is slow.
 *
 * @param query - The query string — duh.
 * @param qlen - The length of the query.
 * @returns the lcp interval of the match.
 */
lcp_interval esa::get_match(const char *query, size_t qlen) const
{
	saidx_t m = left_child(m_size);
	lcp_interval ij = {.l = LCP[m], .i = 0, .j = m_size - 1, .m = m};

	return get_match_from(query, qlen, 0, ij);
}

/** @brief Compute the LCP interval of a query. For a certain prefix length of
 * the query its LCP interval is retrieved from a cache. Hence this is faster
 * than the naive `get_match`. If the cache fails to provide a proper value, we
 * fall back to the standard search.
 *
 * @param query - The query sequence.
 * @param qlen - The length of the query. Should correspond to `strlen(query)`.
 * @returns The LCP interval for the longest prefix.
 */
lcp_interval esa::get_match_cached(const char *query, size_t qlen) const
{
	if (qlen <= CACHE_LENGTH) return get_match(query, qlen);

	ssize_t offset = 0;
	for (size_t i = 0; i < CACHE_LENGTH && offset >= 0; i++) {
		offset <<= 2;
		offset |= char2code(query[i]);
	}

	if (offset < 0) {
		return get_match(query, qlen);
	}

	lcp_interval ij = cache[offset];

	if (ij.i == -1 && ij.j == -1) {
		return get_match(query, qlen);
	}

	return get_match_from(query, qlen, ij.l, ij);
}
