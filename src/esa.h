/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 - 2019 © Fabian Klötzl
 */
/**
 * @file
 * @brief This header contains the declarations for functions in esa.cxx.
 *
 */
#pragma once

#include <memory>

#include <sys/types.h>

#include <divsufsort64.h>
using saidx_t = saidx64_t;

#include "sequence.h"

/**
 * @brief Represents LCP-Intervals.
 *
 * This struct is used to represent LCP-intervals. The member `i` should
 * coincide with the lower bound whereas `j` is the upper bound. Both bounds
 * are inclusive. So if `i == j` the interval contains exactly one element,
 * namely `i`. To represent an empty interval please use `i == j == -1`.
 * Other variants, such as `i == j == -2` can be used to indicate an error.
 * The common prefix length is denoted by l and should always be non-negative.
 * Variables of this type are often called `ij`.
 */
typedef struct {
	/** @brief The common prefix length */
	saidx_t l;
	/** @brief lower bound */
	saidx_t i;
	/** @brief upper bound */
	saidx_t j;
	/** The new middle. */
	saidx_t m;
} lcp_interval;

/** @brief A full text index.
 * Basically, this is a bunch of arrays working together for fast lookups.
 */
class esa
{
	/** Length of the arrays. */
	saidx_t m_size = 0;
	/** The sequence at the basis of the ESA. */
	sequence m_master{};

	/** The Longest Common Prefix array. */
	std::unique_ptr<saidx_t[]> LCP;
	/** The child array. */
	std::unique_ptr<saidx_t[]> CLD;
	/** The First Variant Character array. See Klötzl (2015) for an explanation.
	 */
	std::unique_ptr<char[]> FVC;
	/** A cache to speed up the look-up. */
	std::unique_ptr<lcp_interval[]> cache;

  public:
	/** The Suffix Array */
	std::unique_ptr<saidx_t[]> SA;
	/** The base string */
	std::string S{""};

	/** @brief a reasonable default constructor */
	esa() = default;
	esa(const sequence &);

	lcp_interval get_match(const char *, size_t) const;
	lcp_interval get_match_cached(const char *, size_t) const;

	/** @brief Get length of the arrays.
	 * @returns the length of the arrays.
	 */
	auto size() const noexcept
	{
		return m_size;
	}

  private:
	void init_LCP();
	void init_CLD();
	void init_FVC();
	void init_cache();
	void init_cache_dfs(char *, size_t, lcp_interval);
	void init_cache_fill(char *, size_t, lcp_interval);
	lcp_interval get_match_from(const char *, size_t, saidx_t,
								lcp_interval) const;
	lcp_interval get_interval(lcp_interval ij, char a) const;

	/** @brief Get the right child of a node.
	 * @param idx - Index of the node.
	 * @returns the index of the right child node.
	 */
	auto &right_child(ssize_t idx) const noexcept
	{
		return CLD[idx];
	}

	/** @brief Get the left child of a node.
	 * @param idx - Index of the node.
	 * @returns the index of the left child node.
	 */
	auto &left_child(ssize_t idx) const noexcept
	{
		return CLD[idx - 1];
	}
};

#ifdef DEBUG

char code2char(ssize_t code);

#endif // DEBUG
