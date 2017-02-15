/**
 * @file
 * @brief This header contains the declarations for functions in esa.cxx.
 *
 */
#pragma once

#include <memory>

#include <sys/types.h>

#include <divsufsort.h>

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

class esa
{
	int m_size = 0;
	sequence m_master{};

	std::unique_ptr<saidx_t[]> LCP;
	std::unique_ptr<saidx_t[]> CLD;
	std::unique_ptr<char[]> FVC;
	std::unique_ptr<lcp_interval[]> cache;

  public:
	std::unique_ptr<saidx_t[]> SA;
	std::string S{""};

	esa() = default;
	esa(const sequence &);

	lcp_interval get_match(const char *, size_t) const;
	lcp_interval get_match_cached(const char *, size_t) const;

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

	auto &right_child(ssize_t idx) const noexcept
	{
		return CLD[idx];
	}

	auto &left_child(ssize_t idx) const noexcept
	{
		return CLD[idx - 1];
	}
};

#ifdef DEBUG

char code2char(ssize_t code);

#endif // DEBUG
