#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iterator>

extern gsl_rng *RNG;

class evo_model
{
  protected:
	size_t substitutions = 0;
	size_t homologs = 0;

  public:
	static int hash(char c) noexcept;

	evo_model() = default;

	void account(char a, char b) noexcept;
	void account(const char *stra, const char *strb, size_t length) noexcept;
	void account_rev(const char *stra, const char *strb, size_t b_offset,
					 size_t length) noexcept;
	evo_model &operator+=(const evo_model &other) noexcept;
	size_t total() const noexcept;
	double estimate_raw() const noexcept;
	double estimate_JC() const noexcept;

	static const evo_model &select_by_total(const evo_model &self,
											const evo_model &other)
	{
		return self.total() < other.total() ? other : self;
	}

	evo_model bootstrap() const;
};
