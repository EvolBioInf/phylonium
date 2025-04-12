/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 - 2021 © Fabian Klötzl
 */
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>

/** @brief An "evolutionary model". Basically this class counts mutations. */
class evo_model
{
  protected:
	/** Count of substitutions among homologous nucleotides. */
	size_t substitutions = 0;
	/** Number of homologous nucleotides. */
	size_t homologs = 0;

  public:
	static int hash(char c) noexcept;

	/** @brief Reasonable default constructor. */
	evo_model() = default;

	void account(char a, char b) noexcept;
	void account(const char *stra, const char *strb, size_t length) noexcept;
	void account_rev(const char *stra, const char *strb, size_t b_offset,
					 size_t length) noexcept;
	evo_model &operator+=(const evo_model &other) noexcept;
	size_t total() const noexcept;
	double estimate_raw(bool zero_on_error = false) const noexcept;
	double estimate_JC(bool zero_on_error = false) const noexcept;
	double estimate_ani(bool zero_on_error = false) const noexcept;
	double coverage(size_t length) const noexcept;

	/** @brief Compare two counts of homologous nucleotides by length.
	 * @param self - This count.
	 * @param other - The count to compare to.
	 * @returns the model with the larger count of homologous nucleotides.
	 */
	static const evo_model &select_by_total(const evo_model &self,
											const evo_model &other)
	{
		return self.total() < other.total() ? other : self;
	}

	evo_model bootstrap() const;
};
