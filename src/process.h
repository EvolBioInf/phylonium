/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 - 2019 © Fabian Klötzl
 */
#pragma once

#include <assert.h>
#include <vector>
#include "evo_model.h"
#include "sequence.h"

std::vector<evo_model> process(const sequence &, const std::vector<sequence> &);


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

	auto start() const noexcept
	{
		return index_reference_projected;
	}

	auto end() const noexcept
	{
		return index_reference_projected + length;
	}

	auto start_query() const noexcept
	{
		return index_query;
	}

	auto end_query() const noexcept
	{
		return index_query + length;
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
		if (start() == other.start()) {
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
		return start() < other.start();
	}

	/** @brief Determine the order of two homologies.
	 * @param other - a homology to compare with.
	 * @returns true if this homology ends left of the argument wrt. the
	 * reference genome.
	 */
	bool ends_left_of(const homology &other) const noexcept
	{
		return end() <= other.start();
	}

	homology trim(size_t start, size_t end) const
	{
		if (end <= start) return *this; // invalid input range, ignore.

		// Carefully handle cases where the given range is bigger than *this.
		auto that = *this;
		auto offset = start > this->start() && start < this->end()
						  ? start - this->start()
						  : 0;
		auto drift =
			this->end() > end && end > this->start() ? this->end() - end : 0;
		that.index_reference_projected += offset;

		if (direction == dir::forward) {
			that.index_reference += offset;
			that.index_query += offset;
		} else {
			that.index_reference += drift;
			that.index_query += drift;
		}

		assert(this->length > offset + drift);
		that.length = this->length - offset - drift;
		return that;
	}
};

std::vector<std::vector<homology>> complete_delete(
	const std::vector<std::vector<homology>> &homologies);
