/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 - 2020 © Fabian Klötzl
 */
#pragma once

#include <assert.h>
#include <vector>
#include "evo_model.h"
#include "sequence.h"

std::vector<evo_model> process(const sequence &, const std::vector<sequence> &);

class span
{
	size_t m_pos = 0, m_length = 0;

  public:
	span() = default;
	span(size_t p, size_t l) : m_pos{p}, m_length{l}
	{
	}

	auto start() const noexcept
	{
		return m_pos;
	}

	auto end() const noexcept
	{
		return m_pos + m_length;
	}

	auto length() const noexcept
	{
		return m_length;
	}

	bool starts_left_of(const span &other) const noexcept
	{
		return start() < other.start();
	}

	bool ends_left_of(const span &other) const noexcept
	{
		return end() <= other.start();
	}

	bool overlaps(span other) const noexcept
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

	span trim(size_t start, size_t end) const
	{
		if (end <= start) return *this; // invalid input range, ignore.
		if (!overlaps(span(start, end - start))) return span{};

		auto offset = start > this->start() ? start - this->start() : 0;
		auto drift = this->end() > end ? this->end() - end : 0;
		// Carefully handle cases where the given range is bigger than *this.
		auto that = *this;
		// shorten left by offset, right by drift.

		assert(this->m_length > offset + drift);
		that.m_pos += offset;
		that.m_length = this->m_length - offset - drift;

		return that;
	}

	span trim(span other) const
	{
		return trim(other.start(), other.end());
	}

	span trim_unsafe(size_t start, size_t end) const
	{
		// Carefully handle cases where the given range is bigger than *this.
		auto offset = start > this->start() ? start - this->start() : 0;
		auto drift = this->end() > end ? this->end() - end : 0;

		auto that = *this;

		// shorten left by offset, right by drift.
		that.m_pos += offset;
		that.m_length = this->m_length - offset - drift;

		return that;
	}

	span trim_unsafe(span other) const
	{
		return trim_unsafe(other.start(), other.end());
	}

	span shift_left(size_t offset) const
	{
		auto pos = start() >= offset ? start() - offset : 0;
		return {pos, length()};
	}

	static span from_start_end(size_t s, size_t e)
	{
		return span(s, e - s);
	}
};

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

	std::vector<span> anchors = {};

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

	homology copy_without_anchors() const
	{
		auto that = homology{};
		that.direction = direction;
		that.index_reference = index_reference;
		that.index_reference_projected = index_reference_projected;
		that.index_query = index_query;
		that.length = length;
		return that;
	}

	homology trim(size_t start, size_t end) const
	{
		if (end <= start) return *this; // invalid input range, ignore.

		// Carefully handle cases where the given range is bigger than *this.
		auto that = copy_without_anchors();
		auto offset = start > this->start() && start < this->end()
						  ? start - this->start()
						  : 0;
		auto drift =
			this->end() > end && end > this->start() ? this->end() - end : 0;
		that.index_reference_projected += offset;

		// shorten left by offset, right by drift.

		if (direction == dir::forward) {
			that.index_reference += offset;
			that.index_query += offset;
		} else {
			that.index_reference += drift;
			that.index_query += drift;
		}

		assert(this->length > offset + drift);
		that.length = this->length - offset - drift;

		auto interval = span(offset, that.length);

		that.anchors.clear();
		for (auto anchor : this->anchors) {
			if (anchor.overlaps(interval)) {
				that.anchors.push_back(
					anchor.trim_unsafe(interval).shift_left(offset));
			}
		}

		// auto it = anchors.begin();
		// auto out = std::back_inserter(that.anchors);
		// while (it < anchors.end() && it->ends_left_of(interval)) {
		// 	it++;
		// }
		// if (it < anchors.end() && it->overlaps(interval)) {
		// 	*out = it->trim(interval).shift_left(offset);
		// 	it++, out++;
		// }
		// while(it < anchors.end() && !interval.ends_left_of(*it)){
		// 	*out = it->trim(interval).shift_left(offset);
		// 	it++, out++;
		// }
		// that.anchors.erase(out, )

		return that;
	}
};

std::vector<std::vector<homology>>
complete_delete(const std::vector<std::vector<homology>> &homologies);
