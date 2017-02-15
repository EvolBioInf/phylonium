/*
 * Copyright (C) 2017  Fabian Klötzl
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <cassert>
#include <err.h>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <unistd.h>
#include <vector>

class matrix
{
  public:
	using size_type = size_t;

  protected:
	/// The size of the matrix
	size_type size = 0;
	/// A list of names (of sequences)
	std::vector<std::string> names = {};
	/// The matrix itself
	std::vector<double> values = {};

  public:
	matrix() = default;
	/** @brief Create a new matrix from a set of names and values.
	 *
	 * Take arguments by value, because compilers can optimize this anyway.
	 *
	 * @param _names - The new set of names
	 * @param _values - The new values
	 * @returns the new matrix
	 */
	matrix(std::vector<std::string> _names, std::vector<double> _values)
		: size{_names.size()}, names{std::move(_names)},
		  values{std::move(_values)}
	{
		assert(size == names.size());
		assert(size * size == values.size());
	}

	/** @brief Access an entry by coordinates.
	 *
	 * @param i - the row index
	 * @param j - the column index
	 * @returns a mutable reference to the entry
	 */
	double &entry(size_type i, size_type j)
	{
		return values[i * size + j];
	}

	/** @brief Access an entry by coordinates.
	 *
	 * @param i - the row index
	 * @param j - the column index
	 * @returns an immutable reference to the entry
	 */
	const double &entry(size_type i, size_type j) const
	{
		return values[i * size + j];
	}

	/** @brief Get an iterator to the beginning of a row.
	 *
	 * @param i - the row index
	 * @returns An iterator.
	 */
	auto row(size_type i) noexcept
	{
		return values.begin() + i * size;
	}

	/** @brief Get an iterator to the beginning of a row.
	 *
	 * @param i - the row index
	 * @returns A read-only iterator.
	 */
	auto row(size_type i) const noexcept
	{
		return values.cbegin() + i * size;
	}

	/** @brief Get an iterator to the end of a row.
	 *
	 * @param i - the row index
	 * @returns An iterator.
	 */
	auto row_end(size_type i) noexcept
	{
		return values.begin() + (i + 1) * size;
	}

	/** @brief Get an iterator to the beginning of a row.
	 *
	 * @param i - the row index
	 * @returns A read-only iterator.
	 */
	auto row_end(size_type i) const noexcept
	{
		return values.cbegin() + (i + 1) * size;
	}

	/** @brief Get the size of the matrix
	 *
	 * @returns The size.
	 */
	auto get_size() const noexcept
	{
		return size;
	}

	/** @brief Get the list of names
	 *
	 * @returns Returns a read-only reference of the names.
	 */
	auto get_names() const noexcept -> const std::vector<std::string> &
	{
		return names;
	}

	// defined in matrix.cxx
	std::string to_string() const;
};

/** @brief Sample a distance matrix. Only the names referenced by index via the
 * given range are included in the new submatrix. The implied order of names
 * from the index list is preserved.
 *
 * @param self - The matrix to sample.
 * @param first - An iterator to a list of indices.
 * @param last - An iterator past the list of indices.
 * @returns a new matrix with only the names specified in the index list.
 */
template <typename InputIt>
matrix sample(const matrix &self, InputIt first, InputIt last)
{
	auto new_size = distance(first, last);

	const auto &self_names = self.get_names();
	auto new_names = std::vector<std::string>();
	auto inserter = std::back_inserter(new_names);
	new_names.reserve(new_size);

	// The new names consists of only those given in the index list
	std::transform(first, last, inserter,
				   [&self_names](matrix::size_type index) {
					   return self_names.at(index);
				   });

	auto new_values = std::vector<double>(new_size * new_size);
	auto ret = matrix{new_names, new_values};

	// Copy and rearrange the values
	for (auto it = first; it != last; it++) {
		auto old_it = *it;
		auto new_it = distance(first, it);
		for (auto jk = first; jk != last; jk++) {
			auto old_jk = *jk;
			auto new_jk = distance(first, jk);
			ret.entry(new_it, new_jk) = self.entry(old_it, old_jk);
		}
	}

	return ret;
}

// defined in matrix.cxx
std::vector<matrix> parse_all(const char *const *);
