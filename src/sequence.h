/**
 * @file
 * @brief Functions and structures for DNA sequences
 *
 */
#pragma once

#include <memory>
#include <numeric>
#include <string>
#include <vector>

class sequence
{
	std::string name{};
	std::string nucl{};
	size_t length{0};

  public:
	sequence() = default;
	sequence(std::string, std::string) noexcept;

	std::string to_fasta() const;

	auto size() const noexcept
	{
		return length;
	}

	auto get_name() const
	{
		return name;
	}

	const auto &get_nucl() const
	{
		return nucl;
	}

	auto begin() const
	{
		return nucl.begin();
	}

	auto end() const
	{
		return nucl.end();
	}

	auto c_str() const
	{
		return nucl.c_str();
	}
};

std::string reverse(const std::string &);
double gc_content(const std::string &) noexcept;
std::string filter_nucl(const std::string &);

class genome
{
	std::string name{};
	std::vector<sequence> contigs{};
	size_t joined_length{0};

  public:
	genome() = default;

	genome(std::string _name, std::vector<sequence> _contigs) noexcept
		: name{std::move(_name)},
		  contigs{std::move(_contigs)}
	{
		joined_length = std::accumulate(
			begin(contigs), end(contigs), contigs.size() - 1,
			[](auto sum, const auto &contig) { return sum + contig.size(); });
	}

	std::string &get_name()
	{
		return name;
	}

	const std::string &get_name() const
	{
		return name;
	}

	std::vector<sequence> &get_contigs()
	{
		return contigs;
	}

	const std::vector<sequence> &get_contigs() const
	{
		return contigs;
	}

	auto get_length() const noexcept
	{
		return joined_length;
	}
};

sequence join(const genome &gen);
