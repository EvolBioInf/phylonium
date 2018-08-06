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

/** @brief A class to hold DNA sequences. */
class sequence
{
	/** An identifying name for the sequence. Commonly taken from FASTA header. */
	std::string name{};
	/** The DNA sequence. */
	std::string nucl{};
	/** Length of the DNA sequence. */
	size_t length{0};

  public:
	/** Reasonable default constructor. */
	sequence() = default;
	sequence(std::string, std::string) noexcept;

	std::string to_fasta() const;

	/** @brief Returns the length of the sequence.
	 * @returns the length of the sequence.
	 */
	auto size() const noexcept
	{
		return length;
	}

	/** @brief Return the name.
	 * @returns a copy of the name.
	 */
	auto get_name() const
	{
		return name;
	}

	/** @brief Return the nucleotide sequence.
	 * @returns a constant reference to the nucleotide sequence.
	 */
	const auto &get_nucl() const
	{
		return nucl;
	}

	/** @brief Return an iterator to the beginning of the nucleotide sequence.
	 * @returns an iterator to the beginning of the nucleotide sequence.
	 */
	auto begin() const
	{
		return nucl.begin();
	}

	/** @brief Return an iterator to the end of the nucleotide sequence.
	 * @returns an iterator to the end of the nucleotide sequence.
	 */
	auto end() const
	{
		return nucl.end();
	}

	/** @brief Give access to the raw nucleotide string.
	 * @returns a raw pointer to the nucleotides.
	 */
	auto c_str() const
	{
		return nucl.c_str();
	}
};

std::string reverse(const std::string &);
double gc_content(const std::string &) noexcept;
std::string filter_nucl(const std::string &);

/**
 * @brief A class to hold all the sequences from one genome/FASTA file.
 */
class genome
{
	/** ID for the genome (usually stripped file name). */
	std::string name{};
	/** Set of sequences. */
	std::vector<sequence> contigs{};
	/** Sum of all sequences plus border characters. */
	size_t joined_length{0};

  public:
	/** @brief Reasonable default constructor. */
	genome() = default;

	/** @brief Reasonable constructor.
	 * @param _name - New name of the genome.
	 * @param _contigs - New set of sequences.
	 */
	genome(std::string _name, std::vector<sequence> _contigs) noexcept
		: name{std::move(_name)}, contigs{std::move(_contigs)}
	{
		joined_length = std::accumulate(
			begin(contigs), end(contigs), contigs.size() - 1,
			[](auto sum, const auto &contig) { return sum + contig.size(); });
	}

	/** @brief Return the name.
	 * @returns reference to the name.
	 */
	std::string &get_name()
	{
		return name;
	}

	/** @brief Return the name.
	 * @returns a constant reference to the name.
	 */
	const std::string &get_name() const
	{
		return name;
	}

	/** @brief Return the set of sequences.
	 * @returns reference to the set of sequences.
	 */
	std::vector<sequence> &get_contigs()
	{
		return contigs;
	}

	/** @brief Return the set of sequences.
	 * @returns a constant reference to the set of sequences.
	 */
	const std::vector<sequence> &get_contigs() const
	{
		return contigs;
	}

	/** @brief Return the joined length.
	 * @returns the joined length.
	 */
	auto get_length() const noexcept
	{
		return joined_length;
	}
};

sequence join(const genome &gen);
