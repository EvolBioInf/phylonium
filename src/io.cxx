/**
 * @file
 * @brief This file contains the definitions for various io methods.
 */

#include "io.h"
#include <algorithm>
#include <cmath>
#include <err.h>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <pfasta.h>
#include <string>
#include <unistd.h>
#include <vector>
#include "global.h"
#include "sequence.h"

/** @brief extracts the genome name from a file path
 *
 * We try to be clever about the genome name. Given the file
 * path we extract just the file name. ie. path/file.ext -> file
 * This obviously fails on Windows.
 *
 * @param s_file_name - The file path
 * @returns just the genome name
 */
std::string extract_genome(const std::string &s_file_name)
{
	// find the last path separator
	auto left = s_file_name.rfind('/');
	left = (left == std::string::npos) ? 0 : left + 1;
	// left is the position one of to the right of the path separator

	// find the extension
	auto right = s_file_name.find('.', left);
	right = (right == std::string::npos) ? s_file_name.size() : right;

	// copy only the file name, not its path or extension
	return s_file_name.substr(left, right - left);
}

/**
 * @brief This function reads sequences from a file.
 * @param file_name - The file to read.
 * @param dsa - (output parameter) An array that holds found sequences.
 */
std::vector<sequence> read_fasta(std::string s_file_name, std::string prefix)
{
	std::vector<sequence> sequences{};
	const char *file_name = s_file_name.c_str();

	int file_descriptor =
		s_file_name != "-" ? open(file_name, O_RDONLY) : STDIN_FILENO;

	if (file_descriptor < 0) {
		err(errno, "%s", file_name);
	}

	int l;
	pfasta_file pf;

	if ((l = pfasta_parse(&pf, file_descriptor)) != 0) {
		errx(1, "%s: %s", file_name, pfasta_strerror(&pf));
	}

	pfasta_seq ps;
	while ((l = pfasta_read(&pf, &ps)) == 0) {
		sequences.emplace_back(prefix + ps.name, filter_nucl(ps.seq));
		pfasta_seq_free(&ps);
	}

	if (l < 0) {
		errx(1, "%s: %s", file_name, pfasta_strerror(&pf));
	}

fail:
	pfasta_free(&pf);
	close(file_descriptor);

	return sequences;
}

genome read_genome(std::string file_name)
{
	std::string species{extract_genome(file_name.c_str())};

	return genome{species, read_fasta(file_name, "")};
}

void print_matrix(const std::vector<sequence> &queries,
				  const std::vector<evo_model> &matrix)
{
	auto N = queries.size();

	auto dist_matrix = std::vector<double>(N * N, NAN);
	std::transform(std::begin(matrix), std::end(matrix),
				   std::begin(dist_matrix),
				   [](const evo_model &em) { return em.estimate_JC(); });

	std::cout << N << std::endl;
	for (size_t i = 0; i < N; i++) {
		std::cout << queries[i].get_name();
		for (size_t j = 0; j < N; j++) {
			std::cout << "  " << std::setw(8) << std::setprecision(4)
					  << (i == j ? 0.0 : dist_matrix[i * N + j]);
		}
		std::cout << std::endl;
	}
}

void print_matrix(const sequence &subject, const std::vector<sequence> &queries,
				  const std::vector<evo_model> &matrix)
{
	auto N = queries.size();

	auto M = [&matrix, N = N](size_t i, size_t j) -> const evo_model & {
		return matrix[i * N + j];
	};

	auto dist_matrix = std::vector<double>(N * N, NAN);
	std::transform(std::begin(matrix), std::end(matrix),
				   std::begin(dist_matrix),
				   [](const evo_model &em) { return em.estimate_JC(); });

	std::cout << N << std::endl;
	for (size_t i = 0; i < N; i++) {
		std::cout << queries[i].get_name();
		for (size_t j = 0; j < N; j++) {
			std::cout << "  " << std::setw(8) << std::setprecision(4)
					  << (i == j ? 0.0 : dist_matrix[i * N + j]);
		}
		std::cout << std::endl;
	}

	if (FLAGS & flags::verbose) {
		auto covf = std::ofstream(subject.get_name() + ".abscov");
		covf << "Absolute Coverages:\n";
		for (size_t i = 0; i < N; i++) {
			covf << queries[i].get_name();

			for (size_t j = 0; j < N; j++) {
				// TODO: queries[j].get_length() wrong length!
				covf << "  " << std::setw(8) << std::setprecision(4)
					 << M(i, j).total() /* / ((double)queries[j].size())*/
					;
			}
			covf << std::endl;
		}
	}
}
