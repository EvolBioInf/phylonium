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

	int file_descriptor = open(file_name, O_RDONLY);

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

	pfasta_free(&pf);
	close(file_descriptor);

	return sequences;
}

genome read_genome(std::string file_name)
{
	std::string species{extract_genome(file_name.c_str())};

	return genome{species, read_fasta(file_name, "")};
}

void just_print(const std::vector<sequence> &queries,
				const std::vector<std::string> &names,
				const std::vector<evo_model> &matrix, bool warnings = false)
{

	auto N = names.size();
	auto dist_matrix = std::vector<double>(N * N, NAN);
	std::transform(std::begin(matrix), std::end(matrix),
				   std::begin(dist_matrix),
				   [](const evo_model &em) { return em.estimate_JC(); });

	std::cout << N << std::endl;
	for (size_t i = 0; i < N; i++) {
		std::cout << names[i];
		for (size_t j = 0; j < N; j++) {
			auto index = i * N + j;
			auto dist = (i == j ? 0.0 : dist_matrix[index]);

			std::cout << "  " << std::setw(8) << std::setprecision(4) << dist;

			if (std::isnan(dist) && warnings) {
				const char str[] = {
					"For the two sequences '%s' and '%s' the distance "
					"computation failed and is reported as nan."};
				soft_errx(str, names[i].c_str(), names[j].c_str());
			}

			if (!std::isnan(dist) && i < j && warnings) {
				double coverage1 = matrix[index].coverage(queries[i].size());
				double coverage2 = matrix[index].coverage(queries[j].size());

				if (coverage1 < 0.2 || coverage2 < 0.2) {
					const char str[] = {
						"For the two sequences '%s' and '%s' less than 20%% "
						"homology were found (%f and %f, respectively)."};
					soft_errx(str, names[i].c_str(), names[j].c_str(),
							  coverage1, coverage2);
				}
			}
		}

		std::cout << std::endl;
	}
}

void print_matrix(const std::vector<sequence> &queries,
				  const std::vector<evo_model> &matrix)
{
	auto N = queries.size();
	auto names = std::vector<std::string>(N);
	std::transform(std::begin(queries), std::end(queries), std::begin(names),
				   [&](const sequence &seq) { return seq.get_name(); });

	just_print(queries, names, matrix, true);
	if (BOOTSTRAP) {
		auto neu = std::vector<evo_model>(N * N);
		for (auto k = 0; k < BOOTSTRAP; k++) {
			std::transform(std::begin(matrix), std::end(matrix),
						   std::begin(neu),
						   [](const evo_model &em) { return em.bootstrap(); });
			just_print(queries, names, neu);
		}
	}
}

void print_matrix(const sequence &subject, const std::vector<sequence> &queries,
				  const std::vector<evo_model> &matrix)
{
	auto N = queries.size();

	print_matrix(queries, matrix);

	if (FLAGS & flags::verbose) {
		auto M = [&matrix, N = N](size_t i, size_t j) -> const evo_model & {
			return matrix[i * N + j];
		};
		auto covf = std::ofstream(subject.get_name() + ".abscov");
		covf << "Absolute Coverages:\n";
		for (size_t i = 0; i < N; i++) {
			covf << queries[i].get_name();

			for (size_t j = 0; j < N; j++) {
				covf << "  " << std::setw(8) << std::setprecision(4)
					 << M(i, j).total();
			}
			covf << std::endl;
		}
	}
}
