/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 - 2019 © Fabian Klötzl
 */
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

extern size_t reference_index;

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
	auto right = s_file_name.rfind('.');
	if (right != std::string::npos) {
		auto ext = s_file_name.substr(right);
		if (ext == ".fa" || ext == ".fas" || ext == ".fasta") {
			// pass and strip
		} else {
			// don't strip unkown extention
			right = s_file_name.size();
		}
	} else {
		right = s_file_name.size();
	}

	// copy only the file name, not its path or extension
	return s_file_name.substr(left, right - left);
}

/**
 * @brief This function reads sequences from a file.
 * @param file_name - The file to read.
 * @param prefix - A prefix for the name.
 */
std::vector<sequence> read_fasta(std::string s_file_name, std::string prefix)
{
	std::vector<sequence> sequences{};
	const char *file_name = s_file_name.c_str();

	int file_descriptor = open(file_name, O_RDONLY);

	if (file_descriptor < 0) {
		err(errno, "%s", file_name);
	}

	auto parser = pfasta_init(file_descriptor);
	if (parser.errstr) {
		errx(1, "%s: %s", file_name, parser.errstr);
	}

	while (!parser.done) {
		auto record = pfasta_read(&parser);
		if (parser.errstr) {
			errx(1, "%s: %s", file_name, parser.errstr);
		}

		sequences.emplace_back(prefix + record.name,
							   filter_nucl(record.sequence));
		pfasta_record_free(&record);
	}

	pfasta_free(&parser);
	close(file_descriptor);

	return sequences;
}

genome read_genome(std::string file_name)
{
	std::string species{extract_genome(file_name.c_str())};

	return genome{species, read_fasta(file_name, "")};
}

void print_warnings(const std::vector<sequence> &queries,
					const std::vector<std::string> &names,
					const std::vector<double> &dist_matrix,
					const std::vector<evo_model> &matrix)
{
	auto N = names.size();

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < i; j++) {
			auto index = i * N + j;
			auto dist = dist_matrix[index];

			if (std::isnan(dist)) {
				const char str[] = {
					"For the two sequences '%s' and '%s' the distance "
					"computation failed and is reported as nan."};
				soft_errx(str, names[i].c_str(), names[j].c_str());
			}

			if (!std::isnan(dist)) {
				double cov1 = matrix[index].coverage(queries[i].size());
				double cov2 = matrix[index].coverage(queries[j].size());
				const char str[] = {
					"For the two sequences '%s' and '%s' less than 20%% "
					"homology were found (%f and %f, respectively)."};

				if (cov1 < 0.2 || cov2 < 0.2) {
					soft_errx(str, names[i].c_str(), names[j].c_str(), cov1,
							  cov2);
				}
			}
		}
	}
}

void just_print(const std::vector<std::string> &names,
				const std::vector<double> &dist_matrix)
{
	auto N = names.size();

	// Produce output in PHYLIP distance matrix format
	std::cout << N << std::endl;
	std::cout.precision(4);
	std::cout << std::scientific;

	for (size_t i = 0; i < N; i++) {
		std::cout << names[i];

		for (size_t j = 0; j < N; j++) {
			auto index = i * N + j;
			auto dist = (i == j ? 0.0 : dist_matrix[index]);

			std::cout << "  " << dist;
		}

		std::cout << std::endl;
	}
}

void print_matrix(const std::vector<sequence> &queries,
				  const std::vector<evo_model> &matrix)
{
	auto dist_getter =
		FLAGS & flags::dist_raw
			? [](const evo_model &em) { return em.estimate_raw(); }
		: FLAGS & flags::dist_ani
			? [](const evo_model &em) { return em.estimate_ani(); }
			: [](const evo_model &em) { return em.estimate_JC(); };

	auto N = queries.size();
	auto names = std::vector<std::string>(N);
	std::transform(std::begin(queries), std::end(queries), std::begin(names),
				   [&](const sequence &seq) { return seq.get_name(); });

	auto dist_matrix = std::vector<double>(N * N, NAN);
	std::transform(std::begin(matrix), std::end(matrix),
				   std::begin(dist_matrix), dist_getter);

	// print warnings before distance matrix
	print_warnings(queries, names, dist_matrix, matrix);

	just_print(names, dist_matrix);
	if (BOOTSTRAP) {
		auto neu = std::vector<evo_model>(N * N);
		for (auto k = (size_t)0; k < BOOTSTRAP; k++) {
			std::transform(std::begin(matrix), std::end(matrix),
						   std::begin(neu),
						   [](const evo_model &em) { return em.bootstrap(); });

			std::transform(std::begin(neu), std::end(neu),
						   std::begin(dist_matrix), dist_getter);

			just_print(names, dist_matrix);
		}
	}

	double sum = 0.0;
	size_t counter = 0;
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < i; j++) {
			auto index = i * N + j;
			auto dist = dist_matrix[index];

			if (std::isnan(dist)) continue;
			double cov1 = matrix[index].coverage(queries[i].size());
			double cov2 = matrix[index].coverage(queries[j].size());

			sum += cov1 + cov2;
			counter += 2;
		}
	}

	size_t aln_aligned = 0;
	size_t aln_total = 0;
	for (size_t i = 0; i < N; i++) {
		if (i == reference_index) continue;

		auto index = reference_index * N + i;
		aln_aligned += matrix[index].total();
		aln_total += queries[i].size();
	}

	if (FLAGS & flags::verbose) {
		std::cerr << "avg coverage:\t" << sum / counter << std::endl;
		std::cerr << "alignment:\t" << aln_aligned << "\t" << aln_total << "\t"
				  << aln_aligned / (double)aln_total << std::endl;
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
