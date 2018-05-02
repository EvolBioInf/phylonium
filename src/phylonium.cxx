/**
 * @file
 *
 * This is the main file. It contains functions to parse the commandline
 * arguments, read files etc.
 *
 * @brief The main file
 * @author Fabian Klötzl
 *
 * @section License
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <err.h>
#include <errno.h>
#include <getopt.h>
#include <string.h>
#include <unistd.h>

#include "config.h"
#include "global.h"
#include "io.h"
#include "process.h"
#include "sequence.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int FLAGS = flags::none;
int THREADS = 1;
double RANDOM_ANCHOR_PROP = 0.025;

void usage(int);
void version(void);

int main(int argc, char *argv[])
{
	int version_flag = 0;
	auto reference_names = std::vector<std::string>();

	static struct option long_options[] = {
		{"version", no_argument, &version_flag, 1},
		{"help", no_argument, NULL, 'h'},
		{"verbose", no_argument, NULL, 'v'},
		{"threads", required_argument, NULL, 't'},
		{0, 0, 0, 0}};

#ifdef _OPENMP
	// Use all available processors by default.
	THREADS = omp_get_num_procs();
#endif

	// parse arguments
	while (1) {
		int c = getopt_long(argc, argv, "hr:t:v", long_options, NULL);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: break;
			case 'h': usage(EXIT_SUCCESS); break;
			case 'r': {
				reference_names.push_back(optarg);
				break;
			}
			case 't': {
#ifdef _OPENMP
				errno = 0;
				char *end;
				long unsigned int threads = strtoul(optarg, &end, 10);

				if (errno || end == optarg || *end != '\0') {
					warnx("Expected a number for -t argument, but '%s' was "
						  "given. Ignoring -t argument.",
						  optarg);
					break;
				}

				if (threads > (long unsigned int)omp_get_num_procs()) {
					warnx(
						"The number of threads to be used, is greater then the "
						"number of available processors; Ignoring -t %lu "
						"argument.",
						threads);
					break;
				}

				THREADS = threads;
#else
				warnx(
					"This version of phylonium was built without OpenMP "
					"and thus "
					"does not support multi threading. Ignoring -t argument.");
#endif
				break;
			}
			case 'v':
				FLAGS |= FLAGS & flags::verbose ? flags::extra_verbose
												: flags::verbose;
				break;
			case '?': /* intentional fall-through */
			default: usage(EXIT_FAILURE); break;
		}
	}

	if (version_flag) {
		version();
	}

	argc -= optind;
	argv += optind;

	auto extra_file_names = std::vector<std::string>(argv, argv + argc);

	// add missing reference names to file names
	std::sort(reference_names.begin(), reference_names.end());
	std::sort(extra_file_names.begin(), extra_file_names.end());
	auto split = std::unique(reference_names.begin(), reference_names.end());
	reference_names.erase(split, reference_names.end());
	split = std::unique(extra_file_names.begin(), extra_file_names.end());
	extra_file_names.erase(split, extra_file_names.end());

	auto file_names = std::vector<std::string>();
	file_names.reserve(
		std::max(reference_names.size(), extra_file_names.size()));

	std::set_union(reference_names.begin(), reference_names.end(),
				   extra_file_names.begin(), extra_file_names.end(),
				   std::back_inserter(file_names));
	split = std::unique(file_names.begin(), file_names.end());
	file_names.erase(split, file_names.end());

	if (file_names.size() < 2) {
		if (!isatty(STDIN_FILENO)) {
			// if no files are supplied, read from stdin
			file_names.push_back("-");
		} else {
			// print a helpful message on './phylonium' without args
			usage(EXIT_FAILURE);
		}
	}

	if (reference_names.empty()) {
		errx(1, "no reference given.");
	}

	/*
	if (reference_names.empty()) {
		// pick longest genome
		auto longest_genome_it = max_element(
			begin(genomes), end(genomes), [](const genome &a, const genome &b) {
				return a.get_length() < b.get_length();
			});


		auto ref_file_name_it =
			std::find(file_names.begin(), file_names.end(), ref_name);
	}*/

	// at max `file_names` many files have to be read.
	auto genomes = std::vector<genome>();
	auto queries = std::vector<sequence>();
	genomes.reserve(file_names.size());
	queries.reserve(genomes.size());

	// read all genomes
	std::transform(file_names.begin(), file_names.end(),
				   std::back_inserter(genomes), read_genome);
	std::transform(genomes.begin(), genomes.end(), std::back_inserter(queries),
				   join);

	if (genomes.size() < 2) {
		errx(1, "Less than two genomes given, nothing to compare.");
	}

	auto references = std::vector<std::reference_wrapper<sequence>>();
	for (auto ref_name : reference_names) {
		auto it = std::find(file_names.begin(), file_names.end(), ref_name);
		auto index = it - file_names.begin();
		references.push_back(queries[index]);
	}

	using mat_type = std::vector<evo_model>;
	auto matrices = std::vector<mat_type>();

	for (auto subject : references) {
		auto matrix = process(subject, queries);
		matrices.push_back(matrix);
		// print_matrix(subject, queries, matrix);
	}

	auto super_matrix = mat_type(matrices[0].size());

	for (const auto &matrix : matrices) {
		// for (ssize_t i = 0; i < matrix.size(); i++) {
		// 	super_matrix[i] =
		// 		evo_model::select_by_total(super_matrix[i], matrix[i]);
		// }

		for (ssize_t i = 0; i < matrix.size(); i++) {
			super_matrix[i] =
				evo_model::select_by_total(super_matrix[i], matrix[i]);
		}
	}

	print_matrix(queries, super_matrix);

	return 0;
}

/** @brief
 * Prints the usage to stdout and then exits successfully.
 */
void usage(int status)
{
	const char str[] = {
		"Usage: phylonium [OPTIONS] FILES...\n"
		"\tFILES... can be any sequence of FASTA files.\n"
		"\tUse '-' to force reading from standard input.\n\n"
		"Options:\n"
		"  -r FILE           Add FILE to the list of references\n"
		"  -v, --verbose     Prints additional information\n"
#ifdef _OPENMP
		"  -t, --threads <INT> \n"
		"                    The number of threads to be used; by default, all "
		"available processors are used\n"
#endif
		"  -h, --help        Display this help and exit\n"
		"      --version     Output version information and acknowledgments\n"};

	fprintf(status == EXIT_SUCCESS ? stdout : stderr, "%s", str);
	exit(status);
}

/**
 * This function just prints the version string and then aborts
 * the program.
 */
void version(void)
{
	const char str[] = {
		"phylonium " VERSION "\n"
		"Copyright (C) 2017 - 2018 Fabian Klötzl\n"
		"License GPLv3+: GNU GPL version 3 or later "
		"<http://gnu.org/licenses/gpl.html>\n"
		"This is free software: you are free to change and redistribute it.\n"
		"There is NO WARRANTY, to the extent permitted by law.\n\n"};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}
