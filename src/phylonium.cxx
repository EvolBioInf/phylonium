/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 - 2019 © Fabian Klötzl
 */
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
#include <cmath>
#include <err.h>
#include <errno.h>
#include <functional>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <limits>
#include <string.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "config.h"
#include "global.h"
#include "io.h"
#include "process.h"
#include "sequence.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern gsl_rng *RNG;
double ANCHOR_P_VALUE = 0.025;
int FLAGS = flags::none;
int THREADS = 1;
long unsigned int BOOTSTRAP = 0;
int RETURN_CODE = EXIT_SUCCESS;

// TODO: This is a hack and should be redone at some point.
size_t reference_index = 0;

void usage(int);
void version(void);

using mat_type = std::vector<evo_model>;

size_t pick_second_pass(std::vector<sequence> &sequences,
						const mat_type &matrix);
size_t pick_first_pass(std::vector<sequence> &sequences);
void cleanup_names(const std::string &reference_name,
				   std::vector<std::string> &file_names);

int main(int argc, char *argv[])
{
	RNG = gsl_rng_alloc(gsl_rng_default);
	if (!RNG) {
		err(1, "RNG allocation failed.");
	}

	// seed the random number generator with the current time
	// TODO: enable seeding for reproducibility
	gsl_rng_set(RNG, time(NULL));

	int version_flag = 0;
	bool two_pass = false;
	auto reference_name = std::string{};

	enum { P_AUTO, P_NEVER, P_ALWAYS } progress = P_AUTO;

	static struct option long_options[] = {
		{"2pass", no_argument, NULL, '2'},
		{"bootstrap", required_argument, NULL, 'b'},
		{"complete-deletion", no_argument, NULL, 0},
		{"help", no_argument, NULL, 'h'},
		{"progress", optional_argument, NULL, 0},
		{"threads", required_argument, NULL, 't'},
		{"verbose", no_argument, NULL, 'v'},
		{"version", no_argument, &version_flag, 1},
		{0, 0, 0, 0}};

#ifdef _OPENMP
	// Use all available processors by default.
	THREADS = omp_get_num_procs();
#endif

	// parse arguments
	while (1) {
		int option_index;
		int c =
			getopt_long(argc, argv, "2b:hr:t:v", long_options, &option_index);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: {
				auto arg_str = std::string(long_options[option_index].name);
				if (arg_str == "complete-deletion") {
					FLAGS |= flags::complete_deletion;
					break;
				}
				if (arg_str == "progress") {
					if (!optarg || strcasecmp(optarg, "always") == 0) {
						progress = P_ALWAYS;
					} else if (strcasecmp(optarg, "auto") == 0) {
						progress = P_AUTO;
					} else if (strcasecmp(optarg, "never") == 0) {
						progress = P_NEVER;
					} else {
						warnx("invalid argument to --progress '%s'. Expected "
							  "one of 'auto', 'always', or 'never'.",
							  optarg);
					}
					break;
				}
				break;
			}
			case '2': {
				two_pass = true;
				break;
			}
			case 'b': {
				errno = 0;
				char *end;
				long unsigned int bootstrap = strtoul(optarg, &end, 10);

				if (errno || end == optarg || *end != '\0' || bootstrap == 0) {
					soft_errx(
						"Expected a positive number for -b argument, but '%s' "
						"was given. Ignoring -b argument.",
						optarg);
					break;
				}

				BOOTSTRAP = bootstrap - 1;
				break;
			}
			case 'h': usage(EXIT_SUCCESS); break;
			case 'r': {
				reference_name = optarg;
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

	// determine whether to print a progress bar
	if (progress == P_AUTO) {
		progress = isatty(STDERR_FILENO) ? P_ALWAYS : P_NEVER;
	}
	if (progress == P_ALWAYS) {
		FLAGS |= flags::print_progress;
	}

	argc -= optind;
	argv += optind;

	auto file_names = std::vector<std::string>(argv, argv + argc);

	// add missing reference names to file names
	if (reference_name != "") {
		cleanup_names(reference_name, file_names);
	}

	if (file_names.size() < 2) {
		usage(EXIT_FAILURE);
	}

	// at max `file_names` many files have to be read.
	auto queries = std::vector<sequence>(file_names.size());

	// read all genomes
	#pragma omp parallel for num_threads(THREADS)
	for (size_t i = 0; i < file_names.size(); i++) {
		queries[i] = join(read_genome(file_names[i]));
	}

	// avoid copying sequences
	if (reference_name == "") {
		pick_first_pass(queries);
	} else {
		auto it =
			std::find(file_names.begin(), file_names.end(), reference_name);
		auto index = it - file_names.begin();
		reference_index = index;
	}

	auto matrix = process(queries[reference_index], queries);

	if (two_pass) {
		auto new_reference_index = pick_second_pass(queries, matrix);
		auto new_matrix = process(queries[new_reference_index], queries);

		print_matrix(queries, new_matrix);
	} else {
		print_matrix(queries, matrix);
	}

	return RETURN_CODE;
}

/** @brief Picks the references for the second pass according to some criterion.
 *
 * For the second pass we make an informed choice about the reference(s) using
 * the distances computed in the first pass. There are a number of possible
 * options:
 *
 *  - ingroup
 *  - outgroup
 *  - best coverage
 *
 * At the moment the most central sequence is chosen.
 *
 * @param sequences - The input sequences.
 * @param matrix - The distance matrix from the first pass.
 * @returns a new reference (index).
 */
size_t pick_second_pass(std::vector<sequence> &sequences,
						const std::vector<evo_model> &matrix)
{
	auto ret = std::vector<std::reference_wrapper<sequence>>();

	auto size = sequences.size();

	auto dist_matrix = std::vector<double>(size * size, NAN);
	std::transform(std::begin(matrix), std::end(matrix),
				   std::begin(dist_matrix),
				   [](const evo_model &em) { return em.estimate_JC(true); });

	auto central_value = std::numeric_limits<double>::max();
	auto central_index = (size_t)0;
	for (size_t i = 0; i < size; i++) {
		auto sum = std::accumulate(dist_matrix.begin() + i * size,
								   dist_matrix.begin() + i * size + size, 0.0);

		if (sum < central_value) {
			central_value = sum;
			central_index = i;
		}
	}
	ret.push_back(sequences[central_index]);
	reference_index = central_index;

	return central_index;
}

/** @brief Picks the references for the first pass according to some criterion.
 *
 * For the first pass we make an best-effort choice for a suitable reference if
 * none was supplied by the user. There are a number of possible options:
 *
 *  - size (smallest, medium, largest)
 *  - gc content (medium)
 *  - assembly quality
 *
 * At the moment a sequence of medium length is chosen.
 *
 * @param sequences - The input sequences.
 * @returns a new reference (index).
 */
size_t pick_first_pass(std::vector<sequence> &sequences)
{
	// pick a reference by some criterion.
	auto ret = std::vector<std::reference_wrapper<sequence>>(sequences.begin(),
															 sequences.end());

	std::nth_element(ret.begin(), ret.begin() + ret.size() / 2, ret.end(),
					 [](const sequence &a, const sequence &b) {
						 return a.size() < b.size();
					 });

	auto &reference = ret[ret.size() / 2].get();

	// recover which element the reference was in the original array
	auto it = std::find(sequences.begin(), sequences.end(), reference);
	reference_index = it - sequences.begin();

	if (FLAGS & flags::verbose) {
		std::cerr << "chosen reference: " << reference.get_name() << std::endl;
	}

	return reference_index;
}

void cleanup_names(const std::string &reference_name,
				   std::vector<std::string> &file_names)
{
	file_names.push_back(reference_name);
	std::sort(file_names.begin(), file_names.end());
	auto split = std::unique(file_names.begin(), file_names.end());
	file_names.erase(split, file_names.end());
}

/** @brief Prints the usage and then exits.
 * @param status - The return status.
 * @returns - It doesn't.
 */
void usage(int status)
{
	const char str[] = {
		"Usage: phylonium [OPTIONS] FILES...\n"
		"\tFILES... can be any sequence of FASTA files, each file representing "
		"one genome.\n\n"
		"Options:\n"
		"  -2, --2pass          Enable two-pass algorithm\n"
		"  -b, --bootstrap=N    Print additional bootstrap matrices\n"
		"  --complete-deletion  Delete the whole aligned column in case of "
		"gaps\n"
		"  -r FILE              Set the reference genome\n"
#ifdef _OPENMP
		"  -t, --threads=N      The number of threads to be used; by default, "
		"all available processors are used\n"
#endif
		"  -v, --verbose        Print additional information\n"
		"  -h, --help           Display this help and exit\n"
		"      --version        Output version information and "
		"acknowledgments\n"};
	fprintf(status == EXIT_SUCCESS ? stdout : stderr, "%s", str);
	exit(status);
}

/**
 * @brief This function just prints the version string and then aborts
 * the program.
 */
void version(void)
{
	const char str[] = {
		"phylonium " VERSION "\n"
		"Copyright (C) 2017 - 2019 Fabian Klötzl\n"
		"License GPLv3+: GNU GPL version 3 or later "
		"<http://gnu.org/licenses/gpl.html>\n"
		"This is free software: you are free to change and redistribute it.\n"
		"There is NO WARRANTY, to the extent permitted by law.\n\n"};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}
