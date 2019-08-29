/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 - 2019 © Fabian Klötzl
 */
#include <cstdlib>

enum flags {
	none,
	verbose,
	extra_verbose,
	complete_deletion = 4,
	print_progress = 8,
	print_positions = 16
};
extern int FLAGS;
extern int THREADS;
extern int RETURN_CODE;

extern double ANCHOR_P_VALUE;
extern long unsigned int BOOTSTRAP;

/**
 * @brief This macro is used to print a warning and make the program exit with
 * an failure exit code, later.
 */
#define soft_err(...)                                                          \
	do {                                                                       \
		RETURN_CODE |= EXIT_FAILURE;                                           \
		warn(__VA_ARGS__);                                                     \
	} while (0)

/**
 * @brief This macro is used to print a warning and make the program exit with
 * an failure exit code, later.
 */
#define soft_errx(...)                                                         \
	do {                                                                       \
		RETURN_CODE |= EXIT_FAILURE;                                           \
		warnx(__VA_ARGS__);                                                    \
	} while (0)
