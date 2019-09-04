/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 - 2019 © Fabian Klötzl
 */
/**
 * @file
 * @brief This header contains function declarations for io procedures.
 */
#ifndef _IO_H_
#define _IO_H_

#include <evo_model.h>
#include <string>
#include <vector>
#include "sequence.h"

genome read_genome(std::string);
void print_matrix(const sequence &subject, const std::vector<sequence> &queries,
				  const std::vector<evo_model> &matrix);

void print_matrix(const std::vector<sequence> &queries,
				  const std::vector<evo_model> &matrix);

#endif // _IO_H_
