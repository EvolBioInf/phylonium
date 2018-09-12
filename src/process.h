/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2018 © Fabian Klötzl
 */
#pragma once

#include <vector>
#include "sequence.h"

std::vector<evo_model> process(const sequence &, const std::vector<sequence> &);
