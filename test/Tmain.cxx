#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "global.h"
#include <random>

double ANCHOR_P_VALUE = 0.025;
int FLAGS = flags::none;
int THREADS = 1;
long unsigned int BOOTSTRAP = 0;
int RETURN_CODE = EXIT_SUCCESS;
size_t reference_index = 0;
std::mt19937 prng;
