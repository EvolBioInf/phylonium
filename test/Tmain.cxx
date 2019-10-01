#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "global.h"

double ANCHOR_P_VALUE = 0.025;
int FLAGS = flags::none;
int THREADS = 1;
long unsigned int BOOTSTRAP = 0;
int RETURN_CODE = EXIT_SUCCESS;
size_t reference_index = 0;
