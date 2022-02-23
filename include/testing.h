#pragma once

#include <chrono>
#include <stack>
#include <fstream>
#include <inttypes.h>
#include <iostream>
#include <random>
#include <sstream>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <queue>
#include <mutex>
#include <assert.h>

using namespace std;

#include "hexl/ntt/ntt.hpp"

#include <bitset>
#include <map>
#include <math.h>

#include "core.h"
#include "poly.h"
#include "values.h"
#include "client.h"

#include "common.h"
#include <constants.h>

#include <functional>
#include <iomanip>

void pack(
    MatPoly& result,
    size_t out_n,
    size_t m_conv,
    const vector<MatPoly>& v_ct,
    const vector<MatPoly>& v_W
);
void testPacking(size_t out_n);
void testHighRate(size_t num_expansions, size_t further_dims, size_t idx_target);