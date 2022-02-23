#pragma once

#include <chrono>
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

#include "common.h"
#include <constants.h>

#include <functional>
#include <iomanip>

inline size_t get_bits_per(size_t dim) {
    if (dim == logQ) return 1;
    size_t bits_per = (size_t) floor(logQ / (double)dim) + 1;
    return bits_per;
}

inline size_t min(size_t a, size_t b) {
    return a > b ? b : a;
}

int64_t variance_raw(const vector<int64_t> &x);

double get_log_var(const MatPoly &A, const MatPoly &B, bool print = false, size_t mod = p_i);
vector<int64_t> get_diffs(const MatPoly &A, const MatPoly &B, uint64_t modulus);
void output_diffs(ostream &os, vector<int64_t> diffs);
void uniform_matrix(MatPoly &A);
void buildGadget(MatPoly &G);
MatPoly buildGadget(size_t rows, size_t cols);
void gadget_invert(size_t mx, MatPoly &out, const MatPoly &in, size_t rdim = n1);
MatPoly gadget_invert(size_t mx, size_t m, const MatPoly &in, size_t rdim = n1);
MatPoly mul_over_integers(const MatPoly &a, const MatPoly &b, uint64_t modulus, bool b_is_single_poly = false);
void to_ntt_qprime(MatPoly &a);
void from_ntt_qprime(MatPoly &a);
MatPoly mul_over_qprime(const MatPoly &a, const MatPoly &b);
int64_t inv_mod(int64_t a, int64_t b);