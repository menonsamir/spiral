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
#include <cstring>

#include "hexl/ntt/ntt.hpp"

#ifndef __EMSCRIPTEN__
    // #include <omp.h>
#endif

using namespace std;

#if defined(__AVX2__)
    #define USE_AVX2 1
    #include <immintrin.h>
#endif

#include <constants.h>
#include "values.h"
#include "poly.h"

#ifdef __EMSCRIPTEN__
    #include <emscripten.h>
#endif


extern intel::hexl::NTT *ntt_p;
extern intel::hexl::NTT *ntt_a;
extern intel::hexl::NTT *ntt_b;
extern intel::hexl::NTT *ntt_qprime;

uint64_t read_arbitrary_bits(const uint64_t *p, size_t bit_offs, size_t num_bits);
void write_arbitrary_bits(uint64_t *p, uint64_t val, size_t bit_offs, size_t num_bits);

void naive_product_crt_ntt(
    uint64_t *z, const uint64_t *x, const uint64_t *y,
    size_t rs, size_t ms, size_t cs
);

void build_table();
int64_t sample();

void ntt_forward(uint64_t *operand_overall);
void ntt_inverse(uint64_t *operand_overall);

typedef struct __attribute__((packed, aligned(16))) {
    uint64_t x, y;
} ulonglong2_h;

inline ulonglong2_h umul64wide(uint64_t a, uint64_t b) {
    ulonglong2_h z;
    __uint128_t val = ((__uint128_t)a) * ((__uint128_t)b);
    z.x = (uint64_t)val;
    z.y = (uint64_t)(val >> 64);
    return z;
}