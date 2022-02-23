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
#include "util.h"

#include "common.h"
#include <constants.h>

#include <functional>
#include <iomanip>

uint64_t sample_u64();
void noise(MatPoly &E);
void keygen(MatPoly &S, MatPoly &Sp, MatPoly &sr, size_t n_val = n0);
MatPoly get_fresh_public_key_raw(MatPoly &Sp, long m);
MatPoly enc_scalar(const MatPoly &P, const MatPoly &G, const MatPoly &sigma, size_t mm, size_t num_expansions);
void load_modswitched_into_ct(MatPoly &ct, const uint64_t *modswitched);
void dec_compressed(MatPoly &O_out, const MatPoly &C, const MatPoly &S, uint64_t divisor = 1);
void dec(MatPoly &Z, const MatPoly &S, const MatPoly &C);
MatPoly encryptSimpleRegev(MatPoly sigma, size_t noise_factor = 1);
MatPoly encryptSimpleRegevMatrix(const MatPoly &s, const MatPoly &G_conv, MatPoly sigma_m, size_t noise_factor = 1);
MatPoly encryptSimpleRegevMatrix(const MatPoly &s, const MatPoly &mat_to_enc, size_t noise_factor = 1);
void generate_X_improved(size_t mx, MatPoly &X, MatPoly A, MatPoly G);
MatPoly encryptSimpleRegevMatrix(const MatPoly &s, const MatPoly &mat_to_enc, size_t noise_factor);
void generate_X_improved(size_t mx, MatPoly &X, MatPoly A, MatPoly G);
MatPoly getSkVec(size_t idx);
MatPoly getRegevSample(const MatPoly &s);
void getPublicEncryptions(
    size_t g,
    vector<MatPoly> &X_v,
    vector<MatPoly> &Y_v,
    vector<MatPoly> &W_v,
    vector<MatPoly> &W_exp_v,
    size_t m_conv,
    size_t m_exp,
    bool for_composition = false,
    bool for_exp = false,
    size_t stopround = 0
);
MatPoly getRegevPublicKeyMatrix(const MatPoly &s, size_t m);