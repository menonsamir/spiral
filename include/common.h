#pragma once

// extern uint64_t q_i;
// extern uint64_t qq_i;
#include <immintrin.h>
#include "poly.h"

#if defined(__AVX2__)
    #define USE_AVX2 1
#endif

extern MatPoly S_mp;
extern MatPoly Sp_mp;
extern MatPoly sr_mp;

extern bool nonoise;
extern bool direct_upload;
extern bool ternary;
extern bool random_data;
extern bool show_diff;
extern bool output_err;
extern size_t dummyWorkingSet;
extern size_t max_trials;

extern size_t total_offline_query_size;