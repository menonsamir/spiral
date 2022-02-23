#pragma once

#include <inttypes.h>

#define MODULUS_56_BIT

constexpr long coeff_count_pow_of_2 = 11;
constexpr size_t coeff_count = 1 << coeff_count_pow_of_2;
constexpr size_t crt_count = 2;
constexpr size_t pol_bytes = crt_count * coeff_count * sizeof(uint64_t);
constexpr size_t poly_len = 1 << coeff_count_pow_of_2; // 2048;

constexpr uint32_t p_i = 268369921; // 2^28 - 2^16 + 1
constexpr uint32_t a_i = 1;

constexpr uint64_t q_intermediate = p_i * (uint64_t)a_i;
constexpr size_t bits_to_hold_q_intermediate = 28;
#ifdef MODULUS_56_BIT
    // 56-bit coeffs
    constexpr size_t logQ = 56;
    constexpr uint32_t b_i = 249561089UL;//266334209; // 2^28 - 2^21 - 2^12 + 1
    // constexpr uint64_t a_inv_p_i = 10533UL * a_i;                   // (inverse of a mod p) * a
    // constexpr uint64_t p_inv_a_i = 5853UL * p_i;                    // (inverse of p mod a) * p
    constexpr __uint128_t pa_inv_b_i = 97389680UL * q_intermediate;    // (inverse of pa mod b) * pa
    constexpr __uint128_t b_inv_pa_i = 163640210UL * b_i;              // (inverse of b mod pa) * b
    constexpr uint64_t cr0_Q = 7906011006380390721UL;
    constexpr uint64_t cr1_Q = 275UL;
#else
    // 60-bit coeffs
    constexpr size_t logQ = 60;
    constexpr uint32_t b_i = 2147352577;
    constexpr uint64_t a_inv_p_i = 10533UL * a_i;                  // (inverse of a mod p) * a
    constexpr uint64_t p_inv_a_i = 5853UL * p_i;                 // (inverse of p mod a) * p
    constexpr __uint128_t pa_inv_b_i = 1668642485UL * q_intermediate;  // (inverse of pa mod b) * pa
    constexpr __uint128_t b_inv_pa_i = 112216397UL * b_i;            // (inverse of b mod pa) * b
    constexpr uint64_t cr0_Q = 1215693566780587878UL;
    constexpr uint64_t cr1_Q = 17UL;
#endif
constexpr uint32_t moduli[2] = {p_i, b_i};
constexpr uint64_t ab_i = a_i * (uint64_t)b_i;
constexpr uint64_t Q_i = p_i * (uint64_t)a_i * (uint64_t)b_i;
constexpr __uint128_t Q_i_u128 = Q_i;

// Must be set manually
constexpr uint64_t num_bits_q = 13;
constexpr uint64_t bits_to_hold_q = 14;
constexpr uint64_t packed_offset_1 = 32;
constexpr uint64_t packed_offset_diff = 0; // log (a_i)
constexpr uint64_t packed_offset_2 = packed_offset_1 + packed_offset_diff; // + log (p*a_i)

// the max # of elements of bit width (2 log pa OR 2 log b) that can
// summed and still fit in a 64-bit unsigned integer
// must be a power of 2
// log2((pa)^2 * __) <= 64
// log2((b)^2 * __) <= 64
constexpr size_t max_summed_pa_or_b_in_u64 = 1 << 6;

constexpr uint64_t cr0_p = 16144578669088582089UL;
constexpr uint64_t cr1_p = 68736257792UL;
constexpr uint64_t cr0_b = 10966983149909726427UL;
constexpr uint64_t cr1_b = 73916747789UL;

// constexpr uint64_t cr0_b = 9220205491147943927UL;
// constexpr uint64_t cr1_b = 8589606920UL;


constexpr size_t n0 = 2;
constexpr size_t n1 = 3;
constexpr size_t n2 = 2;
constexpr size_t k_param = n1 - n0;
constexpr bool modswitch_on_server = true;
constexpr size_t base_dim = 2;

constexpr uint64_t qprime_mods[37] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12289, 12289, 61441, 65537, 65537, 520193, 786433, 786433, 3604481, 7340033, 16515073, 33292289, 67043329, 132120577, 268369921, 469762049, 1073479681, 2013265921, 4293918721, 8588886017, 17175674881, 34359214081, 68718428161
};

constexpr size_t query_elems_first = QNUMFIRST;
constexpr size_t query_elems_rest = QNUMREST;
constexpr size_t m2 = TGSW * n1;
constexpr size_t t_GSW = TGSW;
constexpr size_t m_conv = TCONV;
constexpr size_t m_exp = TEXP;
constexpr size_t m_exp_right = TEXPRIGHT;
constexpr uint64_t p_db = PVALUE;
constexpr uint64_t arb_qprime = qprime_mods[QPBITS];
constexpr size_t bits_to_hold_arb_qprime = QPBITS;
#ifdef OUTN
constexpr size_t out_n = OUTN;
#else
constexpr size_t out_n = 4;
#endif
constexpr uint64_t scale_k = Q_i / p_db;