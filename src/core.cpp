#include "core.h"

// constexpr uint64_t q_i = q_const;
// constexpr uint64_t qq_i = qq_const;

/*
(2[inv or for], 2[q or qq], 2[normal or scaled], poly_len)
above tables are stuctured as:              [each row is poly_len 64-bit vals]
    q  inv_root_powers_div_two              [for inv ntt]
    q  scaled_inv_root_powers_div_two
    qq inv_root_powers_div_two
    qq scaled_inv_root_powers_div_two
    q  root_powers                          [for forward ntt]
    q  scaled_root_powers
    qq root_powers
    qq scaled_root_powers
*/

// num_bits must be less than 64
uint64_t read_arbitrary_bits(const uint64_t *p, size_t bit_offs, size_t num_bits) {
    size_t word_off = bit_offs / 64;
    size_t bit_off_within_word = bit_offs % 64;
    if ((bit_off_within_word + num_bits) <= 64) {
        uint64_t val = *(uint64_t *)(&((uint64_t *)p)[word_off]);
        return (val >> bit_off_within_word) & ((1UL << num_bits) - 1);
    } else {
        __uint128_t val = *(__uint128_t *)(&((uint64_t *)p)[word_off]);
        return (val >> bit_off_within_word) & ((1ULL << num_bits) - 1);
    }
}
// num_bits must be less than 64
void write_arbitrary_bits(uint64_t *p, uint64_t val, size_t bit_offs, size_t num_bits) {
    size_t word_off = bit_offs / 64;
    size_t bit_off_within_word = bit_offs % 64;
    val = val & ((1ULL << num_bits) - 1);
    if ((bit_off_within_word + num_bits) <= 64) {
        uint64_t *curVal = (uint64_t *)(&((uint64_t *)p)[word_off]);
        *curVal &= ~(((1L << num_bits) - 1) << bit_off_within_word);
        *curVal |= (val << bit_off_within_word);
    } else {
        __uint128_t *curVal = (__uint128_t *)(&((uint64_t *)p)[word_off]);
        // cout << bitset<64>(*curVal) << " " << endl;
        __uint128_t mask = ~(((1ULL << num_bits) - 1) << bit_off_within_word);
        // cout << bitset<64>(mask) << " &" << endl;
        // cout << "----------------------------------------------------------------" << endl;
        *curVal &= mask;
        // cout << bitset<64>(*curVal) << " w> " << val << endl;
        *curVal |= ((__uint128_t)val << bit_off_within_word);
        // cout << bitset<64>(*curVal) << endl;
        // cout << endl;
    }
}

void naive_product_crt_ntt(
    uint64_t *out, const uint64_t *x, const uint64_t *y,
    size_t rs, size_t ms, size_t cs
) {
    for (size_t r = 0; r < rs; r++) {
        for (size_t c = 0; c < cs; c++) {
            uint64_t *acc = &out[(r * cs + c) * crt_count * poly_len];
            for (size_t z = 0; z < poly_len; z++) {
                acc[z] = 0;
                acc[poly_len + z] = 0;
                acc[2*poly_len + z] = 0;
            }

            for (size_t m = 0; m < ms; m++) {
                for (size_t z = 0; z < poly_len; z++) {
                    // acc += x[r][m] * y[m][c];
                    for (size_t n = 0; n < crt_count; n++) {
                        const uint64_t *x_rm = &x[(r * ms + m) * crt_count * poly_len + n * poly_len];
                        const uint64_t *y_mc = &y[(m * cs + c) * crt_count * poly_len + n * poly_len];
                        acc[n*poly_len + z] += (x_rm[z] * (uint64_t)y_mc[z]);// % moduli[n];// % q_i;
                    }
                }
            }
            
            for (size_t n = 0; n < crt_count; n++) {
                for (size_t z = 0; z < poly_len; z++) {
                    acc[n*poly_len + z] %= moduli[n];
                }
            }
        }
    }
}

#ifdef __EMSCRIPTEN__
const size_t SIZE_OF_RANDOM_BUF = 65536;
uint8_t *random_buf = NULL;
uint8_t *random_buf_cur = NULL;

static int randombytes_emscripten(void *buf, size_t n)
{
    const int ret = EM_ASM_INT({
        var mybuf = new window.Uint8Array($1);
        window.crypto.getRandomValues(mybuf);
        writeArrayToMemory(mybuf, $0);
        return 0;
    }, buf, n);
    return ret;
}

/*
void uniform_matrix(Mat<ZZ_pX> &A)
{
    long rs = A.NumRows();
    long cs = A.NumCols();

    fill_zero(A);

    uint8_t *random_bytes = (uint8_t *)malloc(poly_len * sizeof(uint64_t));
    
    for (size_t r = 0; r < rs; r++)
    {
        for (size_t c = 0; c < cs; c++)
        {
            randombytes_emscripten(random_bytes, poly_len * sizeof(uint64_t));
            
            for (size_t z = 0; z < poly_len; z++)
            {
                ZZ val = ZZFromBytes(&random_bytes[z * sizeof(uint64_t)], sizeof(uint64_t));
                A[r][c][z] = ZZ_p(INIT_VAL, val);
            }
            // random(A[r][c], 4096);
        }
    }

    size_right(A);

    free(random_bytes);
}
*/

class random_device_doubles
{
    public:
    uint32_t operator()()
    {
        if (random_buf == NULL)
        {
            random_buf = (uint8_t *)malloc(SIZE_OF_RANDOM_BUF);
            randombytes_emscripten(random_buf, SIZE_OF_RANDOM_BUF);
            random_buf_cur = random_buf;
        }
        if ((random_buf + SIZE_OF_RANDOM_BUF) - random_buf_cur < sizeof(uint32_t))
        {
            randombytes_emscripten(random_buf, SIZE_OF_RANDOM_BUF);
            random_buf_cur = random_buf;
        }
        uint32_t a = 0;
        memcpy(&a, random_buf_cur, sizeof(uint32_t));
        random_buf_cur += sizeof(uint32_t);
        return a;
    }

    typedef uint32_t result_type;

    static constexpr uint32_t min()
    {
        return 0u;
    }

    static constexpr uint32_t max()
    {
        return 0xFFFFFFFFu;
    }
};

extern "C" {
    int32_t sample_fast()
    {
        return (int32_t)sample();
    }
}

random_device_doubles random_device_dist;
#else

std::random_device random_device_dist;
#endif

constexpr int64_t num_widths = 10;
constexpr double width = 6.4;
constexpr double std_dev = 2.55323059;
static int64_t max_val = 0;
static vector<double> table;
discrete_distribution<int64_t> distribution;

void build_table() {
    max_val = (int64_t)ceil(width * num_widths);
    for (int64_t i = -max_val; i < max_val+1; i++) {
        double p_val = exp(-M_PI * pow(i, 2) / pow(width, 2));
        table.push_back(p_val);
    }

    distribution = discrete_distribution<int64_t>(table.begin(), table.end());
}

// std::normal_distribution<double> distribution(0, 16.042420957638403); // variance > 10

// mt19937_64 gen(7);
mt19937_64 gen(random_device_dist()); // NOT SECURE

int64_t sample() {
    int64_t val = (distribution(gen)) - max_val;
    return val;
}

// q  = 4398046511039           (~2^42)
// q' = 1152921504606846719     (~2^60)

inline unsigned char add_uint64(uint64_t operand1, uint64_t operand2,
                                uint64_t *result) {
    *result = operand1 + operand2;
    return static_cast<unsigned char>(*result < operand1);
}

inline void multiply_uint64_hw64(uint64_t operand1, uint64_t operand2,
                                 uint64_t *hw64) {
    // auto operand1_coeff_right = operand1 & 0x00000000FFFFFFFF;
    // auto operand2_coeff_right = operand2 & 0x00000000FFFFFFFF;
    // operand1 >>= 32;
    // operand2 >>= 32;

    // auto middle1 = operand1 * operand2_coeff_right;
    // uint64_t middle;
    // auto left = operand1 * operand2 +
    //             (static_cast<uint64_t>(add_uint64(
    //                  middle1, operand2 * operand1_coeff_right, &middle))
    //              << 32);
    // auto right = operand1_coeff_right * operand2_coeff_right;
    // auto temp_sum = (right >> 32) + (middle & 0x00000000FFFFFFFF);

    // *hw64 = static_cast<unsigned long long>(left + (middle >> 32) +
    //                                         (temp_sum >> 32));
    *hw64 = (operand1 * (unsigned __int128)operand2) >> 64;
}
// expects operand: 2 x poly_len uint64_ts for coefficients in CRT (mod q and then
// mod qq)

intel::hexl::NTT *ntt_p;
intel::hexl::NTT *ntt_a;
intel::hexl::NTT *ntt_b;
intel::hexl::NTT *ntt_qprime;

void ntt_forward_old(uint64_t *operand_overall);
void ntt_forward(uint64_t *operand_overall) {
    // ntt_p->ComputeForward(operand_overall, operand_overall, 1, 1);
    // ntt_a->ComputeForward(&operand_overall[poly_len], &operand_overall[poly_len], 1, 1);
    ntt_forward_old(operand_overall);
    // ntt_b->ComputeForward(&operand_overall[poly_len], &operand_overall[poly_len], 1, 1);
}

void ntt_forward_old(uint64_t *operand_overall) {
    for (size_t coeff_mod = 0; coeff_mod < 2; coeff_mod++) {
        size_t n = poly_len;

        const uint64_t *forward_tables = &tables[poly_len * 2 * 2 + coeff_mod * poly_len * 2];
        uint64_t *operand = &operand_overall[coeff_mod * poly_len];
        uint32_t modulus_small = coeff_mod == 0 ? p_i : b_i;
        uint32_t two_times_modulus_small = 2 * modulus_small;
        
        for (size_t mm = 0; mm < coeff_count_pow_of_2; mm++) {
            size_t m = 1 << mm;
            size_t t = n >> (mm+1);
            // m*t is always 1024

            for (size_t i = 0; i < m; i++) {
                const uint64_t W = forward_tables[m + i];
                const uint64_t Wprime = forward_tables[poly_len + m + i];

                #ifdef USE_AVX2
                if (t < 4) {
                #endif

                    for (size_t j = 0; j < t; j++) {
                        // The Harvey butterfly: assume X, Y in [0, 2p), and return
                        // X', Y' in [0, 4p). X', Y' = X + WY, X - WY (mod p).
                        
                        uint64_t *p_x = &operand[2 * i * t + j];
                        uint64_t *p_y = &operand[2 * i * t + t + j];
                        uint32_t x = *p_x;
                        uint32_t y = *p_y;

                        uint32_t currX = x - (two_times_modulus_small * static_cast<int>(x >= two_times_modulus_small));

                        uint64_t Q = ((y) * (uint64_t)(Wprime)) >> 32;

                        Q = W * y - Q * modulus_small;
                        *p_x = currX + Q;
                        *p_y = currX + (two_times_modulus_small - Q);
                    }
                
                #ifdef USE_AVX2
                } else { // t >= 4
                    uint64_t *p_base = &operand[2 * i * t];
                    for (size_t j = 0; j < t; j+= 4) {
                        // The Harvey butterfly: assume X, Y in [0, 2p), and return
                        // X', Y' in [0, 4p). X', Y' = X + WY, X - WY (mod p).
                        
                        uint64_t *p_x = &p_base[j];
                        uint64_t *p_y = &p_base[t + j];
                        __m256i x = _mm256_loadu_si256((__m256i const *)p_x);
                        __m256i y = _mm256_loadu_si256((__m256i const *)p_y);

                        // uint32_t currX = x - (two_times_modulus_small * static_cast<int>(x >= two_times_modulus_small));
                        __m256i cmp_val = _mm256_set1_epi64x(two_times_modulus_small);
                        __m256i gt_mask = _mm256_cmpgt_epi64(x, cmp_val);
                        __m256i to_subtract = _mm256_and_si256(gt_mask, cmp_val);
                        __m256i currX = _mm256_sub_epi64(x, to_subtract);

                        // uint32_t Q = ((y) * (uint64_t)(Wprime)) >> 32;
                        __m256i Wprime_vec = _mm256_set1_epi64x(Wprime);
                        __m256i product = _mm256_mul_epu32(y, Wprime_vec);
                        __m256i Q = _mm256_srli_epi64(product, 32);

                        // Q = W * y - Q * modulus_small;
                        __m256i W_vec = _mm256_set1_epi64x(W);
                        __m256i WY = _mm256_mul_epu32(y, W_vec);
                        __m256i modulus_small_vec = _mm256_set1_epi64x(modulus_small);
                        __m256i Q_scaled = _mm256_mul_epu32(Q, modulus_small_vec);
                        __m256i Q_final = _mm256_sub_epi64(WY, Q_scaled);

                        __m256i new_x = _mm256_add_epi64(currX, Q_final);
                        __m256i Q_final_inverted = _mm256_sub_epi64(cmp_val, Q_final);
                        __m256i new_y = _mm256_add_epi64(currX, Q_final_inverted);

                        // *p_x = currX + Q;
                        // *p_y = currX + (two_times_modulus_small - Q);
                        _mm256_storeu_si256 ((__m256i *)p_x, new_x);
                        _mm256_storeu_si256 ((__m256i *)p_y, new_y);
                    }
                }
                #endif
            }
        }

        #ifdef USE_AVX2
            for (size_t i = 0; i < poly_len; i+=4) {
                    __m256i cmp_val1 = _mm256_set1_epi64x(two_times_modulus_small);
                    __m256i x = _mm256_loadu_si256((__m256i const *)(&operand[i]));
                    __m256i gt_mask = _mm256_cmpgt_epi64(x, cmp_val1);
                    __m256i to_subtract = _mm256_and_si256(gt_mask, cmp_val1);
                    x = _mm256_sub_epi64(x, to_subtract);

                    __m256i cmp_val2 = _mm256_set1_epi64x(modulus_small);
                    gt_mask = _mm256_cmpgt_epi64(x, cmp_val2);
                    to_subtract = _mm256_and_si256(gt_mask, cmp_val2);
                    x = _mm256_sub_epi64(x, to_subtract);
                    _mm256_storeu_si256 ((__m256i *)(&operand[i]), x);
            }
        #else
            for (size_t i = 0; i < poly_len; i++) {
                operand[i] -= static_cast<int>(operand[i] >= two_times_modulus_small) * two_times_modulus_small;
                operand[i] -= static_cast<int>(operand[i] >= modulus_small) * modulus_small;
            }
        #endif
        
        // forward_tables = &tables[poly_len * 2 * 2 + poly_len * 2];
        // operand = &operand_overall[poly_len];
        // uint64_t modulus = a_i;
        // uint64_t two_times_modulus = 2 * modulus;

        // for (size_t mm = 0; mm < coeff_count_pow_of_2; mm++) {
        //     size_t m = 1 << mm;
        //     size_t t = n >> (mm+1);
        //     // m*t is always 1024

        //     for (size_t i = 0; i < m; i++) {
        //         const uint64_t W = forward_tables[m + i];
        //         const uint64_t Wprime = forward_tables[poly_len + m + i];

        //         for (size_t j = 0; j < t; j++) {
        //             // The Harvey butterfly: assume X, Y in [0, 2p), and return
        //             // X', Y' in [0, 4p). X', Y' = X + WY, X - WY (mod p).
                    
        //             uint64_t *p_x = &operand[2 * i * t + j];
        //             uint64_t *p_y = &operand[2 * i * t + t + j];
        //             uint64_t x = *p_x;
        //             uint64_t y = *p_y;

        //             uint64_t currX = x - (two_times_modulus * static_cast<int>(x >= two_times_modulus));

        //             uint64_t Q = ((y) * (unsigned __int128)(Wprime)) >> 64;

        //             Q = W * y - Q * modulus;
        //             *p_x = currX + Q;
        //             *p_y = currX + (two_times_modulus - Q);
        //         }
        //     }
        // }

        // #if USE_AVX2
        //     for (size_t i = 0; i < poly_len; i+=4) {
        //         // operand[i] -= static_cast<int>(operand[i] >= two_times_modulus) * two_times_modulus;
        //         // operand[i] -= static_cast<int>(operand[i] >= modulus) * modulus;
        //         __m256i cmp_val1 = _mm256_set1_epi64x(two_times_modulus);
        //         __m256i x = _mm256_loadu_si256((__m256i const *)(&operand[i]));
        //         __m256i gt_mask = _mm256_cmpgt_epi64(x, cmp_val1);
        //         __m256i to_subtract = _mm256_and_si256(gt_mask, cmp_val1);
        //         x = _mm256_sub_epi64(x, to_subtract);

        //         __m256i cmp_val2 = _mm256_set1_epi64x(modulus);
        //         gt_mask = _mm256_cmpgt_epi64(x, cmp_val2);
        //         to_subtract = _mm256_and_si256(gt_mask, cmp_val2);
        //         x = _mm256_sub_epi64(x, to_subtract);
        //         _mm256_storeu_si256 ((__m256i *)(&operand[i]), x);
        //     }
        // #else
        //     for (size_t i = 0; i < poly_len; i++) {
        //         operand[i] -= static_cast<int>(operand[i] >= two_times_modulus) * two_times_modulus;
        //         operand[i] -= static_cast<int>(operand[i] >= modulus) * modulus;
        //     }
        // #endif
    }
}

void ntt_inverse_old(uint64_t *operand_overall);
void ntt_inverse(uint64_t *operand_overall) {
    // ntt_p->ComputeInverse(operand_overall, operand_overall, 1, 1);
    // ntt_a->ComputeInverse(&operand_overall[poly_len], &operand_overall[poly_len], 1, 1);
    ntt_inverse_old(operand_overall);
    // ntt_b->ComputeInverse(&operand_overall[poly_len], &operand_overall[poly_len], 1, 1);
}

void ntt_inverse_old(uint64_t *operand_overall) {
    for (size_t coeff_mod = 0; coeff_mod < 2; coeff_mod++) {
        const uint64_t *inverse_tables = &tables[coeff_mod * poly_len * 2];
        uint64_t *operand = &operand_overall[coeff_mod == 0 ? 0 : poly_len];
        uint64_t modulus = coeff_mod == 0 ? p_i : b_i;

        uint64_t two_times_modulus = 2 * modulus;
        size_t n = poly_len;
        size_t t = 1;

        for (size_t m = n; m > 1; m >>= 1) {
            size_t j1 = 0;
            size_t h = m >> 1;
            for (size_t i = 0; i < h; i++) {
                size_t j2 = j1 + t;
                // Need the powers of  phi^{-1} in bit-reversed order
                const uint64_t W = inverse_tables[h + i];
                const uint64_t Wprime = inverse_tables[poly_len + h + i];

                uint64_t *U = operand + j1;
                uint64_t *V = U + t;
                uint64_t currU;
                uint64_t T;
                uint64_t H;
                for (size_t j = j1; j < j2; j++) {
                    // U = x[i], V = x[i+m]

                    // Compute U - V + 2q
                    T = two_times_modulus - *V + *U;

                    // Cleverly check whether currU + currV >= two_times_modulus
                    currU = *U + *V - (two_times_modulus * static_cast<int>((*U << 1) >= T));

                    // Need to make it so that div2_uint_mod takes values that
                    // are > q.
                    // div2_uint_mod(U, modulusptr, coeff_uint64_count, U);
                    // We use also the fact that parity of currU is same as
                    // parity of T. Since our modulus is always so small that
                    // currU + masked_modulus < 2^64, we never need to worry
                    // about wrapping around when adding masked_modulus.
                    // uint64_t masked_modulus = modulus &
                    // static_cast<uint64_t>(-static_cast<int64_t>(T & 1));
                    // uint64_t carry = add_uint64(currU, masked_modulus, 0,
                    // &currU); currU += modulus &
                    // static_cast<uint64_t>(-static_cast<int64_t>(T & 1));
                    *U++ = (currU + (modulus * static_cast<int>(T & 1))) >> 1;

                    //multiply_uint64_hw64(Wprime, T, &H);
                    H = ((T) * (uint64_t)(Wprime)) >> 32;
                    
                    // effectively, the next two multiply perform multiply
                    // modulo beta = 2**wordsize.
                    *V++ = W * T - H * modulus;
                }
                j1 += (t << 1);
            }
            t <<= 1;
        }

        #ifdef USE_AVX2
            for (size_t i = 0; i < poly_len; i+=4) {
                // operand[i] -= static_cast<int>(operand[i] >= two_times_modulus) * two_times_modulus;
                // operand[i] -= static_cast<int>(operand[i] >= modulus) * modulus;

                __m256i cmp_val1 = _mm256_set1_epi64x(two_times_modulus);
                __m256i x = _mm256_loadu_si256((__m256i const *)(&operand[i]));
                __m256i gt_mask = _mm256_cmpgt_epi64(x, cmp_val1);
                __m256i to_subtract = _mm256_and_si256(gt_mask, cmp_val1);
                x = _mm256_sub_epi64(x, to_subtract);

                __m256i cmp_val2 = _mm256_set1_epi64x(modulus);
                gt_mask = _mm256_cmpgt_epi64(x, cmp_val2);
                to_subtract = _mm256_and_si256(gt_mask, cmp_val2);
                x = _mm256_sub_epi64(x, to_subtract);
                _mm256_storeu_si256 ((__m256i *)(&operand[i]), x);

                // operand[i] = operand[i] % modulus;
                // if (operand[i] >= two_times_modulus) operand[i] -=
                // two_times_modulus; if (operand[i] >= modulus) operand[i] -=
                // modulus;
            }
        #else
            for (size_t i = 0; i < poly_len; i++) {
                operand[i] -= static_cast<int>(operand[i] >= two_times_modulus) * two_times_modulus;
                operand[i] -= static_cast<int>(operand[i] >= modulus) * modulus;
            }
        #endif
    }
}

