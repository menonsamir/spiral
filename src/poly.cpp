// Implements structures and routines for manipulating polynomials 

#include "poly.h"

inline uint64_t cpu_add_u64(uint64_t operand1, uint64_t operand2,
                                uint64_t *result) {
    *result = operand1 + operand2;
    return (*result < operand1) ? 0x01 : 0x00; // detect overflow
}

inline uint64_t barrett_raw_u128(__uint128_t val, uint64_t const_ratio_0, uint64_t const_ratio_1, uint64_t modulus) {
    uint64_t zx = val & (((__uint128_t)1 << 64) - 1);
    uint64_t zy = val >> 64;

    uint64_t tmp1, tmp3, carry;
    ulonglong2_h prod = umul64wide(zx, const_ratio_0);
    carry = prod.y;
    ulonglong2_h tmp2 = umul64wide(zx, const_ratio_1);
    tmp3 = tmp2.y + cpu_add_u64(tmp2.x, carry, &tmp1);
    tmp2 = umul64wide(zy, const_ratio_0);
    carry = tmp2.y + cpu_add_u64(tmp1, tmp2.x, &tmp1);
    tmp1 = zy * const_ratio_1 + tmp3 + carry;
    tmp3 = zx - tmp1 * modulus;

    return tmp3;
}

inline uint64_t barrett_reduction_u128(__uint128_t val) {
    uint64_t reduced_val = barrett_raw_u128(val, cr0_Q, cr1_Q, Q_i);
    reduced_val -= (Q_i)*(static_cast<int>(reduced_val >= Q_i));
    return reduced_val;
}

void multiply(MatPoly& out, const MatPoly& a, const MatPoly& b) {
    assert(a.cols == b.rows);
    assert(out.rows == a.rows);
    assert(out.cols == b.cols);
    assert(a.isNTT);
    assert(b.isNTT);
    assert(out.isNTT);

    // naive_product_crt_ntt(c.data, a.data, b.data, a.rows, a.cols, b.cols);

    size_t rs = a.rows;
    size_t ms = a.cols;
    size_t cs = b.cols;

    for (size_t r = 0; r < rs; r++) {
        for (size_t c = 0; c < cs; c++) {
            uint64_t *acc = &out.data[(r * cs + c) * crt_count * poly_len];
            for (size_t z = 0; z < crt_count * poly_len; z++) {
                acc[z] = 0;
            }

            for (size_t m = 0; m < ms; m++) {
                for (size_t z = 0; z < poly_len; z++) {
                    // acc += x[r][m] * y[m][c];
                    for (size_t n = 0; n < crt_count; n++) {
                        const uint64_t *x_rm = &a.data[(r * ms + m) * crt_count * poly_len + n * poly_len];
                        const uint64_t *y_mc = &b.data[(m * cs + c) * crt_count * poly_len + n * poly_len];
                        #ifdef MODULUS_56_BIT
                            acc[n*poly_len + z] = (acc[n*poly_len + z] + x_rm[z] * (uint64_t)y_mc[z]);// % moduli[n];//  % moduli[n]; WARN: assumes little overflow
                        #else
                            acc[n*poly_len + z] += (x_rm[z] * (uint64_t)y_mc[z]) % moduli[n];
                        #endif
                    }
                }
            }
            
            #pragma unroll
            for (size_t n = 0; n < crt_count; n++) {
                for (size_t z = 0; z < poly_len; z++) {
                    acc[n*poly_len + z] %= moduli[n];
                }
            }
        }
    }
}

void multiply_no_reduce(MatPoly& out, const MatPoly& a, const MatPoly& b) {
    assert(a.cols == b.rows);
    assert(out.rows == a.rows);
    assert(out.cols == b.cols);
    assert(a.isNTT);
    assert(b.isNTT);
    assert(out.isNTT);

    // naive_product_crt_ntt(c.data, a.data, b.data, a.rows, a.cols, b.cols);

    size_t rs = a.rows;
    size_t ms = a.cols;
    size_t cs = b.cols;

    for (size_t r = 0; r < rs; r++) {
        for (size_t c = 0; c < cs; c++) {
            uint64_t *acc = &out.data[(r * cs + c) * crt_count * poly_len];
            for (size_t z = 0; z < crt_count * poly_len; z++) {
                acc[z] = 0;
            }

            for (size_t m = 0; m < ms; m++) {
                for (size_t z = 0; z < poly_len; z++) {
                    // acc += x[r][m] * y[m][c];
                    for (size_t n = 0; n < crt_count; n++) {
                        const uint64_t *x_rm = &a.data[(r * ms + m) * crt_count * poly_len + n * poly_len];
                        const uint64_t *y_mc = &b.data[(m * cs + c) * crt_count * poly_len + n * poly_len];
                        #ifdef MODULUS_56_BIT
                            acc[n*poly_len + z] += (x_rm[z] * (uint64_t)y_mc[z]);
                        #else
                            acc[n*poly_len + z] += (x_rm[z] * (uint64_t)y_mc[z]) % moduli[n];
                        #endif
                    }
                }
            }
        }
    }
}

MatPoly multiply(const MatPoly& a, const MatPoly& b) {
    const MatPoly& ap = a.isNTT ? a : to_ntt(a);
    const MatPoly& bp = b.isNTT ? b : to_ntt(b);
    
    MatPoly c(a.rows, b.cols);
    multiply(c, ap, bp);
    return c;
}

inline uint64_t add_coeff(uint64_t a, uint64_t b, size_t n) {
    //return (a+b) % moduli[n];
    return barrett_coeff(a+b, n);
}

inline void mul_coeff(uint64_t& c, uint64_t a, uint64_t b, size_t n) {
    // c = (a*b) % moduli[n];
    c = barrett_coeff(a*b, n);
}

void add(MatPoly& c, const MatPoly& a, const MatPoly& b) {
    assert(a.rows == b.rows);
    assert(a.cols == b.cols);
    assert(a.rows == c.rows);
    assert(a.cols == c.cols);
    assert(a.isNTT);
    assert(b.isNTT);
    assert(c.isNTT);

    for (size_t i = 0; i < a.rows * a.cols; i++) {
        for (size_t n = 0; n < crt_count; n++) {
            for (size_t z = 0; z < coeff_count; z++) {
                size_t idx = i * crt_count * coeff_count + n * coeff_count + z;
                c.data[idx] = add_coeff(a.data[idx], b.data[idx], n);
            }
        }
    }
}

void add_into(MatPoly& out, const MatPoly& a, const MatPoly& b, size_t t_row, size_t t_col) {
    assert(a.isNTT);
    assert(b.isNTT);
    assert(out.isNTT);
    assert(t_row >= 0);
    assert(t_col >= 0);
    assert(a.rows == out.rows);
    assert(a.cols == out.cols);
    assert(b.rows + t_row <= a.rows);
    assert(b.cols + t_col <= a.cols);

    for (size_t r = 0; r < b.rows; r++) {
        for (size_t c = 0; c < b.cols; c++) {
            for (size_t n = 0; n < crt_count; n++) {
                for (size_t z = 0; z < coeff_count; z++) {
                    size_t r_t = r + t_row;
                    size_t c_t = c + t_col;
                    size_t idx_a = (r_t * a.cols + c_t) * crt_count * coeff_count + n * coeff_count + z;
                    size_t idx_b = (r * b.cols + c) * crt_count * coeff_count + n * coeff_count + z;
                    out.data[idx_a] = add_coeff(a.data[idx_a], b.data[idx_b], n);
                }
            }
        }
    }
}


MatPoly add(const MatPoly& a, const MatPoly& b) {
    MatPoly out(a.rows, a.cols);
    add(out, a, b);
    return out;
}

void mul_by_const(MatPoly& out, const MatPoly& singlePoly, const MatPoly& a) {
    assert(a.rows == out.rows);
    assert(a.cols == out.cols);
    assert(a.isNTT);
    assert(out.isNTT);
    assert(singlePoly.isNTT);
    assert(singlePoly.rows == 1);
    assert(singlePoly.cols == 1);

    for (size_t r = 0; r < a.rows; r++) {
        for (size_t c = 0; c < a.cols; c++) {
            for (size_t n = 0; n < crt_count; n++) {
                for (size_t z = 0; z < coeff_count; z++) {
                    size_t idx = r * a.cols * crt_count * coeff_count
                        + c * crt_count * coeff_count
                        + n * coeff_count
                        + z;

                    mul_coeff(out.data[idx], a.data[idx], singlePoly.data[n * coeff_count + z], n);
                }
            }
        }
    }
}

MatPoly mul_by_const(const MatPoly& singlePoly, const MatPoly& a) {
    MatPoly out(a.rows, a.cols);
    mul_by_const(out, singlePoly, a);
    return out;
}

void vertical_merge(MatPoly& c, const MatPoly& a, const MatPoly& b) {
    assert(a.cols == b.cols);
    assert(a.rows + b.rows == c.rows);
    assert(a.cols == c.cols);
    assert(a.isNTT == c.isNTT);
    assert(a.isNTT == b.isNTT);
    size_t factor = (a.isNTT) ? crt_count : 1; 

    size_t a_coeffs = a.rows * a.cols * factor * coeff_count;

    for (size_t i = 0; i < a_coeffs; i++) {
        c.data[i] = a.data[i];
    }

    for (size_t i = 0; i < b.rows * b.cols * factor * coeff_count; i++) {
        c.data[a_coeffs + i] = b.data[i];
    }
}

void automorph(MatPoly& out, const MatPoly& a, uint64_t t) {
    assert(a.rows == out.rows);
    assert(a.cols == out.cols);
    assert(!a.isNTT);
    assert(!out.isNTT);

    for (size_t r = 0; r < a.rows; r++) {
        for (size_t c = 0; c < a.cols; c++) {
            size_t inp_idx = r * a.cols * coeff_count + c * coeff_count;
            for (size_t i = 0; i < coeff_count; i++) {
                uint64_t num = ((i * t) / coeff_count);
                uint64_t rem = ((i * t) % coeff_count);

                if (num % 2 == 0) {
                    out.data[inp_idx + rem] = a.data[inp_idx + i];
                } else {
                    out.data[inp_idx + rem] = Q_i - a.data[inp_idx + i];
                }
            }
        }
    }
}

MatPoly automorph(const MatPoly& a, uint64_t t) {
    MatPoly out(a.rows, a.cols, false);
    automorph(out, a, t);
    return out;
}

void invert(MatPoly& out, const MatPoly& a) {
    assert(a.rows == out.rows);
    assert(a.cols == out.cols);
    assert(!a.isNTT);
    assert(!out.isNTT);

    for (size_t r = 0; r < a.rows; r++) {
        for (size_t c = 0; c < a.cols; c++) {
            size_t inp_idx = r * a.cols * coeff_count + c * coeff_count;
            for (size_t i = 0; i < coeff_count; i++) {
                out.data[inp_idx + i] = Q_i_u128 - a.data[inp_idx + i];
            }
        }
    }
}

MatPoly invert(const MatPoly& a) {
    MatPoly out(a.rows, a.cols, false);
    invert(out, a);
    return out;
}

void to_ntt_no_reduce(MatPoly& out, const MatPoly& a) {
    assert(a.rows == out.rows);
    assert(a.cols == out.cols);
    assert(!a.isNTT);
    assert(out.isNTT);

    for (size_t r = 0; r < a.rows; r++) {
        for (size_t c = 0; c < a.cols; c++) {
            size_t inp_idx = r * a.cols * coeff_count + c * coeff_count;
            for (size_t n = 0; n < crt_count; n++) {
                for (size_t z = 0; z < coeff_count; z++) {
                    uint64_t val = a.data[inp_idx + z];
                    out.data[crt_count * inp_idx + coeff_count * n + z] = val;
                }
            }
            ntt_forward(&out.data[crt_count * inp_idx]);
        }
    }
}

void to_ntt(MatPoly& out, const MatPoly& a) {
    assert(a.rows == out.rows);
    assert(a.cols == out.cols);
    assert(!a.isNTT);
    assert(out.isNTT);

    for (size_t r = 0; r < a.rows; r++) {
        for (size_t c = 0; c < a.cols; c++) {
            size_t inp_idx = r * a.cols * coeff_count + c * coeff_count;
            for (size_t n = 0; n < crt_count; n++) {
                for (size_t z = 0; z < coeff_count; z++) {
                    uint64_t val = a.data[inp_idx + z];
                    out.data[crt_count * inp_idx + coeff_count * n + z] = barrett_coeff(val, n);//val % moduli[n];
                }
            }
            ntt_forward(&out.data[crt_count * inp_idx]);
        }
    }
}

MatPoly to_ntt(const MatPoly& a) {
    MatPoly res(a.rows, a.cols, true);
    to_ntt(res, a);
    return res;
}

uint64_t crt_compose_pa(uint64_t x, uint64_t y) {
    // uint64_t val_ap = x * a_inv_p_i;
    // val_ap += y * p_inv_a_i;
    // return val_ap % q_intermediate;
    return x;
}

uint64_t crt_compose(uint64_t x, uint64_t y, uint64_t z) {
    // uint64_t val_ap = x * a_inv_p_i;
    // val_ap += y * p_inv_a_i;
    uint64_t val_ap_u64 = x;
    __uint128_t val = val_ap_u64 * b_inv_pa_i;
    val += y * pa_inv_b_i;
    // cout << "(" << x << ", " << y << ", " << z << ") -> " << (uint64_t)(val % Q_i) << endl;
    // return val % Q_i;
    return barrett_reduction_u128(val);
}

uint64_t *scratch;

void from_ntt(MatPoly& out, const MatPoly& a) {
    assert(a.rows == out.rows);
    assert(a.cols == out.cols);
    assert(a.isNTT);
    assert(!out.isNTT);

    for (size_t r = 0; r < a.rows; r++) {
        for (size_t c = 0; c < a.cols; c++) {
            size_t inp_idx = r * a.cols * coeff_count + c * coeff_count;
            memcpy(scratch, &a.data[crt_count * inp_idx], pol_bytes);
            ntt_inverse(scratch);

            for (size_t z = 0; z < coeff_count; z++) {
                uint64_t v_x = scratch[z];
                uint64_t v_y = scratch[coeff_count + z];
                // uint64_t v_z = scratch[2*coeff_count + z];
                uint64_t val = crt_compose(v_x, v_y, 0);
                out.data[inp_idx + z] = val;
            }
        }
    }
}

MatPoly from_ntt(const MatPoly& a) {
    MatPoly res(a.rows, a.cols, false);
    from_ntt(res, a);
    return res;
}

MatPoly single_poly(uint64_t value) {
    MatPoly res(1, 1, false);
    res.data[0] = value;
    return res;
}

void to_simple_crtd(uint64_t *out, const MatPoly& a) {
    assert(!a.isNTT);

    for (size_t i = 0; i < a.rows * a.cols * coeff_count; i++) {
        out[i] = a.data[i];
    }
}

void reduce_mod(MatPoly &a, uint64_t mod) {
    assert(!a.isNTT);

    for (size_t i = 0; i < a.rows * a.cols * coeff_count; i++) {
        a.data[i] = a.data[i] % mod;
    }
}

void cop(MatPoly& out, const MatPoly& a) {
    assert(out.isNTT == a.isNTT);
    assert(out.rows <= a.rows);
    assert(out.cols <= a.cols);

    size_t factor = (out.isNTT) ? crt_count : 1;

    for (size_t r = 0; r < out.rows; r++) {
        for (size_t c = 0; c < out.cols; c++) {
            memcpy(
                &out.data[(r * out.cols + c) * factor * coeff_count],
                &a.data[(r * a.cols + c) * factor * coeff_count],
                factor * coeff_count * sizeof(uint64_t)
            );
        }
    }
}


void cop(
    MatPoly& out, const MatPoly& a, 
    size_t s_row, size_t s_col, 
    size_t t_row, size_t t_col,
    size_t num_row, size_t num_col
) {
    assert(out.isNTT == a.isNTT);
    assert(t_row + num_row <= out.rows);
    assert(t_col + num_col <= out.cols);
    assert(s_row + num_row <= a.rows);
    assert(s_col + num_col <= a.cols);

    size_t factor = (out.isNTT) ? crt_count : 1;

    for (size_t r = 0; r < num_row; r++) {
        for (size_t c = 0; c < num_col; c++) {
            memcpy(
                &out.data[((r + t_row) * out.cols + (c + t_col)) * factor * coeff_count],
                &a.data[((r + s_row) * a.cols + (c + s_col)) * factor * coeff_count],
                factor * coeff_count * sizeof(uint64_t)
            );
        }
    }
}

void place(MatPoly& out, const MatPoly& a, size_t t_row, size_t t_col) {
    assert(out.isNTT == a.isNTT);
    assert(t_row >= 0);
    assert(t_col >= 0);
    assert(a.rows + t_row <= out.rows);
    assert(a.cols + t_col <= out.cols);

    size_t factor = (out.isNTT) ? crt_count : 1;

    for (size_t r = 0; r < a.rows; r++) {
        for (size_t c = 0; c < a.cols; c++) {
            memcpy(
                &out.data[((r + t_row) * out.cols + (c + t_col)) * factor * coeff_count],
                &a.data[(r * a.cols + c) * factor * coeff_count],
                factor * coeff_count * sizeof(uint64_t)
            );
        }
    }
}

void pick(MatPoly& out, const MatPoly& a, size_t t_row, size_t t_col) {
    assert(out.isNTT == a.isNTT);
    assert(t_row >= 0);
    assert(t_col >= 0);
    assert(t_row < a.rows);
    assert(t_col < a.cols);

    size_t factor = (out.isNTT) ? crt_count : 1;

    for (size_t r = 0; r < out.rows; r++) {
        for (size_t c = 0; c < out.cols; c++) {
            memcpy(
                &out.data[(r * out.cols + c) * factor * coeff_count],
                &a.data[((r + t_row) * a.cols + (c + t_col)) * factor * coeff_count],
                factor * coeff_count * sizeof(uint64_t)
            );
        }
    }
}

MatPoly pick(const MatPoly& a, size_t t_row, size_t t_col, size_t num_rows, size_t num_cols) {
    MatPoly out(num_rows, num_cols, false);
    pick(out, a, t_row, t_col);
    return out;
}

bool is_eq(const MatPoly &A, const MatPoly &B) {
    assert(!A.isNTT);
    assert(!B.isNTT);
    assert(A.rows == B.rows);
    assert(A.cols == B.cols);

    // for (size_t i = 0; i < A.rows * A.cols * crt_count * coeff_count; i++) {
    //     if (A.data[i] != B.data[i]) {
    //         cout << i << " " << A.data[i] << ", " << B.data[i] << endl;
    //         exit(1);
    //     }
    // }

    return memcmp(
        A.data,
        B.data,
        A.rows * A.cols * coeff_count * sizeof(uint64_t)
    ) == 0;
}

// void serialize(const MatPoly &mat, fstream &f) {
//     assert(!mat.isNTT);
//     for (size_t r = 0; r < mat.rows; r++) {
//         for (size_t c = 0; c < mat.cols; c++) {
//             for (size_t m = 0; m < coeff_count; m++) {
//                 __uint128_t val = ((__uint128_t*)mat.data)[(r * mat.cols + c) * coeff_count + m];
//                 f.write((char *)&val, 16);
//             }
//         }
//     }
// }

// void deserialize(MatPoly &mat, size_t rs, size_t cs, fstream &f) {
//     mat = MatPoly(rs, cs, false);
//     uint64_t val[2];
//     for (size_t r = 0; r < rs; r++) {
//         for (size_t c = 0; c < cs; c++) {
//             for (size_t m = 0; m < coeff_count; m++) {
//                 f.read((char *)val, 16);
//                 ((__uint128_t*)mat.data)[(r * cs + c) * coeff_count + m] = *((__uint128_t*)&val);
//             }
//         }
//     }
// }

void build_from_constants(MatPoly &mat, vector<vector<uint64_t> > inp) {
    assert(!mat.isNTT);

    for (size_t r = 0; r < mat.rows; r++) {
        for (size_t c = 0; c < mat.cols; c++) {
            mat[r][c * coeff_count] = inp[r][c];
        }
    }
}

void decode_mat(MatPoly &mat, uint64_t *inp) {
    assert(mat.isNTT);

    memcpy(
        mat.data,
        inp,
        mat.rows * mat.cols * crt_count * coeff_count * sizeof(uint64_t)
    );
}

void divide_by_const(MatPoly& out, const MatPoly& a, uint64_t divisor, uint64_t mod) {
    assert(!a.isNTT);
    assert(!out.isNTT);

    for (size_t i = 0; i < a.rows * a.cols * coeff_count; i++) {
        int64_t val = a.data[i] % Q_i;
        // if (divisor == (1<<20)) cout << "!!!!!! " << val << endl;
        // if (val > (Q_i_u128/2UL)) val -= (int64_t)Q_i_u128;
        //assert(val > 0 && val < (q_const * (1<<12)));
        // cout << val << " / " << divisor << " = " << ((int64_t)round(((long double)val) / ((long double)(divisor)))) << endl;
        int64_t scaled_val = (int64_t)round(((long double)val) / ((long double)(divisor)));
        out.data[i] = (scaled_val + mod) % mod;
    }
}

uint64_t rescale(uint64_t a, uint64_t inp_mod, uint64_t out_mod, bool nosign) {
    int64_t inp_val = a % inp_mod;
    if (inp_val >= (inp_mod / 2)) {
        inp_val -= inp_mod;
    }
    int64_t sign = inp_val >= 0 ? 1 : -1;
    // if (nosign) sign = 1;
    __int128_t val = inp_val * (__int128_t)out_mod;
    __int128_t result = (val + sign*((int64_t)inp_mod/2)) / (__int128_t)inp_mod;
    // assert((result > -inp_mod/2) && (result < inp_mod/2));
    result = (result + (inp_mod/out_mod)*out_mod + 2*out_mod) % out_mod;
    assert(result >= 0);
    return (result + out_mod) % out_mod;
}

MatPoly getRescaled(const MatPoly &a, uint64_t inp_mod, uint64_t out_mod) {
    MatPoly b = a;
    for (size_t i = 0; i < a.rows * a.cols * poly_len; i++) {
        uint64_t inp_val = a.data[i] % Q_i;
        uint64_t result = rescale(inp_val, inp_mod, out_mod, true);
        b.data[i] = result;
    }
    return b;
}

void divide_by_const(MatPoly& out, const MatPoly& a, uint64_t numer, uint64_t denom, uint64_t mod) {
    assert(!a.isNTT);
    assert(!out.isNTT);

    for (size_t i = 0; i < a.rows * a.cols * coeff_count; i++) {
        // int64_t inp_val = a.data[i] % denom;
        // if (inp_val > denom / 2) {
        //     inp_val -= denom;
        // }
        // __int128_t val = inp_val * (__uint128_t)numer;
        // int64_t sign = val >= 0 ? 1 : -1;
        // int64_t result = ((val + sign*denom/2) / denom) % mod;
        // result = (result + (denom/mod)*mod + 2*mod) % mod;

        // // int64_t result = (int64_t)round(((long double)inp_val) * ((long double)numer) / ((long double)(denom)));
        // // __uint128_t val = inp_val * (__uint128_t)numer;
        // // uint64_t result = (val + denom/2) / denom;
        // // if (scaled_val < 0) scaled_val += mod;
        // out.data[i] = (result + mod) % mod;
        out.data[i] = rescale(a.data[i], denom, numer);
    }
}