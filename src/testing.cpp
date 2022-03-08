#include "testing.h"

static bool silent = false;
static vector<MatPoly> neg1s_mp_here;

stack<std::chrono::time_point<std::chrono::high_resolution_clock> > time_stack;
void st() {
    time_stack.push(chrono::high_resolution_clock::now());
}
uint64_t rec(string msg = "", bool pop = false) {
    auto last = time_stack.top();
    if (pop) time_stack.pop();
    uint64_t duration_us = std::chrono::duration_cast<std::chrono::microseconds>(chrono::high_resolution_clock::now() - last).count();
    if (!silent) cout << msg << "> " << duration_us << endl;
    return duration_us;
}
uint64_t et(string msg = "") {
    return rec(msg, true);
}

void getExpansionKeySwitchingMatrices(
    vector<MatPoly> &W_exp_v,
    size_t g,
    size_t m_exp
) {
    MatPoly G_exp = buildGadget(1, m_exp);
    MatPoly G_exp_nttd = to_ntt(G_exp);

    MatPoly s0 = sr_mp;
    for (size_t i = 0; i < g; i++) {
        size_t t = (poly_len / (1 << i)) + 1;

        MatPoly tau_s0 = automorph(s0, t);

        MatPoly W_exp_i = encryptSimpleRegevMatrix(s0, multiply(tau_s0, G_exp_nttd));
        W_exp_v.push_back(W_exp_i);
    }
}

void coefficientExpansion(
    vector<MatPoly> &cv_v, // first element is base_dim x 1
    size_t g,
    size_t m_exp,
    const vector<MatPoly> &W_left_v,   // g matrices, each base_dim x m_exp
    const vector<MatPoly> &W_right_v,   // g matrices, each base_dim x m_exp
    size_t max_bits_to_gen_right = 0,
    size_t stopround = 0
) {
    assert(cv_v[0].rows == base_dim);
    assert(cv_v[0].cols == 1);

    MatPoly c(base_dim, 1, false);
    MatPoly c_automorphed(base_dim, 1, false);
    MatPoly c_automorphed_0(1, 1, false);
    MatPoly c_automorphed_1(1, 1, false);
    MatPoly c_automorphed_1_nttd(1, 1);
    MatPoly c_automorphed_1_padded(base_dim, 1);
    MatPoly ginv_c(m_exp, 1, false);
    MatPoly ginv_c_nttd(m_exp, 1);
    MatPoly ginv_c_right(m_exp_right, 1, false);
    MatPoly ginv_c_right_nttd(m_exp_right, 1);

    MatPoly W_times_ginv_c(base_dim, 1);

    for (size_t r = 0; r < g; r++) {
        size_t num_in = 1 << r;
        size_t num_out = 2 * num_in;

        size_t t = (poly_len / (1 << r)) + 1;

        const MatPoly &W_left = W_left_v[r];
        const MatPoly &W_right = W_right_v[r];
        const MatPoly &neg1 = neg1s_mp_here[r];
        
        for (size_t i = 0; i < num_out; i++) {
            if (stopround > 0 && r > stopround && (i % 2) == 1) continue;
            if (stopround > 0 && r == stopround && (i % 2) == 1 && i/2 > max_bits_to_gen_right) continue;

            const MatPoly &W        = (i % 2) == 0 ? W_left : W_right;
            MatPoly &gi_c           = (i % 2) == 0 ? ginv_c : ginv_c_right;
            MatPoly &gi_c_nttd      = (i % 2) == 0 ? ginv_c_nttd : ginv_c_right_nttd;
            size_t gadget_dim       = (i % 2) == 0 ? m_exp : m_exp_right;

            if (i < num_in) mul_by_const(cv_v[num_in + i], neg1, cv_v[i]);
            from_ntt(c, cv_v[i]);
            automorph(c_automorphed, c, t);
            pick(c_automorphed_0, c_automorphed, 0, 0);
            pick(c_automorphed_1, c_automorphed, 1, 0);
            to_ntt(c_automorphed_1_nttd, c_automorphed_1);
            place(c_automorphed_1_padded, c_automorphed_1_nttd, 1, 0);
            gadget_invert(gadget_dim, gi_c, c_automorphed_0, 1);
            to_ntt_no_reduce(gi_c_nttd, gi_c);
            multiply(W_times_ginv_c, W, gi_c_nttd);
            size_t idx = 0;
            for (size_t j = 0; j < base_dim; j++) {
                for (size_t n = 0; n < crt_count; n++) {
                    for (size_t z = 0; z < coeff_count; z++) {
                        cv_v[i].data[idx] = barrett_coeff(cv_v[i].data[idx] + W_times_ginv_c.data[idx] + j*c_automorphed_1_nttd.data[n * coeff_count + z], n);
                        idx++;
                    }
                }
            }
        }
    }
}


void regevToSimpleGsw(
    vector<MatPoly> &v_gsw,
    const vector<MatPoly> &v_inp,
    const MatPoly &V,
    size_t m_conv,
    size_t ell,
    size_t further_dims,
    size_t idx_factor,
    size_t idx_offset
) {
    assert(V.rows == base_dim);
    assert(V.cols == base_dim * m_conv);
    MatPoly ginv_c_inp(base_dim * m_conv, 1, false);
    MatPoly ginv_c_inp_ntt(base_dim * m_conv, 1);
    MatPoly tmp_ct_intt(base_dim, 1, false);
    MatPoly tmp_ct(base_dim, 1);

    for (size_t i = 0; i < further_dims; i++) {
        MatPoly ct(base_dim, base_dim * ell);
        for (size_t j = 0; j < ell; j++) {
            size_t idx_ct = i * ell + j;
            size_t idx_inp = idx_factor*(idx_ct) + idx_offset;
            const MatPoly& c_inp = v_inp[idx_inp];
            place(ct, c_inp, 0, base_dim*j + 1);
            from_ntt(tmp_ct_intt, c_inp);
            gadget_invert(base_dim * m_conv, ginv_c_inp, tmp_ct_intt, base_dim);
            to_ntt(ginv_c_inp_ntt, ginv_c_inp);
            multiply(tmp_ct, V, ginv_c_inp_ntt);
            place(ct, tmp_ct, 0, base_dim*j);
        }
        v_gsw.push_back(ct);
    }
}
MatPoly get_fresh_public_key_raw_arb(MatPoly &Sp, long m, size_t out_n) {
    MatPoly A(k_param, m, false);
    uniform_matrix(A);
    MatPoly E(out_n, m, false);
    noise(E);

    MatPoly A_ntt = to_ntt(A);
    MatPoly B_p(out_n, m);
    multiply(B_p, to_ntt(Sp), A_ntt);
    MatPoly B(out_n, m);
    add(B, to_ntt(E), B_p);
    
    MatPoly A_inv_mp_raw(A_ntt.rows, A_ntt.cols, false);
    invert(A_inv_mp_raw, from_ntt(A_ntt));

    MatPoly P(out_n+1, m, false);
    vertical_merge(P, A_inv_mp_raw, from_ntt(B));

    return P;
}


MatPoly encryptMatrixArb(const MatPoly &G, const MatPoly &A, size_t out_n) {
    assert(G.rows == base_dim);
    assert(A.rows == out_n);
    assert(A.cols == base_dim);

    size_t mx = G.cols;

    MatPoly P = get_fresh_public_key_raw_arb(Sp_mp, mx, out_n);

    MatPoly AG = multiply(A, G);

    MatPoly padding(1, mx);
    MatPoly AG_padded(out_n+1, mx);
    vertical_merge(AG_padded, padding, AG);

    MatPoly result(out_n+1, mx);
    add(result, AG_padded, to_ntt(P));

    return result;
}

MatPoly encryptMatrixArbitrary(const MatPoly &AG, size_t out_n) {
    size_t mx = AG.cols;
    MatPoly P = get_fresh_public_key_raw_arb(Sp_mp, mx, out_n);

    MatPoly padding(1, mx);
    MatPoly AG_padded(out_n+1, mx);
    vertical_merge(AG_padded, padding, AG);

    MatPoly result(out_n+1, mx);
    add(result, AG_padded, to_ntt(P));

    return result;
}

void pack(
    MatPoly& result,
    size_t out_n,
    size_t m_conv,
    const vector<MatPoly>& v_ct,
    const vector<MatPoly>& v_W
) {
    assert(result.rows == (out_n + 1));
    assert(result.cols == out_n);
    assert(result.isNTT);
    assert(v_ct.size() >= out_n * out_n);
    assert(v_W.size() == out_n);
    assert(v_ct[0].rows == base_dim);
    assert(v_ct[0].cols == 1);
    assert(!v_ct[0].isNTT);
    assert(v_W[0].rows == (out_n+1));
    assert(v_W[0].cols == m_conv);
    assert(v_W[0].isNTT);

    MatPoly v_int(out_n + 1, 1);
    MatPoly ginv(m_conv, 1, false);
    MatPoly ginv_nttd(m_conv, 1);
    MatPoly prod(out_n + 1, 1);
    MatPoly ct_1(1, 1, false);
    MatPoly ct_2(1, 1, false);
    MatPoly ct_2_ntt(1, 1);

    for (size_t c = 0; c < out_n; c++) {
        v_int = MatPoly(out_n + 1, 1);
        for (size_t r = 0; r < out_n; r++) {
            const MatPoly& W = v_W[r];
            const MatPoly& ct = v_ct[r * out_n + c];
            pick(ct_1, ct, 0, 0);
            pick(ct_2, ct, 1, 0);
            to_ntt(ct_2_ntt, ct_2);
            gadget_invert(m_conv, ginv, ct_1, 1);
            to_ntt(ginv_nttd, ginv);
            multiply(prod, W, ginv_nttd);
            add_into(v_int, v_int, ct_2_ntt, 1 + r, 0);
            add(v_int, v_int, prod);
        }
        place(result, v_int, 0, c);
    }
}

void testPacking(size_t out_n) {
    S_mp = MatPoly(out_n, out_n+1, false);
    Sp_mp = MatPoly(out_n, 1, false);
    sr_mp = MatPoly(1, 1, false);
    keygen(S_mp, Sp_mp, sr_mp, out_n);

    MatPoly result(out_n + 1, out_n);

    MatPoly pt(1, 1, false);
    uniform_matrix(pt);
    reduce_mod(pt, p_db);

    // generate cts
    vector<MatPoly> v_ct;
    for (size_t i = 0; i < out_n*out_n; i++) {
        MatPoly ct = encryptSimpleRegev(pt);
        v_ct.push_back(from_ntt(ct));
    }
    cout << "done w gen" << endl;

    // generate Ws
    size_t m_conv = logQ;
    MatPoly s0 = sr_mp;
    MatPoly s0_1 = getSkVec(0);
    MatPoly G_conv = buildGadget(base_dim, base_dim * m_conv);

    vector<MatPoly> v_W;
    for (size_t i = 0; i < out_n; i++) {
        MatPoly msg(out_n, base_dim, false);
        place(msg, s0_1, i, 0);
        MatPoly W = encryptMatrixArb(G_conv, msg, out_n);
        v_W.push_back(W);
    }
    cout << "done w setup" << endl;

    auto start = chrono::high_resolution_clock::now();
    pack(result, out_n, m_conv, v_ct, v_W);
    double duration = (chrono::high_resolution_clock::now() - start).count() / 1000;
    cout << "> " << duration << endl;    

    MatPoly pt_res = multiply(S_mp, result);
    MatPoly pt_res_elem = pick(from_ntt(pt_res), 0, 0, 1, 1);
    cout << "packing test? " << is_eq(pt_res_elem, pt) << endl;
    // cout << "> ";
    // for (size_t i = 0; i < poly_len; i++) {
    //     cout << pt_res_elem.data[i] << " ";
    // }
    // cout << endl;
    exit(1);
}

uint64_t *alloc_u64(size_t words) {
    return (uint64_t *)aligned_alloc(64, words * sizeof(uint64_t));
}

// TODO: efficiency
void multiplyQueryByDatabaseDim1(
    vector<MatPoly> &out,
    const vector<MatPoly> &db,              // vector of plaintexts
    const vector<MatPoly> &v_firstdim       // 2^nu_1 simple Regev ciphertexts
) {
    MatPoly prod(base_dim, 1);
    size_t dim_first = v_firstdim.size();
    size_t dim_rest = db.size() / dim_first;
    for (size_t i = 0; i < dim_rest; i++) {
        MatPoly &sum = out[i];
        for (size_t j = 0; j < dim_first; j++) {
            multiply(prod, v_firstdim[j], db[j * dim_rest + i]);
            add(sum, sum, prod);
        }
    }
}

uint64_t *convertDb(const vector<MatPoly> &db, size_t dim0, size_t num_per) {
    size_t pt_rows = db[0].rows;
    size_t pt_cols = db[0].cols;
    uint64_t *db_buf = (uint64_t *)malloc(db.size() * pt_rows * pt_cols * poly_len * sizeof(uint64_t));
    // b': i c n z j m
    for (size_t i = 0; i < db.size(); i++) {
        size_t ii = i % num_per;
        size_t j = i / num_per;
        for (size_t m = 0; m < pt_rows; m++) {
            for (size_t c = 0; c < pt_cols; c++) {
                uint64_t *BB = &(db[i].data)[(m * pt_cols + c) * crt_count * coeff_count];
                for (size_t z = 0; z < poly_len; z++) {
                    size_t idx = z  * (num_per * pt_cols * dim0 * pt_rows) +
                                 ii * (pt_cols * dim0 * pt_rows) + 
                                 c  * (dim0 * pt_rows) +
                                 j  * (pt_rows) + 
                                 m;
                    
                    db_buf[idx] = BB[z] | (BB[poly_len + z] << 32);
                }
            }
        }
    }
    return db_buf;
}

void reorientCiphertextsDim1(uint64_t *out, const vector<MatPoly> &v_firstdim, size_t dim0, size_t idx_factor) {
    size_t ct_rows = v_firstdim[0].rows;
    size_t ct_cols = v_firstdim[0].cols;
    for (size_t j = 0; j < dim0; j++) {
        size_t jp = j * idx_factor;
        for (size_t r = 0; r < ct_rows; r++) {
            for (size_t m = 0; m < ct_cols; m++) {
                for (size_t z = 0; z < poly_len; z++) {
                    size_t idx_a_in  = r * (ct_cols*crt_count*poly_len)
                                     + m * (crt_count*poly_len);
                    size_t idx_a_out = z * (dim0*ct_cols*ct_rows)
                                     + j * (ct_cols*ct_rows)
                                     + m * (ct_rows)
                                     + r;

                    out[idx_a_out] = (v_firstdim[jp].data[idx_a_in + z] % p_i) | ((v_firstdim[jp].data[idx_a_in + poly_len + z] % b_i) << 32);
                }
            }
        }
    }
}

void fastMultiplyQueryByDatabaseDim1(
    vector<MatPoly> &out,
    const uint64_t *db,
    const uint64_t *v_firstdim,
    size_t dim0, size_t num_per
) {
    // v_firstdim: 
    //     poly_len, dim0, ct.rows, ct.cols
    // db: 
    //     poly_len, num_per, pt.cols, dim0, pt.rows
    size_t ct_rows = base_dim;
    size_t ct_cols = 1;
    size_t pt_rows = 1;
    size_t pt_cols = 1;

    #if defined(__AVX512F__) && !defined(NO_CRT)
        assert(dim0 * ct_rows >= max_summed_pa_or_b_in_u64);

        for (size_t z = 0; z < poly_len; z++) {
            size_t idx_a_base = z * (ct_cols * dim0 * ct_rows);
            size_t idx_b_base = z * (num_per * pt_cols * dim0 * pt_rows);
            if (random_data) idx_b_base = (z % dummyWorkingSet) * (num_per * pt_cols * dim0 * pt_rows);

            for (size_t i = 0; i < num_per; i++) {
                for (size_t c = 0; c < pt_cols; c++) {
                    size_t uroll_max = max_summed_pa_or_b_in_u64;
                    size_t inner_limit = uroll_max;
                    size_t outer_limit = dim0 * ct_rows / inner_limit;

                    uint64_t sums_out_n0_u64_acc[4] = { 0, 0, 0, 0 };
                    uint64_t sums_out_n2_u64_acc[4] = { 0, 0, 0, 0 };

                    for (size_t o_jm = 0; o_jm < outer_limit; o_jm++) {
                        __m512i sums_out_n0 = _mm512_setzero_si512();
                        __m512i sums_out_n2 = _mm512_setzero_si512();

                        // #pragma unroll 16
                        for (size_t i_jm = 0; i_jm < inner_limit/8; i_jm++) {
                            size_t jm = o_jm * inner_limit + (8*i_jm);

                            uint64_t b_inp_1 = db[idx_b_base++];
                            uint64_t b_inp_2 = db[idx_b_base++];
                            uint64_t b_inp_3 = db[idx_b_base++];
                            uint64_t b_inp_4 = db[idx_b_base++];
                            __m512i b  = _mm512_set_epi64(b_inp_4, b_inp_4, b_inp_3, b_inp_3, b_inp_2, b_inp_2, b_inp_1, b_inp_1);

                            const uint64_t *v_a = &v_firstdim[idx_a_base + jm];

                            __m512i a = _mm512_load_si512(v_a);
                            __m512i a_lo = a;
                            __m512i a_hi_hi = _mm512_srli_epi64(a, packed_offset_2);
                            __m512i b_lo = b;
                            __m512i b_hi_hi = _mm512_srli_epi64(b, packed_offset_2);
                            
                            sums_out_n0 = _mm512_add_epi64(sums_out_n0, _mm512_mul_epu32(a_lo, b_lo));
                            sums_out_n2 = _mm512_add_epi64(sums_out_n2, _mm512_mul_epu32(a_hi_hi, b_hi_hi));
                        }
                        
                        // reduce here, otherwise we will overflow
                        alignas(64) uint64_t sums_out_n0_u64[8];
                        alignas(64) uint64_t sums_out_n2_u64[8];
                        _mm512_store_si512(sums_out_n0_u64, sums_out_n0);
                        _mm512_store_si512(sums_out_n2_u64, sums_out_n2);

                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n0_u64[idx] + sums_out_n0_u64[4 + idx];
                            sums_out_n0_u64_acc[idx] = (sums_out_n0_u64_acc[idx] + val) % p_i;
                        }
                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n2_u64[idx] + sums_out_n2_u64[4 + idx];
                            sums_out_n2_u64_acc[idx] = (sums_out_n2_u64_acc[idx] + val) % b_i;
                        }
                    }

                    for (size_t idx = 0; idx < 4; idx++) {
                        sums_out_n0_u64_acc[idx] %= p_i;
                        sums_out_n2_u64_acc[idx] %= b_i;
                    }

                    // output n0
                    size_t n = 0;
                    size_t idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i].data[idx_c] = barrett_coeff(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2], 0);
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i].data[idx_c] = barrett_coeff(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3], 0);

                    // output n1
                    n = 1;
                    idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i].data[idx_c] = barrett_coeff(sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2], 1);
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i].data[idx_c] = barrett_coeff(sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3], 1);
                }
            }
        }
    #elif defined(__AVX2__) && !defined(NO_CRT)
        assert(dim0 * ct_rows >= max_summed_pa_or_b_in_u64);

        for (size_t z = 0; z < poly_len; z++) {
            size_t idx_a_base = z * (ct_cols * dim0 * ct_rows);
            size_t idx_b_base = z * (num_per * pt_cols * dim0 * pt_rows);
            if (random_data) idx_b_base = (z % dummyWorkingSet) * (num_per * pt_cols * dim0 * pt_rows);

            for (size_t i = 0; i < num_per; i++) {
                for (size_t c = 0; c < pt_cols; c++) {
                    size_t uroll_max = max_summed_pa_or_b_in_u64;
                    size_t inner_limit = uroll_max;
                    size_t outer_limit = dim0 * ct_rows / inner_limit;

                    uint64_t sums_out_n0_u64_acc[4] = { 0, 0, 0, 0 };
                    uint64_t sums_out_n2_u64_acc[4] = { 0, 0, 0, 0 };

                    for (size_t o_jm = 0; o_jm < outer_limit; o_jm++) {
                        __m256i sums_out_n0 = _mm256_setzero_si256();
                        __m256i sums_out_n2 = _mm256_setzero_si256();

                        #pragma unroll 16
                        for (size_t i_jm = 0; i_jm < inner_limit/4; i_jm++) {
                            size_t jm = o_jm * inner_limit + (4*i_jm);

                            uint64_t b_inp_1 = db[idx_b_base++];
                            uint64_t b_inp_2 = db[idx_b_base++];
                            __m256i b  = _mm256_set_epi64x(b_inp_2, b_inp_2, b_inp_1, b_inp_1);

                            const uint64_t *v_a = &v_firstdim[idx_a_base + jm];

                            __m256i a = _mm256_load_si256((const __m256i *)v_a);
                            __m256i a_lo = a;
                            __m256i a_hi_hi = _mm256_srli_epi64(a, packed_offset_2);
                            __m256i b_lo = b;
                            __m256i b_hi_hi = _mm256_srli_epi64(b, packed_offset_2);
                            
                            sums_out_n0 = _mm256_add_epi64(sums_out_n0, _mm256_mul_epu32(a_lo, b_lo));
                            sums_out_n2 = _mm256_add_epi64(sums_out_n2, _mm256_mul_epu32(a_hi_hi, b_hi_hi));
                        }
                        
                        // reduce here, otherwise we will overflow
                        alignas(64) uint64_t sums_out_n0_u64[4];
                        alignas(64) uint64_t sums_out_n2_u64[4];
                        _mm256_store_si256((__m256i *)sums_out_n0_u64, sums_out_n0);
                        _mm256_store_si256((__m256i *)sums_out_n2_u64, sums_out_n2);

                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n0_u64[idx];
                            sums_out_n0_u64_acc[idx] = (sums_out_n0_u64_acc[idx] + val) % p_i;
                        }
                        for (size_t idx = 0; idx < 4; idx++) {
                            uint64_t val = sums_out_n2_u64[idx];
                            sums_out_n2_u64_acc[idx] = (sums_out_n2_u64_acc[idx] + val) % b_i;
                        }
                    }

                    for (size_t idx = 0; idx < 4; idx++) {
                        sums_out_n0_u64_acc[idx] %= p_i;
                        sums_out_n2_u64_acc[idx] %= b_i;
                    }

                    // output n0
                    size_t n = 0;
                    size_t idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i].data[idx_c] = barrett_coeff(sums_out_n0_u64_acc[0]+sums_out_n0_u64_acc[2], 0);
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i].data[idx_c] = barrett_coeff(sums_out_n0_u64_acc[1]+sums_out_n0_u64_acc[3], 0);

                    // output n1
                    n = 1;
                    idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i].data[idx_c] = barrett_coeff(sums_out_n2_u64_acc[0]+sums_out_n2_u64_acc[2], 1);
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i].data[idx_c] = barrett_coeff(sums_out_n2_u64_acc[1]+sums_out_n2_u64_acc[3], 1);
                }
            }
        }
    #else
        for (size_t z = 0; z < poly_len; z++) {
            size_t idx_a_base = z * (ct_cols * dim0 * ct_rows);
            size_t idx_b_base = z * (num_per * pt_cols * dim0 * pt_rows);
            if (random_data) idx_b_base = (z % dummyWorkingSet) * (num_per * pt_cols * dim0 * pt_rows);

            for (size_t i = 0; i < num_per; i++) {
                for (size_t c = 0; c < pt_cols; c++) {
                    // size_t idx_b = idx_b_base + i * (pt_cols * dim0 * n0) + c * (dim0 * n0);
                    
                    __uint128_t sums_out_n0_0 = 0, sums_out_n0_1 = 0;
                    __uint128_t sums_out_n1_0 = 0, sums_out_n1_1 = 0;

                    // #pragma unroll 16
                    for (size_t jm = 0; jm < dim0 * pt_rows; jm++) {
                        uint64_t b = db[idx_b_base++];
                        
                        const uint64_t *v_a = &v_firstdim[idx_a_base + jm * ct_rows];
                        uint64_t v_a0 = v_a[0];
                        uint64_t v_a1 = v_a[1];

                        uint32_t b_lo = b;
                        uint32_t b_hi = b >> 32L;

                        uint32_t v_a0_lo = v_a0;
                        uint32_t v_a0_hi = v_a0 >> 32L;

                        uint32_t v_a1_lo = v_a1;
                        uint32_t v_a1_hi = v_a1 >> 32L;
                        
                        // do n0
                        sums_out_n0_0 += (v_a0_lo) * (uint64_t)b_lo;
                        sums_out_n0_1 += (v_a1_lo) * (uint64_t)b_lo;

                        // do n1
                        sums_out_n1_0 += ((uint64_t)v_a0_hi) * b_hi;
                        sums_out_n1_1 += ((uint64_t)v_a1_hi) * b_hi;
                    }

                    // output n0
                    size_t n = 0;
                    size_t idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i].data[idx_c] = sums_out_n0_0 % p_i;
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i].data[idx_c] = sums_out_n0_1 % p_i;

                    // output n1
                    n = 1;
                    idx_c = c * (crt_count*poly_len) + n * (poly_len) + z;
                    out[i].data[idx_c] = sums_out_n1_0 % b_i;
                    idx_c += (pt_cols*crt_count*poly_len);
                    out[i].data[idx_c] = sums_out_n1_1 % b_i;
                }
            }
        }
    #endif
}


void foldCiphertextsDim1(
    vector<MatPoly> &v_cts,
    const vector<MatPoly> &v_folding,
    const vector<MatPoly> &v_folding_neg
) {
    size_t further_dims = ceil(log2(v_cts.size()));
    size_t ell = v_folding[0].cols / base_dim;
    MatPoly ginv_c(base_dim * ell, 1, false);
    MatPoly ginv_c_ntt(base_dim * ell, 1);
    MatPoly prod(base_dim, 1);
    MatPoly sum(base_dim, 1);

    size_t num_per = v_cts.size();
    for (size_t cur_dim = 0; cur_dim < further_dims; cur_dim++) {
        num_per = num_per / 2;
        for (size_t i = 0; i < num_per; i++) {
            const MatPoly &ct1 = v_cts[i];
            gadget_invert(base_dim * ell, ginv_c, ct1, base_dim);
            to_ntt(ginv_c_ntt, ginv_c);
            multiply(prod, v_folding_neg[further_dims - 1 - cur_dim], ginv_c_ntt);
            const MatPoly &ct2 = v_cts[num_per + i];
            gadget_invert(base_dim * ell, ginv_c, ct2, base_dim);
            to_ntt(ginv_c_ntt, ginv_c);
            multiply(sum, v_folding[further_dims - 1 - cur_dim], ginv_c_ntt);
            add(sum, sum, prod);
            v_cts[i] = from_ntt(sum);
        }
    }
}

void print_summary_testing(
    size_t num_items,
    size_t item_size_b,
    size_t total_offline_size_b,
    size_t total_query_size_b,
    size_t total_resp_size_b,
    size_t time_expansion_main,
    size_t time_conversion,
    size_t time_first_multiply,
    size_t time_folding,
    size_t time_packing,
    size_t time_key_gen,
    size_t time_query_gen,
    size_t time_decoding
) {
    size_t total_time_overall = time_first_multiply + time_folding + time_packing + time_expansion_main + time_conversion;
    double throughput_MB_per_s = ((double)num_items * item_size_b) / ((double)total_time_overall);
    double rate = ((double)item_size_b) / ((double)total_resp_size_b);

    cout << "ScalToMat took (CPU·us): 0" << endl;
    cout << "RegevToGSW took (CPU·us): 0" << endl;
    cout << "Expansion took (CPU·us): 0" << endl;
    cout << fixed << setprecision(0);
    cout << "Database" << endl;
    cout << endl;
    cout << "                       Number of items: " << num_items << endl;
    cout << "                             Item size: " << item_size_b << endl;
    cout << "Communication" << endl;
    cout << endl;
    cout << "         Total offline query size (b): " << total_offline_size_b << endl;
    cout << "          Total online query size (b): " << total_query_size_b << endl;
    cout << "                    Response size (b): " << total_resp_size_b  << endl;
    cout << fixed << setprecision(4);
    cout << "                                Rate : " << rate  << endl;
    cout << fixed << setprecision(0);
    cout << endl;
    cout << endl;
    cout << "Database-independent computation" << endl;
    cout << endl;
    cout << "              Main expansion  (CPU·us): " << time_expansion_main << endl;
    cout << "                   Conversion (CPU·us): " << time_conversion << endl;
    cout << "                        Total (CPU·us): " << (time_expansion_main + time_conversion) << endl;
    cout << endl;
    cout << "Database-dependent computation" << endl;
    cout << endl;
    cout << "     First dimension multiply (CPU·us): " << time_first_multiply << endl;
    cout << "                      Folding (CPU·us): " << time_folding << endl;
    cout << "                      Packing (CPU·us): " << time_packing << endl;
    cout << "                        Total (CPU·us): " << (time_first_multiply + time_folding + time_packing) << endl;
    cout << "                   Throughput (MB / s): " << throughput_MB_per_s << endl;
    cout << endl;
    cout << "Client computation" << endl;
    cout << endl;
    cout << "               Key generation (CPU·us): " << time_key_gen << endl;
    cout << "             Query generation (CPU·us): " << time_query_gen << endl;
    cout << "                     Decoding (CPU·us): " << time_decoding << endl;
    cout << endl;
}

inline size_t bits_to_bytes(size_t bits) {
    return (size_t)round((double)bits / 8.0);
}

static size_t time_expansion_main = 0;
static size_t time_conversion = 0;
static size_t time_first_multiply = 0;
static size_t time_folding = 0;
static size_t time_packing = 0;
static size_t time_key_gen = 0;
static size_t time_query_gen = 0;
static size_t time_decoding = 0;

static size_t time_db_gen = 0;

static size_t total_offline_size_b = 0;

void print_summary_testing(
    size_t nu_1,
    size_t nu_2,
    size_t logp,
    size_t logqprime,
    size_t out_n
) {
    size_t total_query_size_b = bits_to_bytes(((1<<nu_1) + 2*nu_2*t_GSW) * poly_len * logQ);
    if (query_elems_first == 1) total_query_size_b = bits_to_bytes(poly_len * logQ);
    size_t num_bits_resp = out_n * out_n * poly_len * (logp + 2);
    num_bits_resp += out_n * poly_len * logqprime;
    size_t total_resp_size_b = bits_to_bytes(num_bits_resp);

    size_t num_items = (1<<nu_1)*(1<<nu_2);
    size_t item_size_b = bits_to_bytes(out_n * out_n * poly_len * logp);

    print_summary_testing(
        num_items,
        item_size_b,
        total_offline_size_b,
        total_query_size_b,
        total_resp_size_b,
        time_expansion_main,
        time_conversion,
        time_first_multiply,
        time_folding,
        time_packing,
        time_key_gen,
        time_query_gen,
        time_decoding
    );
}

MatPoly recorrect_item(const MatPoly& pt) {
    MatPoly pt_corr(pt.rows, pt.cols, false);
    for (size_t pol = 0; pol < pt.rows * pt.cols * poly_len; pol++) {
        int64_t val = (int64_t) pt.data[pol];
        if ((val < 0) || (val >= Q_i)) {
            cout << pol << " " << val << endl;
        }
        assert((val >= 0) && (val < Q_i));
        if (val >= (Q_i / 2)) {
            val = val - Q_i;
        }
        if (val < 0) {
            val += p_db;
        }
        assert(val >= 0 && val < p_db);
        pt_corr.data[pol] = val;
    }
    return pt_corr;
}

void set_neg1s() {
    size_t num_expansions_max = 16;
    for (size_t i = 0; i < num_expansions_max; i++) {
        size_t idx = poly_len - (1 << i);
        MatPoly ng1(1, 1, false);
        ng1.data[idx] = 1; 
        MatPoly ng1_inv(1, 1, false);
        invert(ng1_inv, ng1);
        neg1s_mp_here.push_back(to_ntt(ng1_inv));
    }
}

static void add_pub_param(const MatPoly &mat) {
    total_offline_size_b += mat.rows * mat.cols * poly_len * logQ / 8;
}

static void add_pub_param(const vector<MatPoly> &vec) {
    for (size_t i = 0; i < vec.size(); i++) {
        add_pub_param(vec[i]);
    }
}

void testHighRate(size_t num_expansions, size_t further_dims, size_t idx_target) {
    cout << "Using n=" << out_n << endl;
    silent = true;

    size_t trials = out_n * out_n;

    size_t num_expanded = 1 << num_expansions;
    size_t total_n      = num_expanded * (1 << further_dims);
    size_t idx_dim0     = idx_target / (1 << further_dims);
    size_t idx_further  = idx_target % (1 << further_dims);
    
    size_t ell = t_GSW;
    size_t bits_per = get_bits_per(ell);

    size_t logp = (size_t)(ceil(log2(p_db)));

    size_t dim0 = 1 << num_expansions;
    size_t num_per = total_n / dim0;

    bool do_expansion = query_elems_first == 1;
    size_t num_bits_to_gen = ell * further_dims + num_expanded;
    size_t g = (size_t) ceil(log2((double)( num_bits_to_gen )));
    size_t stopround = (size_t) ceil(log2((double)( ell * further_dims )));
    size_t pts_to_gen = random_data ? min(1<<16, total_n): total_n;
    dummyWorkingSet = poly_len * pts_to_gen / total_n;

    set_neg1s();

    // === allocations ===
    // keygen / querygen
    MatPoly G_conv = buildGadget(base_dim, base_dim * m_conv);
    MatPoly g_vec = buildGadget(1, m_conv);
    MatPoly g_vec_ntt = to_ntt(g_vec);
    vector<MatPoly> v_W;

    // possible expansion
    vector<MatPoly> v_firstdim_and_folding;
    vector<MatPoly> v_W_exp_left;
    vector<MatPoly> v_W_exp_right;
    MatPoly V(base_dim, base_dim * m_conv);
    MatPoly single_query_ct(base_dim, 1);

    // first dim
    vector<MatPoly> v_firstdim;
    uint64_t *v_firstdim_raw = alloc_u64(num_expanded * base_dim * 1 * poly_len);
    
    // gsw encoding
    MatPoly gadget(base_dim, base_dim * ell, false);
    buildGadget(gadget);
    MatPoly gadget_ntt = to_ntt(gadget);

    // folding
    vector<MatPoly> v_folding;
    vector<MatPoly> v_folding_neg;
    MatPoly ct_gsw_inv(base_dim, base_dim * ell, false);
    vector<MatPoly> v_result_ct;

    // packing + switching
    MatPoly packed_ct(out_n + 1, out_n);
    MatPoly total_resp(out_n + 1, out_n, false);

    // === db generation ===
    st();
    vector<vector<MatPoly> > v_db;
    vector<uint64_t *> v_db_raw;
    for (size_t trial = 0; trial < trials; trial++) {
        MatPoly pt(1, 1, false);
        MatPoly pt_encd(1, 1);
        v_db.push_back(vector<MatPoly>());
        for (size_t i = 0; i < pts_to_gen; i++) {
            if (!random_data || (i == 0)) {
                if (random_data) {
                    random_device rd;
                    uniform_int_distribution<uint64_t> dist(numeric_limits<uint64_t>::min(), numeric_limits<uint64_t>::max());
                    uint64_t val = dist(rd) % (p_db/4);
                    pt.data[0] = val;
                } else {
                    uniform_matrix(pt);
                    reduce_mod(pt, p_db);
                }
                for (size_t pol = 0; pol < 1 * 1 * poly_len; pol++) {
                    int64_t val = (int64_t) pt.data[pol];
                    assert(val >= 0 && val < p_db);
                    if (val >= (p_db / 2)) {
                        val = val - (int64_t)p_db;
                    }
                    if (val < 0) {
                        val += Q_i;
                    }
                    assert(val >= 0 && val < Q_i);
                    pt.data[pol] = val;
                }
                to_ntt(pt_encd, pt);
            }
            v_db[trial].push_back(pt_encd);
        }
        if (random_data) {
            size_t pt_rows = v_db[trial][0].rows;
            size_t pt_cols = v_db[trial][0].cols;
            uint64_t *db_buf = (uint64_t *)malloc(total_n * pt_rows * pt_cols * dummyWorkingSet * sizeof(uint64_t));
            // b': i c n z j m
            for (size_t i = 0; i < total_n; i++) {
                size_t ii = i % num_per;
                size_t j = i / num_per;
                for (size_t m = 0; m < pt_rows; m++) {
                    for (size_t c = 0; c < pt_cols; c++) {
                        uint64_t *BB = &(v_db[trial][0].data)[(m * pt_cols + c) * crt_count * coeff_count];
                        for (size_t z = 0; z < dummyWorkingSet; z++) {
                            size_t idx = z  * (num_per * pt_cols * dim0 * pt_rows) +
                                        ii * (pt_cols * dim0 * pt_rows) + 
                                        c  * (dim0 * pt_rows) +
                                        j  * (pt_rows) + 
                                        m;
                            
                            db_buf[idx] = BB[z] | (BB[poly_len + z] << 32);
                        }
                    }
                }
            }
            v_db_raw.push_back(db_buf);
            if (trial >= max_trials) break;
        } else {
            v_db_raw.push_back(convertDb(v_db[trial], dim0, num_per));
        }
    }
    time_db_gen = et();

    // === key generation ===
    st();
    S_mp = MatPoly(out_n, out_n+1, false);
    Sp_mp = MatPoly(out_n, 1, false);
    sr_mp = MatPoly(1, 1, false);
    keygen(S_mp, Sp_mp, sr_mp, out_n);
    MatPoly s0 = sr_mp;
    MatPoly s0_ntt = to_ntt(s0);
    MatPoly s0_1 = getSkVec(0);
    MatPoly Sp_mp_nttd_qprime(out_n, 1, false);
    Sp_mp_nttd_qprime = Sp_mp;
    to_ntt_qprime(Sp_mp_nttd_qprime);

    // generate Ws
    for (size_t i = 0; i < out_n; i++) {
        MatPoly AG(out_n, g_vec_ntt.cols);
        MatPoly prod = mul_by_const(s0_ntt, g_vec_ntt);
        place(AG, prod, i, 0);
        MatPoly W = encryptMatrixArbitrary(AG, out_n);
        v_W.push_back(W);
    }
    add_pub_param(v_W);
    if (do_expansion) {
        MatPoly s0_sq_ntt = multiply(s0_ntt, s0_ntt);
        getExpansionKeySwitchingMatrices(v_W_exp_left, g, m_exp);
        getExpansionKeySwitchingMatrices(v_W_exp_right, stopround + 1, m_exp_right);

        // Generate V
        MatPoly sigma(1, 1, false);
        for (size_t i = 0; i < base_dim * m_conv; i++) {
            if (i % 2 == 0) {
                uint64_t val = G_conv.data[i * poly_len];
                from_ntt(sigma, multiply(s0_sq_ntt, single_poly(val)));
            } else {
                uint64_t val = G_conv.data[(base_dim * m_conv + i) * poly_len];
                from_ntt(sigma, multiply(s0_ntt, single_poly(val)));
            }
            MatPoly ct = encryptSimpleRegev(sigma);
            place(V, ct, 0, i);
        }

        add_pub_param(v_W_exp_left);
        add_pub_param(v_W_exp_right);
        add_pub_param(V);
    }
    time_key_gen += et();
    
    // === query generation ===
    cout << "query: (";
    cout << idx_dim0 << " ";
    for (size_t i = 0; i < further_dims; i++) {
        uint64_t bit = (idx_further & (1 << i)) >> i;
        cout << bit << " ";
    }
    cout << ")" << endl;

    st();
    if (!do_expansion) {
        // generate cts for nu_1
        for (size_t i = 0; i < num_expanded; i++) {
            MatPoly ct = encryptSimpleRegev(single_poly(i == idx_dim0 ? scale_k : 0));
            v_firstdim.push_back(ct);
        }

        // generate cts for rest of dimensions
        for (size_t i = 0; i < further_dims; i++) {
            uint64_t bit = (idx_further & (1 << i)) >> i;
            MatPoly ct_gsw(base_dim, base_dim * ell);
            for (size_t j = 0; j < ell; j++) {
                uint64_t val = (1UL << (bits_per * j)) * bit;
                MatPoly ct = encryptSimpleRegev(single_poly(val));
                place(ct_gsw, ct, 0, base_dim*j + 1);
                MatPoly prod = multiply(s0, single_poly(val));
                ct = encryptSimpleRegev(from_ntt(prod));
                place(ct_gsw, ct, 0, base_dim*j);
            }
            v_folding.push_back(ct_gsw);
        }

        // reorient (could be sped up)
        reorientCiphertextsDim1(v_firstdim_raw, v_firstdim, dim0, 1);
    } else {
        // pack query into single ciphertext
        MatPoly sigma(1, 1, false);
        sigma.data[2*idx_dim0] = scale_k;
        for (size_t i = 0; i < further_dims; i++) {
            uint64_t bit = (idx_further & (1 << i)) >> i;
            for (size_t j = 0; j < ell; j++) {
                uint64_t val = (1UL << (bits_per * j)) * bit;
                sigma.data[2*(i*ell + j) + 1] = val;
            }
        }

        uint64_t inv_2_g_first = inv_mod(1 << g, Q_i);
        uint64_t inv_2_g_rest = inv_mod(1 << (stopround+1), Q_i);
        for (size_t i = 0; i < coeff_count/2; i++) {
            sigma.data[2*i]   = (sigma.data[2*i] * (__uint128_t)inv_2_g_first) % Q_i;
            sigma.data[2*i+1] = (sigma.data[2*i+1] * (__uint128_t)inv_2_g_rest) % Q_i;
        }
        single_query_ct = encryptSimpleRegev(sigma);
    }
    time_query_gen += et();

    // === conversion ===
    st();
    if (do_expansion) {
        // unpack query
        v_firstdim_and_folding.push_back(single_query_ct);
        for (size_t i = 0; i < (1<<g)-1; i++) {
            v_firstdim_and_folding.emplace_back(base_dim, 1);   
        }
        coefficientExpansion(v_firstdim_and_folding, g, m_exp, v_W_exp_left, v_W_exp_right, ell * further_dims, stopround);

        // rearrange output regev ciphertexts
        reorientCiphertextsDim1(v_firstdim_raw, v_firstdim_and_folding, dim0, 2);
    }
    time_expansion_main += et();
    st();
    if (do_expansion) {
        // convert GSW ciphertexts
        regevToSimpleGsw(v_folding, v_firstdim_and_folding, V, m_conv, ell, further_dims, base_dim, 1);
    }
    // invert (technically a query expansion step)
    for (size_t i = 0; i < further_dims; i++) {
        invert(ct_gsw_inv, from_ntt(v_folding[i]));
        MatPoly ct_gsw_neg(base_dim, base_dim * ell);
        add(ct_gsw_neg, gadget_ntt, to_ntt(ct_gsw_inv));
        v_folding_neg.push_back(ct_gsw_neg);
    }
    time_conversion += et();

    vector<MatPoly> v_out;
    vector<MatPoly> v_out_raw;
    for (size_t i = 0; i < (1 << further_dims); i++) {
        MatPoly zero(base_dim, 1);
        MatPoly zero_raw(base_dim, 1, false);
        v_out.push_back(zero);
        v_out_raw.push_back(zero_raw);
    }
    MatPoly result_ct(base_dim, 1, false);
    // === processing ===
    for (size_t trial = 0; trial < trials; trial++) {
        size_t cur_trial = random_data ? (trial % max_trials) : trial;

        st();
        // multiplyQueryByDatabaseDim1(v_out, v_db[trial], v_firstdim);
        fastMultiplyQueryByDatabaseDim1(v_out, v_db_raw[cur_trial], v_firstdim_raw, dim0, num_per);
        time_first_multiply += et();

        st();
        // ntt invert all input
        for (size_t i = 0; i < (1 << further_dims); i++) {
            from_ntt(v_out_raw[i], v_out[i]);
        }
        foldCiphertextsDim1(v_out_raw, v_folding, v_folding_neg);
        result_ct = v_out_raw[0];
        v_result_ct.push_back(result_ct);
        time_folding += et();
    }

    // === pack ===
    st();
    pack(packed_ct, out_n, m_conv, v_result_ct, v_W);

    // modulus switch all but first row to 'p'
    size_t qprime_bits = bits_to_hold_arb_qprime;
    uint64_t arb_qprime_h = qprime_mods[qprime_bits];

    uint64_t q_1 = 4*p_db;

    MatPoly ct_inp = from_ntt(packed_ct);
    MatPoly first_row = pick(ct_inp, 0, 0, 1, ct_inp.cols);
    MatPoly first_row_sw = getRescaled(first_row, Q_i, arb_qprime_h);
    MatPoly rest_rows = pick(ct_inp, 1, 0, ct_inp.rows - 1, ct_inp.cols);
    MatPoly rest_rows_sw = getRescaled(rest_rows, Q_i, q_1);

    place(total_resp, first_row_sw, 0, 0);
    place(total_resp, rest_rows_sw, 1, 0);
    time_packing += et();

    // ~~transmit~~

    st();
    MatPoly first_row_decoded = pick(total_resp, 0, 0, 1, total_resp.cols);
    MatPoly rest_rows_decoded = pick(total_resp, 1, 0, total_resp.rows - 1, total_resp.cols);

    // MatPoly first_row_scaled = getRescaled(first_row_decoded, arb_qprime_h, Q_i);
    // MatPoly rest_rows_scaled = getRescaled(rest_rows_decoded, q_1, p_db);
    // MatPoly s_prod = mul_over_integers(Sp_mp, first_row_decoded, arb_qprime_h);
    // MatPoly s_prod_sw = getRescaled(s_prod, arb_qprime_h, p_db);
    to_ntt_qprime(first_row_decoded);
    MatPoly s_prod = mul_over_qprime(Sp_mp_nttd_qprime, first_row_decoded);
    from_ntt_qprime(s_prod);
    for (size_t i = 0; i < s_prod.rows * s_prod.cols * poly_len; i++) {
        int64_t val_first = s_prod.data[i];
        if (val_first >= arb_qprime_h/2) val_first -= arb_qprime_h;
        int64_t val_rest = rest_rows_decoded.data[i];
        if (val_rest >= q_1/2) val_rest -= q_1;

        uint64_t denom = arb_qprime_h * (q_1/p_db);

        int64_t r = val_first * q_1;
        r +=  val_rest * arb_qprime_h;
        // divide r by arb_qprime_h, rounding
        int64_t sign = r >= 0 ? 1 : -1;
        __int128_t val = r;
        __int128_t result = (r + sign*((int64_t)denom/2)) / (__int128_t)denom;
        result = (result + (denom/p_db)*p_db + 2*p_db) % p_db;

        s_prod.data[i] = (uint64_t)result;
    }
    // MatPoly ct_out(out_n+1, out_n, false);
    // place(ct_out, first_row_scaled, 0, 0);
    // place(ct_out, rest_rows_scaled, 1, 0);

    // MatPoly Z = multiply(S_mp, ct_out);
    MatPoly Z_uncrtd = s_prod;//from_ntt(Z);
    // divide_by_const(Z_uncrtd, Z_uncrtd, p_db, Q_i, p_db);
    time_decoding += et();

    // MatPoly top_left = pick(Z_uncrtd, 0, 0, 1, 1);
    MatPoly corrected_items(out_n, out_n, false);
    for (size_t i = 0; i < out_n; i++) {
        for (size_t j = 0; j < out_n; j++) {
            size_t idx_fetch = random_data ? ((i * out_n + j) % max_trials) : (i * out_n + j);
            size_t idx_get_target = random_data ? 0 : idx_target;
            MatPoly corrected_item = recorrect_item(from_ntt(v_db[idx_fetch][idx_get_target]));
            place(corrected_items, corrected_item, i, j);
        }
    }
    double log_var = get_log_var(Z_uncrtd, corrected_items, show_diff, p_db);
    cout << "variance of diff with modswitching: 2^" << log_var << endl;
    cout << "Is correct? : " << is_eq(Z_uncrtd, corrected_items) << endl;
    size_t max_ctr = 0;
    for (size_t i = 0; i < Z_uncrtd.rows * Z_uncrtd.cols * poly_len; i++) {
        if (Z_uncrtd.data[i] != corrected_items.data[i]) {
            cout << i << ": " << Z_uncrtd.data[i] << " " << corrected_items.data[i] << endl;
            if (++max_ctr >= 10) {
                cout << "... more" << endl;
                break;
            }
        }
    }

    print_summary_testing(
        num_expansions,
        further_dims,
        logp,
        bits_to_hold_arb_qprime,
        out_n
    );
}