#include "client.h"

MatPoly S_mp;
MatPoly Sp_mp;
MatPoly sr_mp;

uint64_t sample_u64() {
    if (nonoise) return 0;
    return (sample() + Q_i) % Q_i;
}

void noise(MatPoly &E) {
    assert(!E.isNTT);

    size_t rs = E.rows;
    size_t cs = E.cols;

    for (size_t i = 0; i < E.rows * E.cols * poly_len; i++) {
        E.data[i] = sample_u64();
    }
}

void keygen(MatPoly &S, MatPoly &Sp, MatPoly &sr, size_t n_val) {
    size_t rs = S.rows; // n0
    size_t cs = S.cols; // n1

    for (size_t m = 0; m < poly_len; m++) {
        uint64_t val = ternary ? ((rand() % 3) - 1) : sample_u64();
        sr.data[m] = val % Q_i;
    }

    for (size_t r = 0; r < rs; r++) {
        for (size_t c = 0; c < cs - n_val; c++) {
            for (size_t m = 0; m < poly_len; m++) {
                uint64_t val = ternary ? ((rand() % 3) - 1) : sample_u64();
                S.data[r * cs * poly_len + c * poly_len + m] = val % Q_i;
                Sp.data[r * Sp.cols * poly_len + c * poly_len + m] = val % Q_i;
            }
        }
    }
    for (size_t r = 0; r < rs; r++) {
        for (size_t c = 0; c < n_val; c++) {
            S.data[r * cs * poly_len + (c + 1) * poly_len] = r == c ? 1 : 0;
        }
    }
}

MatPoly get_fresh_public_key_raw(MatPoly &Sp, long m) {
    MatPoly A(k_param, m, false);
    uniform_matrix(A);
    MatPoly E(n0, m, false);
    noise(E);

    MatPoly A_ntt = to_ntt(A);
    MatPoly B_p(n0, m);
    multiply(B_p, to_ntt(Sp), A_ntt);
    MatPoly B(n0, m);
    add(B, to_ntt(E), B_p);
    
    MatPoly A_inv_mp_raw(A_ntt.rows, A_ntt.cols, false);
    invert(A_inv_mp_raw, from_ntt(A_ntt));

    MatPoly P(n1, m, false);
    vertical_merge(P, A_inv_mp_raw, from_ntt(B));

    return P;
}

MatPoly enc_scalar(const MatPoly &P, const MatPoly &G, const MatPoly &sigma, size_t mm, size_t num_expansions) {
    assert(!P.isNTT);
    assert(G.isNTT);

    MatPoly prod(G.rows, mm);
    mul_by_const(prod, to_ntt(sigma), G);

    MatPoly padding(k_param, mm);
    MatPoly padded_prod(G.rows + k_param, mm);
    vertical_merge(padded_prod, padding, prod);

    // multiply by 2^num_expansions, so that div by 2 causes no rounding later
    MatPoly P_mul_by_const(P.rows, P.cols);
    mul_by_const(P_mul_by_const, to_ntt(single_poly((long)(((long)1) << num_expansions))), to_ntt(P));

    MatPoly C(G.rows + k_param, mm);
    add(C, P_mul_by_const, padded_prod);

    return C;
}

void load_modswitched_into_ct(MatPoly &ct, const uint64_t *modswitched) {
    assert(!ct.isNTT);

    size_t rs = n1;
    size_t cs = n2;

    size_t bit_offs = 0;
    size_t bit_width = bits_to_hold_arb_qprime;
    // cout << "POSTTT:";
    for (size_t r = 0; r < rs; r++) {
        for (size_t c = 0; c < cs; c++) {
            for (size_t m = 0; m < poly_len; m++) {
                // uint64_t result = modswitched[r * cs * poly_len + c * poly_len + m] = result;
                uint64_t result = read_arbitrary_bits(modswitched, bit_offs, bit_width);
                // if (result > (q_i * switch_factor / 2)) result = Q_i_u128 + result - (q_i * switch_factor);
                // cout << result << " ";
                bit_offs += bit_width;
                ct.data[r * cs * poly_len + c * poly_len + m] = result;
            }
        }
    }
    // cout << endl;
}

// C should be 3x3
void dec_compressed(MatPoly &O_out, const MatPoly &C, const MatPoly &S, uint64_t divisor) {
    assert(C.isNTT);
    assert(S.isNTT);
    
    // apply secret key
    MatPoly O(S.rows, C.cols);
    multiply(O, S, C);
    MatPoly O_raw = from_ntt(O);

    // do division here if we are trying to avoid round-off error

    if (divisor > 1) {
        // reduce_mod(O_raw, q_const * divisor);
        divide_by_const(O_raw, O_raw, divisor, p_db);
    }

    reduce_mod(O_raw, p_db);

    assert(n2 == n0);

    O_out = O_raw;
}

void dec(MatPoly &Z, const MatPoly &S, const MatPoly &C) {
    assert(Z.isNTT);
    assert(S.isNTT);
    assert(C.isNTT);

    // apply secret key
    multiply(Z, S, C);//C_s);
}

MatPoly getRegevSample(const MatPoly &s) { // outputs n0 x 1
    assert(s.rows == 1);
    assert(s.cols == 1);
    MatPoly a(1, 1, false);
    uniform_matrix(a);
    MatPoly e(1, 1, false);
    noise(e);
    MatPoly a_inv(1, 1, false);
    invert(a_inv, a);
    MatPoly b = multiply(a, s);
    add(b, b, to_ntt(e));

    MatPoly P(n0, 1);
    place(P, to_ntt(a_inv), 0, 0);
    place(P, b, 1, 0);
    return P;
}

MatPoly getRegevPublicKeyMatrix(const MatPoly &s, size_t m) {
    MatPoly C(n0, m);

    for (size_t i = 0; i < m; i++) {
        MatPoly P = getRegevSample(s);
        place(C, P, 0, i);
    }

    return C;
}

MatPoly encryptSimpleRegev(MatPoly sigma, size_t noise_factor) {
    // uses s_1 as the secret key (just a single element)
    MatPoly s = sr_mp;
    
    MatPoly P = getRegevSample(s);
    if (noise_factor > 1) {
        P = mul_by_const(to_ntt(single_poly(noise_factor)), P);
    }

    MatPoly sigma_padded(n0, 1, false);
    place(sigma_padded, sigma, 1, 0);

    MatPoly result(n0, 1);
    add(result, P, to_ntt(sigma_padded));

    return result;
}

MatPoly encryptSimpleRegevMatrix(const MatPoly &s, const MatPoly &G_conv, MatPoly sigma_m, size_t noise_factor) {
    assert(G_conv.rows == n0);
    assert(s.rows == 1);
    assert(s.cols == 1);

    size_t m_conv = G_conv.cols;
    MatPoly P = getRegevPublicKeyMatrix(s, m_conv);
    if (noise_factor > 1) {
        P = mul_by_const(to_ntt(single_poly(noise_factor)), P);
    }

    MatPoly prod = multiply(sigma_m, G_conv);
    MatPoly padded_prod(n0, m_conv);
    place(padded_prod, prod, 1, 0);

    MatPoly result(n0, m_conv);
    add(result, P, padded_prod);

    return result;
}

MatPoly encryptSimpleRegevMatrix(const MatPoly &s, const MatPoly &mat_to_enc, size_t noise_factor) {
    assert(mat_to_enc.isNTT);
    assert(mat_to_enc.rows == 1);
    assert(s.rows == 1);
    assert(s.cols == 1);

    size_t m = mat_to_enc.cols;
    MatPoly P = getRegevPublicKeyMatrix(s, m);
    if (noise_factor > 1) {
        P = mul_by_const(to_ntt(single_poly(noise_factor)), P);
    }

    MatPoly padded_prod(n0, m);
    place(padded_prod, mat_to_enc, 1, 0);

    MatPoly result(n0, m);
    add(result, P, padded_prod);

    return result;
}

void generate_X_improved(
    size_t mx,
    MatPoly &X,     // n1 x mx
    MatPoly A,      // n0 x n0
    MatPoly G       // n0 x mx
) {
    assert(X.isNTT);
    assert(A.isNTT);
    assert(G.isNTT);
    assert(G.cols == mx);
    assert(X.cols == mx);
    assert(A.rows == n0);
    assert(A.cols == n0);
    MatPoly P = get_fresh_public_key_raw(Sp_mp, mx);

    MatPoly AG = multiply(A, G);

    MatPoly padding(n1-n0, mx);
    MatPoly AG_padded(n1, mx);
    vertical_merge(AG_padded, padding, AG);

    MatPoly result(n1, mx);
    add(result, AG_padded, to_ntt(P));

    X = result;
}

MatPoly getSkVec(size_t idx) {
    MatPoly s0 = sr_mp;
    MatPoly s0_1(1, n0, false);
    place(s0_1, s0, 0, 0);
    place(s0_1, single_poly(1), 0, 1);
    return s0_1;
}

void getPublicEncryptions(
    size_t g,
    vector<MatPoly> &X_v,
    vector<MatPoly> &Y_v,
    vector<MatPoly> &W_v,
    vector<MatPoly> &W_exp_v,
    size_t m_conv,
    size_t m_exp,
    bool for_composition,
    bool for_exp,
    size_t stopround
) {
    MatPoly G_conv = buildGadget(n0, n0 * m_conv);
    MatPoly G_exp = buildGadget(1, m_exp);
    MatPoly G_exp_nttd = to_ntt(G_exp);

    for (size_t i = 0; i < (stopround == 0 ? g : stopround); i++) {
        size_t t = (poly_len / (1 << i)) + 1;

        MatPoly s0 = sr_mp;
        MatPoly tau_s0 = automorph(s0, t);
        size_t noise_factor = 1;// << (g - i);

        MatPoly W_exp_i = encryptSimpleRegevMatrix(s0, multiply(tau_s0, G_exp_nttd), noise_factor);
        W_exp_v.push_back(W_exp_i);
    }

    if (for_exp) return;

    // MatPoly A(n0, n0, false);
    // memcpy(A[0], Sp_mp[0], poly_len * sizeof(uint64_t));
    // A[0][1 * poly_len] = 1;
    // MatPoly B(n0, n0, false);
    // if (for_composition) {
    //     memcpy(B[1], Sp_mp[0], poly_len * sizeof(uint64_t));
    // } else {
    //     memcpy(B[1], Sp_mp[1], poly_len * sizeof(uint64_t));
    // }
    // B[1][1 * poly_len] = 1;

    // MatPoly Xa(n1, n0 * m_conv);
    // MatPoly Xb(n1, n0 * m_conv);
    // generate_X_improved(n0 * m_conv, Xa, to_ntt(A), to_ntt(G_conv));
    // generate_X_improved(n0 * m_conv, Xb, to_ntt(B), to_ntt(G_conv));
    // X_v.resize(n0);
    // X_v[0] = Xa;
    // X_v[1] = Xb;

    // MatPoly s0_1 = getSkVec(0);

    // for (size_t i = 0; i < n0; i++) {
    //     MatPoly s_i = pick(Sp_mp, i, 0, 1, 1);
    //     MatPoly Y_pt = multiply(s_i, s0_1);

    //     MatPoly Y_i = encryptSimpleRegevMatrix(s_i, G_conv, Y_pt);
    //     Y_v.push_back(Y_i);
    //     total_offline_query_size += Y_i.rows * Y_i.cols * coeff_count * logQ / 8;
    // }

    // for (size_t i = 1; i < n0; i++) {
    //     MatPoly s_i = pick(Sp_mp, i, 0, 1, 1);

    //     MatPoly W_i = encryptSimpleRegevMatrix(s_i, G_conv, s0_1);
    //     W_v.push_back(W_i);
    //     total_offline_query_size += W_i.rows * W_i.cols * coeff_count * logQ / 8;
    // }
}