#include "util.h"

size_t dummyWorkingSet = 128;
size_t max_trials = 16;

int64_t variance_raw(const vector<int64_t> &x) {
    double sum = 0;
    for (size_t m = 0; m < x.size(); m++) {
        int64_t v = x[m];
        sum += (double) v;
    }
    double mean = sum / (double) x.size();
    double var = 0;
    for (size_t m = 0; m < x.size(); m++) {
        double v = (double) x[m];
        double v2 = v - mean;
        var += v2 * v2;
    }
    var = var / (double) x.size();
    return (int64_t)var;
}

double get_log_var(const MatPoly &A, const MatPoly &B, bool print, size_t mod) {
    assert(!A.isNTT);
    assert(!B.isNTT);
    assert(A.rows == B.rows);
    assert(A.cols == B.cols);

    __uint128_t sum = 0;
    __uint128_t sum_real = 0;
    size_t count = 0;

    vector<int64_t> v;
    if (print) cout << ">>>";
    for (size_t i = 0; i < A.rows * A.cols * poly_len; i++) {
        int64_t a = (int64_t) (A.data[i] % mod);
        int64_t b = (int64_t) (B.data[i] % mod);
        // if (print) cout << a << "," << b << " ";
        if (a >= mod / 2)
            a -= mod;
        if (b >= mod / 2)
            b -= mod;
        int64_t diff = abs((long)(a+b)) < abs((long)(a - b)) ? (a+b) : a - b;
        // if (print) cout << a << "," << b << " (" << abs(diff) << ") ";
        if (print) cout << (diff) << ", ";
        sum += abs(diff);
        sum_real += (a - b > mod / 2) ? a - b - mod : a - b;
        count += 1;

        v.push_back(diff);
    }
    if (print) cout << endl;
    double var = (double) variance_raw(v);
    double lg = var != 0 ? (log(var) / log(2)) : 0;
    return lg;
}

vector<int64_t> get_diffs(const MatPoly &A, const MatPoly &B, uint64_t modulus) {
    vector<int64_t> v;
    for (size_t i = 0; i < A.rows * A.cols * poly_len; i++) {
        int64_t a = (int64_t) (A.data[i] % modulus);
        int64_t b = (int64_t) (B.data[i] % modulus);
        if (a >= modulus / 2) a -= modulus;
        if (b >= modulus / 2) b -= modulus;
        int64_t diff = abs((long)(a+b)) < abs((long)(a-b)) ? (a+b) : a - b;
        v.push_back(diff);
    }
    return v;
}

void output_diffs(ostream &os, vector<int64_t> diffs) {
    for (size_t i = 0; i < diffs.size(); i++) {
        os << diffs[i] << " ";
    }
}

random_device rd;
mt19937_64 gen2(rd()); // NOT SECURE
uniform_int_distribution<uint64_t> dist(numeric_limits<uint64_t>::min(), numeric_limits<uint64_t>::max());

void uniform_matrix(MatPoly &A) {
    assert(!A.isNTT);

    for (size_t i = 0; i < A.rows * A.cols * poly_len; i++) {
        A.data[i] = dist(gen2) % Q_i;
    }
}

void buildGadget(MatPoly &G) {
    assert(!G.isNTT);
    size_t nx = G.rows;
    size_t m = G.cols;

    assert(m % nx == 0);

    size_t num_elems = m / nx;
    size_t bits_per = get_bits_per(num_elems);
    uint64_t bound = 1UL << bits_per;

    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < num_elems; j++) {
            if (bits_per * j >= 64) continue;
            G[i][(i + j * nx) * poly_len] = 1UL << (bits_per * j);
        }
    }
}

MatPoly buildGadget(size_t rows, size_t cols) {
    MatPoly G(rows, cols, false);
    buildGadget(G);
    return G;
}

void gadget_invert( size_t mx, MatPoly &out, const MatPoly &in, size_t rdim) {
    size_t m = out.cols;
    assert(in.rows == rdim);
    assert(in.cols == m);
    assert(out.rows == mx);
    assert(!out.isNTT);
    assert(!in.isNTT);
    assert(mx % rdim == 0);

    size_t num_elems = mx / rdim;
    size_t bits_per = get_bits_per(num_elems);
    uint64_t mask = (1UL << bits_per) - 1;

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < rdim; j++) {
            // decompose this value into num_elems values
            
            for (size_t z = 0; z < poly_len; z++) {
                uint64_t val = in[j][i * poly_len + z];

                for (size_t k = 0; k < num_elems; k++) {
                    size_t row = j + k * rdim;
                    size_t bit_offs = min(k * bits_per, 64);
                    uint64_t piece = (val >> bit_offs) & mask;
                    
                    out[row][i * poly_len + z] = piece;
                }
            }
        }
    }
}

MatPoly gadget_invert(size_t mx, size_t m, const MatPoly &in, size_t rdim) {
    MatPoly out(mx, m, false);
    gadget_invert(mx, out, in, rdim);
    return out;
}

void poly_mul_integers(int64_t *res, const uint64_t *a, const uint64_t *b, uint64_t modulus) {
    for (size_t i = 0; i < poly_len; i++) {
        for (size_t j = 0; j < poly_len; j++) {
            size_t idx = (i + j) % poly_len;
            
            int64_t a_val = a[i];
            // cout << a_val << " -> ";
            if (a_val >= Q_i/2) a_val -= (int64_t)Q_i;
            a_val = (a_val + (Q_i/modulus)*modulus + 2*modulus) % modulus;
            // cout << a_val << endl;
            int64_t b_val = b[j];
            // if (b_val > Q_i/2) b_val -= (int64_t)Q_i;
            // b_val = (b_val + (Q_i/modulus)*modulus + 2*modulus) % modulus;

            assert(a_val >= 0);
            assert(b_val >= 0);

            uint64_t prod = (a_val * (__uint128_t) b_val) % modulus;

            if ((i + j) < poly_len) {
                res[idx] = (res[idx] + prod + modulus) % modulus;
            } else {
                res[idx] = (res[idx] + (modulus - prod) + modulus) % modulus;
            }
        }
    }
}

MatPoly mul_over_integers(const MatPoly &a, const MatPoly &b, uint64_t modulus, bool b_is_single_poly) {
    assert(!a.isNTT);
    assert(!b.isNTT);
    if (!b_is_single_poly) assert(a.cols == b.rows);
    MatPoly res(a.rows, b.cols, false);

    if (b_is_single_poly) res = MatPoly(a.rows, a.cols, false);

    int64_t *scratch = (int64_t *)calloc(poly_len, sizeof(int64_t));

    for (size_t r = 0; r < a.rows; r++) {
        for (size_t c = 0; c < b.cols; c++) {
            for (size_t m = 0; m < a.cols; m++) {
                uint64_t *a_p   = &a.data[r * a.cols * poly_len + m * poly_len];
                uint64_t *b_p   = &b.data[m * b.cols * poly_len + c * poly_len];
                if (b_is_single_poly) b_p = b.data;
                uint64_t *res_p = &res.data[r * res.cols * poly_len + c * poly_len];
                for (size_t i = 0; i < poly_len; i++) {
                    scratch[i] = 0;
                }
                poly_mul_integers(scratch, a_p, b_p, modulus);
                for (size_t i = 0; i < poly_len; i++) {
                    res_p[i] = (res_p[i] + scratch[i]) % modulus;
                }
            }
        }
    }

    free(scratch);

    return res;
}

void poly_mul_qprime(int64_t *res, const uint64_t *a, const uint64_t *b) {
    for (size_t i = 0; i < poly_len; i++) {
        res[i] = (a[i] * (__uint128_t)b[i]) % arb_qprime;
    }
}

// WARNING: does not change NTT property
void to_ntt_qprime(MatPoly &a) {
    assert(!a.isNTT);
    for (size_t r = 0; r < a.rows; r++) {
        for (size_t c = 0; c < a.cols; c++) {
            uint64_t *a_p   = &a.data[r * a.cols * poly_len + c * poly_len];
            for (size_t z = 0; z < poly_len; z++) {
                int64_t a_val = a_p[z];
                if (a_val >= Q_i/2) a_val -= (int64_t)Q_i;
                a_val = (a_val + (Q_i/arb_qprime)*arb_qprime + 2*arb_qprime) % arb_qprime;
                a_p[z] = a_val;
            }
            ntt_qprime->ComputeForward(a_p, a_p, 1, 1);
        }
    }
}

void from_ntt_qprime(MatPoly &a) {
    assert(!a.isNTT);
    for (size_t r = 0; r < a.rows; r++) {
        for (size_t c = 0; c < a.cols; c++) {
            uint64_t *a_p   = &a.data[r * a.cols * poly_len + c * poly_len];
            ntt_qprime->ComputeInverse(a_p, a_p, 1, 1);
        }
    }
}

MatPoly mul_over_qprime(const MatPoly &a, const MatPoly &b) {
    assert(!a.isNTT);
    assert(!b.isNTT);
    assert(a.cols == b.rows);
    MatPoly res(a.rows, b.cols, false);

    int64_t *scratch = (int64_t *)calloc(poly_len, sizeof(int64_t));

    for (size_t r = 0; r < a.rows; r++) {
        for (size_t c = 0; c < b.cols; c++) {
            for (size_t m = 0; m < a.cols; m++) {
                uint64_t *a_p   = &a.data[r * a.cols * poly_len + m * poly_len];
                uint64_t *b_p   = &b.data[m * b.cols * poly_len + c * poly_len];
                uint64_t *res_p = &res.data[r * res.cols * poly_len + c * poly_len];
                for (size_t i = 0; i < poly_len; i++) {
                    scratch[i] = 0;
                }
                poly_mul_qprime(scratch, a_p, b_p);
                for (size_t i = 0; i < poly_len; i++) {
                    res_p[i] = (res_p[i] + scratch[i]) % arb_qprime;
                }
            }
        }
    }

    free(scratch);

    return res;
}

int64_t inv_mod(int64_t a, int64_t b) // https://rosettacode.org/wiki/Modular_inverse#C
{
	int64_t b0 = b, t, q;
	int64_t x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}
