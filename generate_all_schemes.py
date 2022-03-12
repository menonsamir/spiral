import itertools
import math
import pickle
import argparse
from multiprocessing import Pool
from timeit import default_timer as timer
import sys

p_err_bits = 40.0
real_q = 66974689739603969

n = 2
d = 2048
ks = ["p", "q", "t_GSW", "t_exp", "t_exp_right", "t_conv", "factor", "nu_1,nu_2"]
exp_lut_keys = set([(16, 6), (6, 9), (7, 3), (16, 9), (20, 7), (3, 7), (2, 5), (8, 5), (5, 8), (28, 5), (10, 8), (6, 7), (5, 5), (10, 7), (16, 3), (12, 6), (20, 9), (20, 4), (56, 9), (3, 2), (2, 6), (8, 2), (4, 5), (7, 5), (56, 6), (7, 8), (12, 8), (20, 3), (8, 9), (56, 3), (10, 3), (16, 7), (12, 2), (3, 6), (2, 2), (8, 6), (28, 6), (10, 9), (6, 4), (5, 4), (10, 4), (16, 4), (12, 7), (20, 5), (3, 5), (2, 7), (8, 3), (4, 6), (5, 7), (7, 4), (12, 4), (7, 6), (56, 7), (4, 8), (2, 8), (28, 8), (6, 2), (12, 9), (56, 4), (28, 2), (12, 3), (7, 7), (3, 9), (2, 3), (8, 7), (4, 2), (28, 7), (6, 5), (5, 3), (10, 5), (16, 5), (6, 8), (16, 8), (20, 6), (3, 4), (2, 4), (8, 4), (5, 9), (4, 7), (28, 4), (6, 6), (5, 6), (10, 6), (16, 2), (12, 5), (20, 8), (56, 8), (4, 9), (3, 3), (2, 9), (4, 4), (28, 9), (6, 3), (7, 2), (56, 5), (28, 3), (7, 9), (20, 2), (3, 8), (8, 8), (56, 2), (4, 3), (5, 2), (10, 2)])
def calc_fast(
    n=2,
    d=2048,
    p=12289,
    p_db=12289,
    q_prime=int(2**15.5),
    m_pt=1,
    q=2**56,
    sigma=6.4,
    t_conv=56,
    t_exp=56,
    t_exp_right=0,
    t_GSW=9,
    nu_1=7,
    nu_2=7,
    separate=False,
    du_first_dim=False,
    kinda_direct_upload=False,
    direct_upload=False,
    ternary=False,
    C=5
):
    z_GSW = math.ceil(q**(1.0/(t_GSW)))
    m_GSW = (n+1)*t_GSW
    z_exp = math.ceil(q**(1.0/(t_exp)))
    z_conv = math.ceil(q**(1.0/(t_conv)))
    B = (C * sigma)
    if ternary:
        B = 1
        assert(q <= 2**54)
    if t_exp_right == 0:
        t_exp_right = t_exp
    z_exp_right = math.ceil(q**(1.0/t_exp_right))

    num_exp_reg = nu_1 if separate else nu_1+1
    num_exp_reg += m_pt-1
    noise_scale_GSW = 4**(math.ceil(math.log(t_GSW*nu_2, 2))) if separate else 4*(t_GSW*nu_2+1)**2

    if du_first_dim:
        num_exp_reg = 0 if du_first_dim == True else (nu_1 - du_first_dim)
    sigma_hat_regev_2 = 4**(num_exp_reg)*sigma**2*(1 + d*t_exp*z_exp**2/3)
    if du_first_dim == True:
        sigma_hat_regev_2 = sigma**2 
    sigma_regev_2 = sigma_hat_regev_2 + d*t_conv*(z_conv**2)*(sigma**2)/4.0

    sigma_hat_GSW_2 = noise_scale_GSW*sigma**2*(1+t_exp_right*d*z_exp_right**2/3)
    if kinda_direct_upload:
        sigma_hat_GSW_2 = sigma**2
    sigma_GSW_2 = sigma_hat_GSW_2 * d * B**2 + t_conv*d*sigma**2*z_conv**2/2
    if direct_upload:
        sigma_GSW_2 = sigma**2

    sigma_0_2 = 2**nu_1*n*d*m_pt*(p_db**(1/m_pt)/2)**2*(sigma_regev_2) 
    sigma_rest = nu_2*d*m_GSW*z_GSW**2/2*(sigma_GSW_2)
    sigma_r_2 = sigma_0_2 + sigma_rest

    return sigma_r_2

# def get_p_err_fast(p, q_prime, q, s_e, sigma=6.4):
#     p = float(p)
#     q_prime = float(q_prime)
#     q = float(q)
#     s_e = float(s_e)
#     sigma = float(sigma)
    
#     q_mod_p = 1
#     q_prime_mod_p = 1
#     modswitch_adj = (1.0/2.0)*(q_prime_mod_p + (q_prime/q)*q_mod_p + 1)
#     thresh = (q_prime / (2.0*p)) - q_prime_mod_p - modswitch_adj
#     assert thresh > 0

#     s_round_2 = (sigma**2)*d/4
#     numer = float(-math.pi * (thresh)**2)
#     denom = float(s_e*(q_prime/float(q))**2 + (s_round_2))
#     p_single_err_log = math.log(2) + (numer / denom)
#     pr_err_log = p_single_err_log + math.log(2 * 2 * 2048)
#     return (pr_err_log)*math.log(math.e, 2)

def calc_fast_highrate(
    n=2,
    d=2048,
    p=12289,
    p_db=12289,
    q_prime=int(2**15.5),
    m_pt=1,
    q=2**56,
    sigma=6.4,
    t_conv=56,
    t_exp=56,
    t_exp_right=0,
    t_GSW=9,
    nu_1=7,
    nu_2=7,
    separate=False,
    du_first_dim=False,
    kinda_direct_upload=False,
    direct_upload=False,
    ternary=False,
    C=5
):
    true_n = n
    n = 1

    z_GSW = math.ceil(q**(1.0/(t_GSW)))
    m_GSW = (n+1)*t_GSW
    z_conv = math.ceil(q**(1.0/(t_conv)))
    z_exp = math.ceil(q**(1.0/(t_exp)))
    z_exp_right = math.ceil(q**(1.0/t_exp_right))

    num_exp_reg = nu_1+1

    sigma_regev_2 = sigma**2 
    sigma_GSW_2 = sigma**2

    if not kinda_direct_upload:
        noise_scale_GSW = 4**(math.ceil(math.log(t_GSW*nu_2, 2))+1)
        sigma_regev_2 = 4**(num_exp_reg)*sigma**2*(1 + d*t_exp*z_exp**2/3)
        sigma_GSW_2 = noise_scale_GSW*sigma**2*(1+t_exp_right*d*z_exp_right**2/3)
        sigma_GSW_2 = sigma_GSW_2 * d * (C*sigma)**2 + t_conv*d*sigma**2*z_conv**2/2

    sigma_0_2 = 2**nu_1*n*d*(p_db/2)**2*(sigma_regev_2) 
    sigma_rest = nu_2*d*m_GSW*z_GSW**2/2*(sigma_GSW_2)
    sigma_r_2 = sigma_0_2 + sigma_rest

    sigma_packing_2 = (d * true_n * (t_conv))*(sigma**2)*(z_conv**2)/4

    return sigma_r_2 + sigma_packing_2

p_mod_table = {
    17: 131072,
    18: 262144,
    19: 524288,
    20: 1048576,
    21: 2097152,
    22: 4194304,
    23: 8388592,
    24: 16777184,
    25: 33554332,
    26: 67108804,
    27: 134217608,
    28: 268435216,
    29: 536742296,
    30: 1073612276,
}
p_mod_table_pval = {(1 << key): value for (key, value) in p_mod_table.items()}
for i in range(1, 16+1):
    p_mod_table_pval[1<<i] = 1<<i
get_real_p = lambda x: p_mod_table_pval[x]

def get_p_err_fast_highrate(p, q_prime, q, s_e, n=2, sigma=6.4):
    old_p = p
    p = float(get_real_p(int(p)))
    q_prime = float(q_prime)
    q = float(q)
    s_e = float(s_e)
    sigma = float(sigma)

    q_mod_p = q % p
    # explanation:
    # we assume q_1 = 4*p
    # then the threshold for correctness is:
    # thresh = 1/2 - ((1/4)*((1/2)*((q_1 mod p) + (q_1/q)*(q mod p) + 2)))
    #        = 1/2 - (1/8)*(4*p/q)*(q mod p) - 1/4
    #        = 1/4 - (1/8)*((4*p)*q_mod_p/q)
    modswitch_adj = (1.0/8.0)*((4*p)*q_mod_p/q)
    thresh = (1/4) - modswitch_adj
    assert thresh > 0 and thresh <= (1.0/4.0), (thresh, p, old_p)

    s_round_2 = (sigma**2)*d/4
    numer = float(-math.pi * (thresh)**2)
    denom = float(s_e*(p/q)**2 + (s_round_2)*(p/q_prime)**2)

    p_single_err_log = math.log(2) + (numer / denom)
    pr_err_log = p_single_err_log + math.log(n * n * 2048)
    return (pr_err_log)*math.log(math.e, 2)

def simul(c_sel, streaming):
    c = dict(zip(ks, c_sel))
    nu_1 = c["nu_1,nu_2"][0]
    if not streaming and (c["t_exp"], nu_1) not in exp_lut_keys:
        return None
    
    qprime = c["p"]*(2**20)
    nu_2 = c["nu_1,nu_2"][1]
    
    q = c["q"]    
    params = {
        "p": c["p"], 
        "p_db": c["p"], 
        "q_prime": qprime,
        "q": c["q"], 
        "nu_1": nu_1, 
        "nu_2": nu_2, 
        "t_GSW": c["t_GSW"],
        "t_conv": c["t_conv"],
        "t_exp": c["t_exp"],
        "t_exp_right": c["t_exp_right"]
    }
    
    if streaming:
        params["kinda_direct_upload"] = True
        params["du_first_dim"] = True

    p_err = 9999
    s_e = calc_fast(**params)
    p_err = get_p_err_fast_highrate(c["p"], qprime, c["q"], s_e)
    if p_err > (-p_err_bits):
        return None
    
    qp_factor_bits = 8
    while qp_factor_bits <= 20:
        qprime = c["p"]*(2**qp_factor_bits) - c["p"] + 1
        params["q_prime"] = qprime
        
        s_e = calc_fast(**params)
        p_err = get_p_err_fast_highrate(c["p"], qprime, c["q"], s_e)
        if p_err <= (-p_err_bits):
            break
        qp_factor_bits += 0.1
    
    cost = 0
    c["p_err"] = p_err
    c["q_prime"] = params["q_prime"]
    c["s_e"] = math.log(s_e,2)
    return (cost, c, params)

def simul_highrate(c_sel, streaming):
    c = dict(zip(ks + ["n"], c_sel))
    nu_1 = c["nu_1,nu_2"][0]
    qprime = c["p"]*(2**20)
    nu_2 = c["nu_1,nu_2"][1]
        
    q = c["q"]    
    params = {
        "p": c["p"], 
        "p_db": c["p"], 
        "q_prime": qprime,
        "q": c["q"], 
        "nu_1": nu_1, 
        "nu_2": nu_2, 
        "t_GSW": c["t_GSW"],
        "t_conv": c["t_conv"],
        "t_exp": c["t_exp"],
        "t_exp_right": c["t_exp_right"],
        "n": c["n"]
    }

    if streaming:
        params["kinda_direct_upload"] = True
    
    p_err = 9999
    s_e = calc_fast_highrate(**params)
    p_err = get_p_err_fast_highrate(c["p"], qprime, c["q"], s_e, n=c["n"])
    if p_err > (-p_err_bits):
        return None
    # if nu_1 == 10 and c["p"] == 2**20 and c["n"]==4:
    #     print(params, math.log2(s_e), p_err)
    
    qp_factor_bits = 6
    while qp_factor_bits <= 20:
        qprime = c["p"]*(2**qp_factor_bits)
        params["q_prime"] = qprime
        
        s_e = calc_fast_highrate(**params)
        p_err = get_p_err_fast_highrate(c["p"], qprime, c["q"], s_e, n=c["n"])
        if p_err <= (-p_err_bits):
            break
        qp_factor_bits += 0.1

    cost = 0
    c["p_err"] = p_err
    c["q_prime"] = params["q_prime"]
    c["s_e"] = math.log(s_e,2)
    return (cost, c, params)

def simul_highrate_normal(c_sel):
    return simul_highrate(c_sel, False)

def simul_highrate_stream(c_sel):
    return simul_highrate(c_sel, True)

def simul_stream(c_sel):
    return simul(c_sel, True)

def simul_normal(c_sel):
    return simul(c_sel, False)

# L = []
# for t_exp in list(range(4, 56+1)):
#     L.append(simul(c_sel, False)
# return L

def get_regular_choices():
    poss_nus = []
    for j1 in range(2,10+1):
        for j2 in range(2,13+1):
            if j1+j2 < 10:
                continue
            poss_nus.append((j1,j2))
    choices = {
        "p": [2**i for i in range(2,15+1)],
        "q": [real_q],
        "t_GSW": list(range(2, 56+1)),
        "t_exp": [0],#list(range(4, 56+1)),
        "t_exp_right": [56],
        "t_conv": [2,4,8,16,32,56],
        "factor": [1],
        "nu_1,nu_2": poss_nus
    }
    do_streaming = False
    many_choices = []
    for t_exp in [2,4,8,16,32,56]:
        cop = choices.copy()
        cop["t_exp"] = [t_exp]
        many_choices.append(cop)
    return (many_choices, do_streaming, False)

# def get_regular_choices():
#     choices = {
#         "p": [512],
#         "q": [real_q],
#         "t_GSW": [8],
#         "t_exp": [16],
#         "t_exp_right": [56],
#         "t_conv": [4],
#         "factor": [1],
#         "nu_1,nu_2": [(9,7)]
#     }
#     do_streaming = False
#     return (choices, do_streaming)

def get_streaming_choices():
    poss_nus = []
    for j1 in range(2,14):
        for j2 in range(2,14):
            if j1+j2 < 10:
                continue
            poss_nus.append((j1,j2))
    choices = {
        "p": [2**i for i in range(2,20+1)],
        "q": [real_q],
        "t_GSW": list(range(2, 56+1)),
        "t_exp": [0],
        "t_exp_right": [56],
        "t_conv": [2,4,8,16,32,56],
        "factor": [1],
        "nu_1,nu_2": poss_nus
    }
    do_streaming = True
    many_choices = []
    for t_exp in [2,4,8,16,32,56]:
        cop = choices.copy()
        cop["t_exp"] = [t_exp]
        many_choices.append(cop)
    return (many_choices, do_streaming, False)

def get_highrate_choices():
    poss_nus = []
    for j1 in range(2,10+1):
        for j2 in range(2,14):
            if j1+j2 < 10:
                continue
            poss_nus.append((j1,j2))
    choices = {
        "p": [2**i for i in range(2,20+1)],
        "q": [real_q],
        "t_GSW": list(range(2, 56+1)),
        "t_exp": [2,4,8,16,32,56],
        "t_exp_right": [56],
        "t_conv": [2,4,8,16,32,56],
        "factor": [1],
        "nu_1,nu_2": poss_nus,
        "n": [0]
    }
    do_streaming = False
    many_choices = []
    for n in [2,4,8,12]:
        cop = choices.copy()
        cop["n"] = [n]
        many_choices.append(cop)
    return (many_choices, do_streaming, True)

def get_highrate_streaming_choices():
    poss_nus = []
    for j1 in range(2,14):
        for j2 in range(2,14):
            if j1+j2 < 10:
                continue
            poss_nus.append((j1,j2))
    choices = {
        "p": [2**i for i in range(10,30+1)],
        "q": [real_q],
        "t_GSW": list(range(2, 10+1)),
        "t_exp": [56],
        "t_exp_right": [56],
        "t_conv": [56],
        "factor": [1],
        "nu_1,nu_2": poss_nus,
        "n": [0]
    }
    do_streaming = True
    many_choices = []
    for n in [4,5,6,7,8,9,10,11,12]:
        cop = choices.copy()
        cop["n"] = [n]
        many_choices.append(cop)
    return (many_choices, do_streaming, True)

def clean(x):
    newitem = {}
    for k in x[2]:
        if k == 'q_prime':
            newitem['q_prime_bits'] = int(math.ceil(math.log(float(x[2][k]), 2)))
        elif k == 'q' or k == 'p_db':
            continue
        else:
            newitem[k] = int(x[2][k])
    newitem["s_e"] = x[1]["s_e"]
    return newitem

def prod(L):
    x = 1
    for i in L:
        x *= i
    return x

def perform_search(many_choices, do_streaming, do_high_rate):
    # perms = list(itertools.product(*[choices[k] for k in ks]))
    full_res = []
    for choices in many_choices:
        # print(".", len(choices))
        kss = ks
        if do_high_rate:
            kss = kss + ["n"]
        perms = itertools.product(*[choices[k] for k in kss])
        print(prod([len(i) for i in [choices[k] for k in kss]]))
        pool = Pool(processes=8)
        orig_res = []
        if do_high_rate:
            if do_streaming:
                # orig_res = list(map(simul_highrate_stream, perms))
                orig_res = pool.imap_unordered(simul_highrate_stream, perms, chunksize=256)
            else:
                # orig_res = list(map(simul_highrate_normal, perms))
                orig_res = pool.imap_unordered(simul_highrate_normal, perms, chunksize=256)
        elif do_streaming:
            orig_res = pool.imap_unordered(simul_stream, perms, chunksize=256)
        else:
            orig_res = pool.imap_unordered(simul_normal, perms, chunksize=256)
        pool.close()
        pool.join()

        res = list(map(clean, filter(lambda x: x is not None, orig_res)))
        full_res.extend(res)
    
    print("found", len(full_res))
    fname = 'all_params_streaming.pkl' if do_streaming else 'all_params.pkl'
    if do_high_rate:
        fname = 'all_params_highrate_streaming.pkl' if do_streaming else 'all_params_highrate.pkl'
    fh = open(fname, 'wb')
    pickle.dump(full_res, fh)
    fh.close()

def parse_all_args():
    parser = argparse.ArgumentParser(description='Generate scheme data.')
    parser.add_argument('--stream', help='use the streaming variant', action='store_true')
    parser.add_argument('--high-rate', help='use the high rate variant', action='store_true')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_all_args()
    
    result = None
    start = timer()
    if args.high_rate:
        if args.stream:
            perform_search(*get_highrate_streaming_choices())
        else:
            perform_search(*get_highrate_choices())
    else:
        if args.stream:
            perform_search(*get_streaming_choices())
        else:
            perform_search(*get_regular_choices())
    end = timer()
    print(f"Took {end - start} s") 