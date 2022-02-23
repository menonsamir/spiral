import pickle
import math
import os
import sys
import json
import pprint
import subprocess
import random
import re
import argparse

parser = argparse.ArgumentParser(description='Run Spiral.')
parser.add_argument('--stream',
                    help='measure performance for streaming', action='store_true')
parser.add_argument('--direct-upload',
                    help='perform direct uploading (no query expansion)', action='store_true')
parser.add_argument('--high-rate', '--pack',
                    help='use the packing variant', action='store_true')
parser.add_argument('--optimize-for',
                    help='optimize for a certain property', type=str, default="")
parser.add_argument('--set-dims',
                    help='set nu_1 and nu_2', type=str, default="")
parser.add_argument('--show-output',
                    help='show output', action='store_true')
parser.add_argument('--show-top',
                    help='show other top candidates', action='store_true')
parser.add_argument('--quiet',
                    help='quiet output', action='store_true')
parser.add_argument('--corr', '--explicit-db',
                    help='check correctness using an explicit database representation', action='store_true')
parser.add_argument('--skip-cmake',
                    help='don\'t use cmake', action='store_true')
parser.add_argument('--vcpkg-root',
                    help='root directory of vcpkg', type=str, default="")
parser.add_argument('--dry-run',
                    help='don\'t run the scheme', action='store_true')
parser.add_argument('--skip-make',
                    help='don\'t re-make', action='store_true')
parser.add_argument('--addtl-flags',
                    help='add additional flags to make')
parser.add_argument('--max-query-size',
                    help='maximum acceptable query size, in bytes', type=int, default=0)
parser.add_argument('--max-param-size',
                    help='maximum acceptable public parameter size, in bytes', type=int, default=0)
parser.add_argument('--max-total-query-size',
                    help='maximum acceptable query + public parameter size, in bytes', type=int, default=0)
parser.add_argument('--ignore-expansion',
                    help='ignore query expansion costs', action='store_true')
parser.add_argument('--analyze-deviation',
                    help='compute absolute and relative errors in cost model', action='store_true')
parser.add_argument('--build-exp-lut',
                    help='build the lookup table for expansion', action='store_true')
parser.add_argument('--build-fdim-lut',
                    help='build the lookup table for expansion', action='store_true')
parser.add_argument('--reduce-params',
                    help='reduce the parameter space', action='store_true')
parser.add_argument('targetnum', metavar='logN', type=int,
                    help='log2 of the number of items')
parser.add_argument('itemsize', metavar='itemsize', type=int,
                    help='item size in bytes')
parser.add_argument('--trials', metavar='trials', type=int, nargs='?', const=1, default=1,
                    help='number of trials to run over')
args = parser.parse_args()

use_cmake = not args.skip_cmake
toolchain_path = ""
if use_cmake:
    vcpkg_root = args.vcpkg_root
    if not vcpkg_root or len(vcpkg_root) == 0:
        vcpkg_root = os.path.join(os.getenv("HOME"), "vcpkg")
    rel_toolchain_path = "scripts/buildsystems/vcpkg.cmake"
    toolchain_path = os.path.join(vcpkg_root, rel_toolchain_path)

targetnum = 1 << args.targetnum
itemsize = args.itemsize
streaming = args.stream
high_rate = args.high_rate 
show_output = args.show_output
direct_upload = args.direct_upload 
addtl_flags = args.addtl_flags or ""
exp_lut_fname = 'exp_lut.json' if not high_rate else 'exp_lut_highrate.json'
fdim_lut_fname = 'fdim_lut.json' if not high_rate else 'fdim_lut_highrate.json'
optimize_for = args.optimize_for

exp_lut = {}
if not args.build_exp_lut:
    fh = open(exp_lut_fname, 'r')
    obj = json.load(fh)
    fh.close()
    exp_lut = obj

fdim_lut = {}
if high_rate and direct_upload:
    fh = open(fdim_lut_fname, 'r')
    fdim_lut = json.load(fh)
    fh.close()

usd_per_us = 5.41666667e-12
usd_per_byte = 9e-11 #0.054 / (3600000000)
min_q_prime_bits = 14

d = 2048
logq = 56

def get_query_size(params):
    if direct_upload:
        if high_rate:
            return (2**(params["nu_1"]) + 2*params["t_GSW"]*params["nu_2"]) * d*logq/8
        else:
            return (2**(params["nu_1"]) + params["t_GSW"]*params["nu_2"]) * d*logq/8
    else:
        return d*logq/8

def get_pp_sz_here(params):
    if args.high_rate:
        packing_sz = params["n"]*params["t_conv"]*(params["n"]+1)
        exp_mats = params["nu_1"]*2*(params["t_exp"]) + int(math.ceil(math.log2(params["nu_1"]*params["t_GSW"])))*2*(params["t_exp_right"])
        v_sz = 2 * 2 * params["t_conv"]
        total_elems = packing_sz
        if not args.direct_upload:
            total_elems += exp_mats + v_sz
        return int(total_elems * d * logq / 8)
    else:
        exp_mats = params["nu_1"]*2*(params["t_exp"]) + int(math.ceil(math.log2(params["nu_1"]*params["t_GSW"])))*2*(params["t_exp_right"])
        w_sz = 3 * 2 * params["t_conv"]
        v_sz = 3 * 2 * params["t_conv"]
        total_elems = w_sz + v_sz
        if not args.direct_upload:
            total_elems += exp_mats
        return int(total_elems * d * logq / 8)

def calc_cost(
    n=2,
    d=2048,
    p=12289,
    p_db=12289,
    q_prime_bits=16,
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
    print_params=False,
    experimental=False,
    ternary=False,
    C=5,
    factor=1,
    s_e=None
):    
    folding_time_us = 1000*(33 + 29.6*t_GSW)*(2**nu_2 / 2**6)
    
    lr_fd = (0, 619.13591337, 9.25842148) # from linear regression on 2**nu_1, 2**nu_2, 2**(nu_1+nu_2)
    firstdim_time_us = lr_fd[0]*(2**nu_1) + lr_fd[1]*(2**nu_2) + lr_fd[2]*(2**(nu_1+nu_2))

    exp_us = 0
    t_conv_for_4 = 4
    comp_us = 185451 * (2**nu_1 / 2**9) * (t_conv / t_conv_for_4)
    conv_us = 93709 * (nu_2 * t_GSW / (40)) * (t_conv / t_conv_for_4) #??
    if (kinda_direct_upload and du_first_dim):
        exp_us = 0
    else:
        exp_us = exp_lut[str((nu_1, 6, t_exp))][u'exp_us']#[u'exp_specific_us']

    if args.ignore_expansion:
        exp_us = 0
        comp_us = 0
        conv_us = 0

    total_us = exp_us + comp_us + conv_us + factor*(firstdim_time_us + folding_time_us)

    if q_prime_bits < min_q_prime_bits:
        q_prime_bits = min_q_prime_bits

    total_bytes = factor * ((2*2*2048*(math.ceil(math.log2(4*p))) + 2*2048*(q_prime_bits))/8)
    
    total_cost = total_us*usd_per_us + usd_per_byte*total_bytes
    if streaming:
        total_us = firstdim_time_us + folding_time_us
        total_cost = total_us*usd_per_us + usd_per_byte*total_bytes
        total_cost = total_cost / (4*2048*math.log(p, 2)/8.0)
    
    if streaming:
        return (total_cost, total_us, total_bytes)
    else:
        return (total_cost, total_us, total_bytes, (exp_us, comp_us, conv_us, firstdim_time_us, folding_time_us))

def calc_cost_highrate(
    n=2,
    d=2048,
    p=12289,
    p_db=12289,
    q_prime_bits=16,
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
    print_params=False,
    experimental=False,
    ternary=False,
    C=5,
    factor=1,
    s_e=None
):    
    folding_time_per_us = 0.5*1000*(11.7 + 18.9*t_GSW)*(2**nu_2 / 2**6)
    folding_time_us = factor * (n * n) * folding_time_per_us
    
    if str((nu_1, nu_2)) in fdim_lut:
        firstdim_time_per_us = fdim_lut[str((nu_1, nu_2))]["fdim_us"]/4
    else:
        lr_fd = (0, 73.41112501, 1.985572062) # from linear regression on 2**nu_1, 2**nu_2, 2**(nu_1+nu_2)
        firstdim_time_per_us = lr_fd[0]*(2**nu_1) + lr_fd[1]*(2**nu_2) + lr_fd[2]*(2**(nu_1+nu_2))

    firstdim_time_us = factor * (n * n) * firstdim_time_per_us

    pack_us = factor * (3230) * (n * n)

    exp_us = 0
    comp_us = 0
    conv_us = 0

    if kinda_direct_upload or streaming:
        exp_us = 0
    else:
        exp_us = exp_lut[str((nu_1, 6, t_exp))][u'exp_us']#[u'exp_specific_us']

    total_us = exp_us + comp_us + conv_us + firstdim_time_us + folding_time_us + pack_us

    if q_prime_bits < min_q_prime_bits:
        q_prime_bits = min_q_prime_bits
    total_bytes = factor * ((n*n*2048*(math.ceil(math.log2(4*p))) + n*2048*(q_prime_bits))/8)
    
    total_cost = total_us*usd_per_us + usd_per_byte*total_bytes
    if streaming:
        total_us = firstdim_time_us + folding_time_us + pack_us
        total_cost = total_us*usd_per_us + usd_per_byte*total_bytes
        total_cost = total_cost / (n*n*2048*math.log(p, 2)/8.0)
    
    return (total_cost, total_us, total_bytes, (exp_us, comp_us, conv_us, pack_us/factor, firstdim_time_us/factor, folding_time_us/factor))

def shortcut_optim(cost_res, itemsize, targetnum):
    cost = 0
    if optimize_for == "rate":
        cost = -itemsize / cost_res[2]
    elif optimize_for == "tput":
        cost = -(itemsize * targetnum) / (cost_res[1])
    else:
        assert False
    return (cost, cost_res)

all_ks = ['p', 'q_prime_bits', 'nu_1', 'nu_2', 't_GSW', 't_conv', 't_exp', 't_exp_right']
if high_rate:
    all_ks.append("n")
def apply_factor(itemsize, targetnum, x):
    if type(x) == tuple:
        x = dict(zip(all_ks, x))
    n = 2
    if high_rate:
        n = x["n"]
    base_item_size = float(n*n*2048*math.log(x["p"],2)/8)
    fact = math.ceil(itemsize / base_item_size)
    c_fn = calc_cost_highrate if high_rate else calc_cost
    cost_res = c_fn(factor=fact, **x)
    if len(optimize_for) > 0:
        cost_res = shortcut_optim(cost_res, itemsize, targetnum)
    return (cost_res, fact, x) 

def pred(itemsize, targetnum, x):
    if type(x) == tuple:
        x = dict(zip(all_ks, x))
    n = 2
    if high_rate:
        n = x["n"]
    if args.set_dims:
        parts = args.set_dims.split(",")
        req_nu_1 = int(parts[0])
        req_nu_2 = int(parts[1])
        if x["nu_1"] != req_nu_1 or x["nu_2"] != req_nu_2:
            return False
    base_item_size = float(n*n*2048*math.log(x["p"],2)/8)
    fact = math.ceil(itemsize / base_item_size)

    mindbsize = targetnum * itemsize
    dbsz = fact * base_item_size * float(2**(x["nu_1"]+x["nu_2"]))
    p1 = dbsz >= mindbsize
    if streaming:
        p1 = p1 and (2**(x["nu_1"]+x["nu_2"]) >= targetnum) 
    if args.max_query_size > 0:
        p1 = p1 and get_query_size(x) <= args.max_query_size
    if args.max_param_size > 0:
        p1 = p1 and get_pp_sz_here(x) <= args.max_param_size
    if args.max_total_query_size > 0:
        p1 = p1 and (get_query_size(x)+get_pp_sz_here(x)) <= args.max_total_query_size
    if high_rate and not streaming:
        p1 = p1 and (str((x["nu_1"], 6, x["t_exp"])) in exp_lut) 
    if high_rate and direct_upload:
        p1 = p1 and (str((x["nu_1"], x["nu_2"])) in fdim_lut)
    return p1

param_f = lambda t: f"TEXP={t[0]} TEXPRIGHT={t[1]} TCONV={t[2]} TGSW={t[3]} QPBITS={t[4]} PVALUE={t[5]} QNUMFIRST={t[6]} QNUMREST={t[7]} OUTN={t[8]} ADDTL_CXX_FLAGS={addtl_flags}"
cmd_mk = lambda t: f"make clean && PARAMSET=PARAMS_DYNAMIC " + param_f(t)+ " make spiral -j4"
if use_cmake:
    cmd_clean = "rm -rf build/"
    cmd_cmake = f"cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE={toolchain_path}"

    # hack for compatibility
    if 'avx2' in addtl_flags:
        cmd_cmake += " -DNOAVX2=1"
    if 'avx512' in addtl_flags:
        cmd_cmake += " -DNOAVX512=1"
    if 'CRT' in addtl_flags:
        cmd_cmake += " -DNOCRT=1"
    
    cmd_mk_pref = "cmake --build build -j4 -- PARAMSET=PARAMS_DYNAMIC "
    cmd_copy = "cp build/spiral ./spiral"
    cmd_mk = lambda t: " && ".join([cmd_clean, cmd_cmake, cmd_mk_pref + param_f(t), cmd_copy])

def make_for(params):
    ksp = ["t_exp", "t_exp_right", "t_conv", "t_GSW", "q_prime_bits", "p"]
    vals = [params[k] for k in ksp]
    if not direct_upload:
        vals = vals + [1, 0]
    else:
        vals = vals + [1<<params["nu_1"], params["t_GSW"]*params["nu_2"]]
    if high_rate:
        vals = vals + [params["n"]]
    else:
        vals = vals + [2]
    cmd = cmd_mk(vals)
    if not args.quiet:
        print(cmd)
    s = subprocess.check_output(cmd, shell=True)
    return s

command_base = "./spiral"
command_opts_incorr = "a --random-data"
command_opts_corr = "a"

def run_spiral(nu_1, nu_2, opts):
    total_N = 2**(nu_1 + nu_2)
    ridx = random.randint(0, total_N-1)
    cmd_str = f"{command_base} {nu_1} {nu_2} {ridx} {opts}"
    if high_rate:
        cmd_str += " --high-rate"
    s = subprocess.check_output(cmd_str, shell=True)
    return s.decode('utf8')

def analyze_spiral(s, factor=1):
    exp_re = r"\s+Main expansion.*:\s+([0-9]+)"
    exp_specific_re = r"\s+Expansion took.*:\s+([0-9e\+\.]+)"
    conv_re = r"\s+Conversion.*:\s+([0-9]+)"
    scaltomat_re = r"\s+ScalToMat took.*:\s+([0-9]+)"
    regtogsw_re = r"\s+RegevToGSW took.*:\s+([0-9]+)"
    fdim_re = r"\s+First dimension multiply.*:\s+([0-9]+)"
    fold_re = r"\s+Folding.*:\s+([0-9]+)"
    pack_re = r"\s+Packing.*:\s+([0-9]+)"
    query_gen_re = r"\s+Query generation.*:\s+([0-9]+)"
    key_gen_re = r"\s+Key generation.*:\s+([0-9]+)"
    decoding_re = r"\s+Decoding.*:\s+([0-9]+)"
    
    resp_sz_re = r"\s+Response size.*:\s+([0-9]+)"
    query_sz_re = r"\s+online query size.*:\s+([0-9]+)"
    off_query_sz_re = r"\s+offline query size.*:\s+([0-9]+)"
    corr_re = r"\s+Is correct?.*:\s+([0-9])"
    
    
    exp_us = int(re.search(exp_re, s).group(1))
    exp_specific_us = int(float(re.search(exp_specific_re, s).group(1)))
    conv_us = int(re.search(conv_re, s).group(1))
    scaltomat_us = int(re.search(scaltomat_re, s).group(1))
    regtogsw_us = int(re.search(regtogsw_re, s).group(1))
    fdim_us = factor*int(re.search(fdim_re, s).group(1))
    fold_us = factor*int(re.search(fold_re, s).group(1))
    pack_us = 0
    if high_rate:
        pack_us = factor*int(re.search(pack_re, s).group(1))
    
    query_gen_us = int(re.search(query_gen_re, s).group(1))
    key_gen_us = int(re.search(key_gen_re, s).group(1))
    decoding_us = int(re.search(decoding_re, s).group(1))
    
    resp_sz_b = factor*int(re.search(resp_sz_re, s).group(1))
    query_sz_b = int(re.search(query_sz_re, s).group(1))
    off_query_sz_b = int(re.search(off_query_sz_re, s).group(1))
    if not high_rate:
        query_sz_b = query_sz_b / 2 # to use seed trick
    is_corr = int(re.search(corr_re, s).group(1)) == 1
    
    total_us = exp_us + conv_us + fdim_us + fold_us + pack_us
    return {
        "exp_us": exp_us,
        "exp_specific_us": exp_specific_us,
        "conv_us": conv_us,
        "scaltomat_us": scaltomat_us,
        "regtogsw_us": regtogsw_us,
        "fdim_us": fdim_us,
        "fold_us": fold_us,
        "pack_us": pack_us,
        "total_us": total_us,
        "query_gen_us": query_gen_us,
        "key_gen_us": key_gen_us,
        "decoding_us": decoding_us,
        "resp_sz": resp_sz_b,
        "query_sz": query_sz_b,
        "param_sz": off_query_sz_b,
        "is_corr": is_corr
    }
def run_and_analyze_spiral(nu_1=9, nu_2=6, factor=1, opts=command_opts_corr):
    s = run_spiral(nu_1, nu_2, opts)
    if show_output:
        print(s)
    return analyze_spiral(s, factor)

if args.build_exp_lut:
    base_params = {
        "p": 1024,
        "q_prime_bits": 25,
        "nu_1": 0,
        "nu_2": 0,
        "t_GSW": 10,
        "t_conv": 8,
        "t_exp": 0,
        "t_exp_right": 56,
        "query_size": 14336.0
    }
    if high_rate:
        base_params["n"] = 4
    poss_t_exp = [2,4,8,16,32,56]
    exp_table = {}
    nu_2 = 6
    for nu_1 in range(6,12+1):
        for t_exp in poss_t_exp:
            base_params["nu_1"] = nu_1
            base_params["nu_2"] = nu_2
            base_params["t_exp"] = t_exp
            make_for(base_params)
            run_s = run_and_analyze_spiral(nu_1, nu_2, 1, command_opts_incorr)
            exp_table[str((nu_1, nu_2, t_exp))] = run_s
            print(".")
    fh = open(exp_lut_fname, 'w')
    obj = json.dump(exp_table, fh)
    fh.close()
    sys.exit(0)

if args.build_fdim_lut:
    base_params = {
        "kinda_direct_upload": True,
        "n": 2,
        "p": 65536,
        "q_prime_bits": 27,
        "nu_1": 0,
        "nu_2": 0,
        "t_GSW": 3,
        "t_conv": 56,
        "t_exp": 56,
        "t_exp_right": 56,
        "query_size": 58892288.0
    }
    # for nu_1 in range(2,9+1):
    #     for nu_2 in range(2,11+1):
    fdim_table = {}
    for nu_1 in range(5,12+1):
        for nu_2 in range(2,12+1):
            if (nu_1+nu_2 < 10) or (nu_1+nu_2 > 20):
                continue
            base_params["nu_1"] = nu_1
            base_params["nu_2"] = nu_2
            make_for(base_params)
            try:
                run_s = run_and_analyze_spiral(nu_1, nu_2, 1, command_opts_incorr)
            except KeyboardInterrupt:
                sys.exit(0)
            except:
                print("failed for", nu_1, nu_2)
                continue
            fdim_table[str((nu_1, nu_2))] = run_s
            print(".")
    fh = open(fdim_lut_fname, 'w')
    obj = json.dump(fdim_table, fh)
    fh.close()
    sys.exit(0)

if streaming and not args.quiet:
    print("streaming...")
    itemsize = 1

fname = 'all_params.pkl' if not direct_upload else 'all_params_streaming.pkl'
if high_rate:
    fname = 'all_params_highrate.pkl' if not direct_upload else 'all_params_highrate_streaming.pkl'
fh = open(fname, 'rb')
orig_res = pickle.load(fh)
fh.close()

if not args.quiet:
    print('Loaded parameter space.')

keyf = lambda x: x[0]

def apply_factor_curry(x):
    return apply_factor(itemsize, targetnum, x)

res = list(filter(lambda x: pred(itemsize, targetnum, x), orig_res))
res = list(map(apply_factor_curry, res))
res.sort(key=keyf)
paramset = res[0]

if not args.quiet:
    print('Found parameters:')
    pprint.pprint(paramset)

def get_rate(r):
    return itemsize / (r[0][2])
def get_tput(r):
    return (itemsize * targetnum) / (r[0][1])

factor = paramset[1]
params = paramset[2]
if params["q_prime_bits"] < min_q_prime_bits:
    params["q_prime_bits"] = min_q_prime_bits
params["query_size"] = get_query_size(params)

if not args.skip_make:
    make_output = make_for(params)
    # print(make_output)

if args.dry_run:
    sys.exit(0)

opts = command_opts_corr if args.corr else command_opts_incorr
runs_l = []
for i in range(args.trials):
    run_s = run_and_analyze_spiral(params["nu_1"], params["nu_2"], factor, opts=opts)
    runs_l.append(run_s)
    if not args.quiet:
        pprint.pprint(run_s)

avgd_runs = {k: sum([r[k] for r in runs_l])/len(runs_l) for k in runs_l[0].keys()}
this_n = 2
if high_rate:
    this_n = params["n"]
avgd_runs["item_sz"] = factor * (this_n*this_n*2048*math.log(params["p"], 2)/8)
avgd_runs["dbsize"] = factor * (2**(params["nu_1"] + params["nu_2"]))*(this_n*this_n*2048*math.log(params["p"], 2)/8)
avgd_runs["params"] = params

# note: this calculation is throughput *excluding* expansion costs; only useful for streaming case
avgd_runs["tput"] = (avgd_runs["dbsize"]) / (avgd_runs["fdim_us"] + avgd_runs["fold_us"] + avgd_runs["pack_us"])

avgd_runs["rate"] = itemsize / avgd_runs["resp_sz"]
avgd_runs["cost"] = usd_per_us * avgd_runs["total_us"] + usd_per_byte * avgd_runs["resp_sz"]
print(json.dumps(avgd_runs))

if args.analyze_deviation:
    print("factor", factor)
    if len(optimize_for) > 0:
        paramset = list(paramset)
        paramset[0] = paramset[0][1]
    if high_rate:
        keys = ["exp_us", "scaltomat_us", "regtogsw_us", "pack_us", "fdim_us", "fold_us", "total_us"]
    else:
        keys = ["exp_us", "scaltomat_us", "regtogsw_us", "fdim_us", "fold_us", "total_us"]
    predicted_times = list(paramset[0][3])
    if high_rate:
        predicted_times[-3] *= factor
    predicted_times[-2] *= factor
    predicted_times[-1] *= factor
    predicted_times = predicted_times + [sum(predicted_times)]
    actual_times = [avgd_runs[k] for k in keys]

    abs_err = [abs(a-b) for a,b in zip(predicted_times, actual_times)]
    rel_err = [abs(a-b)/float(a) if a != 0 else 0 for a,b in zip(predicted_times, actual_times)]

    table = [
        ["predicted_times"] + predicted_times,
        ["actual_times"] + actual_times,
        ["abs_err"] + abs_err,
        ["rel_err"] + rel_err
    ]
    for row in [["name"] + keys] + table:
        print(" ".join(list(map(str, row))))
