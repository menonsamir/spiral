import pickle
import os
import sys
import subprocess
from tabulate import tabulate
from run_scheme import run_system_tr
import argparse
import math
from util import get_pp_size
from timeit import default_timer as timer

systems_table     = ["sealpir", "fastpir", "onionpir", "spiral", "spiralstream", "nopriv"]
headers_table     = ["SealPIR", "FastPIR", "MulPIR", "OnionPIR", "Spiral", "SpiralStream", "NoPriv"]

systems_asympcomp = ["spiral", "spiral-pack", "sealpir", "fastpir", "onionpir"]
headers_asympcomp = ["Spiral", "Spiral-Pack", "SealPIR", "FastPIR", "OnionPIR"]
scenarios_all_asympcomp = {
    "asympcomp": [(i, 10000) for i in range(10, 20+1, 2)],
    "asympcomplarge": [(i, 100000) for i in range(10, 20+1, 2)]
}

# excluding "spiralstreambar"
scenarios_streaming = list(range(10, 20+1, 2))
scenarios_streaming_disp = list(range(12, 20+1, 4))
systems_streaming = ["spiral", "spiralstream", "spiralstreambar", "spiral-pack", "spiralstream-pack", "sealpir", "fastpir", "onionpir"]
headers_streaming = ["Spiral", "SpiralStream", "SpiralStreamBar", "Spiral-Pack", "SpiralStream-Pack", "SealPIR", "FastPIR", "OnionPIR"]

scenarios_table = [
    (20, 256),
    (18, 30000),
    (14, 100000)
]

scenarios_ubench = [(i, 100000) for i in range(10, 20+1, 2)]
headers_ubench    = ["KeyGen", "QueryGen", "Expansion", "FirstDim", "Folding", "Decoding"]
keys_ubench = ["key_gen_us", "query_gen_us", "exp_us", "fdim_us", "fold_us", "decoding_us"]

scenarios_ablation = [
    (20, 256),
    (18, 30000)
]

systems_packingcomp = ["best",          "spiral", "spiralstream", "spiral-pack", "spiralstream-pack"]
headers_packingcomp = ["Best previous", "Spiral", "SpiralStream", "Spiral-Pack", "SpiralStream-Pack"]

scenarios_packingcomp = [
    (20, 256),
    (18, 30000),
    (14, 100000)
]
best_sys_for_scenario_pc = {
    (20, 256): "fastpir",
    (18, 30000): "onionpir",
    (14, 100000): "onionpir"
}

scenarios_limits = [
    (20, 256),
    (18, 30000),
    (14, 1000000)
]
systems_limits = ["spiralstream", "spiralstream-pack"]
keys_limits = ["rate", "tput", "param_sz", "query_sz", "resp_sz", "total_us"]
keys_limits_disp = ["rate", "tput", "param_sz", "query_sz"]

max_query_sizes_mb = [1, 2, 5, 10, 20, 30, 40, 50, 60, 70]#[1, 5, 10, 20, 30]
max_predicates = ["query", "param", "total-query"]
scenarios_maxtotalquery = [(14, 100000, line, int(i*1000000)) for i in max_query_sizes_mb for line in max_predicates]
systems_maxtotalquery = ["spiral", "spiral-pack", "spiralstream", "spiralstream-pack"]

usd_per_us = 5.41666667e-12
usd_per_byte = 9e-11 #0.054 / (3600000000)
download_speed_B_s = 3576250 # 28.61 Mbps
upload_speed_B_s   = 1047500 #  8.38 Mbps

# def gen_table(, trials=1):
#     scenarios = [
#         (20, 256),
#         (, 256),
#     ]

def save_file(figurename, results):
    fh = open(figurename + '_results.pkl', 'wb')
    pickle.dump(results, fh)
    fh.close()

def load_file(figurename, path=""):
    fname = figurename + '_results.pkl'
    if path and len(path) > 0:
        fname = os.path.join(path, fname)
    fh = open(fname, 'rb')
    result = pickle.load(fh)
    fh.close()
    return result

# Calculations
def get_total_size(scenario):
    return (2**(scenario[0]))* scenario[1]

def check_corr(res):
    if "is_corr" in res:
        return res["is_corr"]
    return True

# Formatting
def get_key_tex(key):
    key_to_tex = {
        "param_sz": "{\\bf Param. Size}",
        "query_sz": "{\\bf Query Size}",
        "resp_sz": "{\\bf Response Size}",
        "total_us": "{\\bf Computation}",
        "rate": "{\\bf Rate}",
        "tput": "{\\bf Throughput}",
        "cost": "{\\bf Server Cost}"
    }
    return key_to_tex[key]

def get_div(sz):
    if sz < 1000:
        return 1
    elif sz < 1000000:
        return 1e3
    elif sz < 1e9:
        return 1e6
    else:
        return 1e9

def get_suf(sz):
    if sz < 1000:
        return "B"
    elif sz < 1000000:
        return "KB"
    elif sz < 1e9:
        return "MB"
    else:
        return "GB"

def get_sz_str(sz, fmt="{:.2f}", force_fmt=False):
    suf = get_suf(sz)
    div = get_div(sz)
    s = str(sz)
    if force_fmt or sz >= 1e9:
        s = fmt.format(sz/div)
    else:
        s = str(int(round(sz/div)))
    return (s, suf)

def get_pp_row(results, scenarios, systems):
    row = ["", get_key_tex("param_sz")]
    param_szs = [max([get_pp_size(system, results[str(scenario)][system]) for scenario in scenarios]) for system in systems]
    min_sz = min(param_szs)
    best_pp_system = systems[param_szs.index(min_sz)]

    for system in systems:
        sz_str = None
        if "spiral" in system:
            szs_over_all_scenarios = [get_pp_size(system, results[str(scenario)][system]) for scenario in scenarios]
            min_sz_str = get_sz_str(min(szs_over_all_scenarios), fmt="{:.1f}")
            max_sz_str = get_sz_str(max(szs_over_all_scenarios), fmt="{:.1f}")
            if min_sz_str[1] == max_sz_str[1]:
                sz_str = min_sz_str[0]+"--"+max_sz_str[0]+" "+min_sz_str[1]
            else:
                sz_str = min_sz_str[0]+" "+min_sz_str[1]+"--"+max_sz_str[0]+" "+max_sz_str[1]
        else:
            sz = get_pp_size(system)
            sz_str = " ".join(get_sz_str(sz, fmt="{:.1f}"))
            if system == best_pp_system:
                sz_str = "\\emphcell{"+sz_str+"}"
        row.append(sz_str)
    return row

def get_db_tex(scenario):
    if type(scenario) == tuple:
        N = str(scenario[0])
        num, suf = get_sz_str(scenario[1])
        return "$\mathbf{2^{"+N+"} \\times "+num+"} \\textbf{"+suf+"}$"
    else:
        N = str(scenario)
        return "$\mathbf{2^{"+N+"}}$"

def get_row_start(figurename, scenario, key):
    k_db, k_tot_sz = ("total_us", "rate") if figurename == "table" else ("resp_sz", "total_us")
    if figurename == "streaming":
        k_db, k_tot_sz = "query_sz", None
    if key == k_db:
        return get_db_tex(scenario)
    elif key == k_tot_sz:
        num, suf = get_sz_str(get_total_size(scenario), fmt="{:.1f}")
        return "{\\bf ("+num+" "+suf+")}"
    return ""

def metric_str_output(key, value, fmt="{:.2f}", force_fmt=False):
    if key in ["param_sz", "query_sz", "resp_sz"]:
        return " ".join(get_sz_str(value, fmt=fmt, force_fmt=force_fmt))
    elif key == "total_us":
        return "{:.2f} s".format(value/1000000)
    elif key == "rate":
        return "{:.4f}".format(value)
    elif key == "tput":
        return " ".join(get_sz_str(value*1e6)) + "/s"
    elif key == "cost":
        return "\\${:.6f}".format(value)
    else:
        return str(value)

def finish_output_str(figurename, table, headers):
    s = tabulate(table, headers, tablefmt="latex_raw", disable_numparse=True)
    lines = s.split("\n")
    lines_out = []
    for i in range(len(lines)):
        suf = ""
        if figurename == "table":
            if "Public" in lines[i]:
                suf = " \\midrule"
            elif "Computation" in lines[i]:
                suf = " \\cmidrule{2-8}"
            elif "Cost" in lines[i] and i != len(lines) - 3:
                suf = " \\midrule"
        elif figurename == "packingcomp":
            if "Computation" in lines[i]:
                suf = " \\cmidrule{2-7}"
            elif "Throughput" in lines[i] and i != len(lines) - 3:
                suf = " \\midrule"
        elif figurename == "limits" and i >= 4:
            idx = i - 4
            if idx % 2 == 1 and "times" in lines[i] and i != len(lines) - 3:
                suf = " \\midrule"
        elif figurename == "streaming":
            if "Throughput" in lines[i] and i != len(lines) - 3:
                suf = " \\midrule"
        lines[i] = lines[i] + suf
    return "\n".join(lines)

obj_fns = {
    "param_sz": min,
    "query_sz": min,
    "resp_sz": min,
    "total_us": min,
    "cost": min,
    "rate": max,
    "tput": max
}

def to_num(x):
    try:
        return float(x)
    except ValueError:
        return math.nan

def emph_str(s):
    return "\\emphcell{"+s+"}"

def emph_data(key, raw_data):
    # print(key, raw_data)
    raw_data = list(map(to_num, raw_data))
    best_val = obj_fns[key](raw_data)

    idxs_of_best = [i for i in range(len(raw_data)) if raw_data[i]==best_val]
    outs = [metric_str_output(key, d) for d in raw_data]
    for idx_of_best in idxs_of_best:
        outs[idx_of_best] = emph_str(outs[idx_of_best])
    return outs

def insert_blank(row, idx_ins):
    return row[:idx_ins] + ["-"] + row[idx_ins:]

# Main benchmarking + display
def gen_asympcomp(figurename, trials=1):
    results = {}
    scenarios = scenarios_all_asympcomp[figurename]
    for scenario in scenarios:
        results[str(scenario)] = {}
        for system in systems_asympcomp:
            result = run_system_tr(system, scenario[0], scenario[1], False, trials=trials)
            results[str(scenario)][system] = result
    save_file(figurename, results)
    return results

def display_asympcomp(figurename, results):
    scenarios = scenarios_all_asympcomp[figurename]
    headers = ["Database"] + headers_asympcomp
    table = [[1<<scenario[0]] + [int(results[str(scenario)][system]["total_us"]) for system in systems_asympcomp] for scenario in scenarios]
    print(tabulate(table, headers, tablefmt="plain"))
    assert all([check_corr(results[str(scenario)][system]) for system in systems_asympcomp for scenario in scenarios])

def gen_streaming(figurename, trials=1):
    results = {}
    scenarios = scenarios_streaming
    for scenario in scenarios:
        results[str(scenario)] = {}
        for system in systems_streaming:
            if system == "spiralstreambar":
                result = run_system_tr("spiralstream", scenario, 1, True, cmd_extras="--max-query-size 16000000", trials=trials)
                results[str(scenario)][system] = result
            elif system == "spiralstream":
                result = run_system_tr("spiralstream", scenario, 1, True, cmd_extras="--max-query-size 33000000", trials=trials)
                results[str(scenario)][system] = result
            else:
                result = run_system_tr(system, scenario, 1, True, cmd_extras="--max-query-size 33000000", trials=trials)
                results[str(scenario)][system] = result
    save_file(figurename, results)
    return results

def display_streaming_graph(figurename, results):
    # print(results['14']['spiral'])
    scenarios = scenarios_streaming
    headers = ["Database"] + [h.replace("-","") for h in headers_streaming]
    table = [[1<<scenario] + ["{:.2f}".format(float(results[str(scenario)][system]["tput"])) for system in systems_streaming] for scenario in scenarios]
    print(tabulate(table, headers, tablefmt="plain"))

def postprocess_results_table(results, scenarios, systems):
    streaming = type(scenarios[0]) != tuple
    for scenario in scenarios:
        for system in systems:
            dat = results[str(scenario)][system]
            resp_sz = dat["resp_sz"]
            item_sz = dat["item_sz"] if streaming else scenario[1]

            results[str(scenario)][system]["rate"] = item_sz / resp_sz
            results[str(scenario)][system]["param_sz"] = get_pp_size(system, dat)
            if not streaming:
                total_us = dat["total_us"]
                results[str(scenario)][system]["tput"] = int(round(get_total_size(scenario) / total_us)) if total_us > 0 else "-"
                results[str(scenario)][system]["cost"] = get_cost(total_us, resp_sz)

def display_streaming(figurename, results):
    systems = ["spiral", "spiral-pack", "spiralstream", "spiralstream-pack", "fastpir", "onionpir"]
    headers = ["Spiral", "SpiralPack",  "SpiralStream", "SpiralStreamPack",  "FastPir", "OnionPir"]
    scenarios = scenarios_streaming_disp
    postprocess_results_table(results, scenarios, systems)
    headers = ["Database", "Metric"] + headers
    table = []
    for scenario in scenarios:
        for key in ["param_sz", "query_sz", "rate", "tput"]:
            for system in systems:
                assert check_corr(results[str(scenario)][system])
            row = [get_row_start(figurename, scenario, key), get_key_tex(key)]
            raw_data = [results[str(scenario)][system][key] for system in systems]
            data = emph_data(key, raw_data)
            table.append(row + data)
    print(finish_output_str(figurename, table, headers))

def gen_table(figurename, trials=1):
    results = {}
    scenarios = scenarios_table
    for scenario in scenarios:
        results[str(scenario)] = {}
        for system in systems_table:
            result = run_system_tr(system, scenario[0], scenario[1], False, trials=trials)
            results[str(scenario)][system] = result
    save_file(figurename, results)
    return results

def get_cost(total_us, resp_sz):
    return usd_per_us*total_us + resp_sz*usd_per_byte

def get_cost_r(res):
    return get_cost(res["total_us"], res["resp_sz"])

def display_table(figurename, results):
    scenarios = scenarios_table
    systems_excl_nopriv = list(filter(lambda x: x != "nopriv", systems_table))
    postprocess_results_table(results, scenarios, systems_excl_nopriv)
    headers = ["Database", "Metric"] + headers_table
    blank_idx = 4
    table = []
    table.append(insert_blank(get_pp_row(results, scenarios, systems_excl_nopriv), blank_idx))
    for scenario in scenarios:
        for key in ["query_sz", "resp_sz", "total_us", "rate", "tput", "cost"]:
            row = [get_row_start(figurename, scenario, key), get_key_tex(key)]
            assert all([check_corr(results[str(scenario)][system]) for system in systems_excl_nopriv])
            raw_data = [results[str(scenario)][system][key] for system in systems_excl_nopriv]
            data = emph_data(key, raw_data)
            table.append(insert_blank(row + data, blank_idx))
    print(finish_output_str(figurename, table, headers))

def gen_ubench(figurename, trials=1):
    results = {}
    scenarios = scenarios_ubench
    system = "spiral"
    for scenario in scenarios:
        results[str(scenario)] = {}
        result = run_system_tr(system, scenario[0], scenario[1], False, trials=trials)
        results[str(scenario)][system] = result
    save_file(figurename, results)
    return results

def display_ubench(figurename, results):
    scenarios = scenarios_ubench
    system = "spiral"
    headers = ["Database"] + headers_ubench
    table = [[1<<scenario[0]] + [int(results[str(scenario)][system][k]) for k in keys_ubench] for scenario in scenarios]
    print(tabulate(table, headers, tablefmt="plain"))
    assert all([check_corr(results[str(scenario)][system]) for scenario in scenarios])

def gen_ablation(figurename, trials=1):
    results = {}
    scenarios = scenarios_ablation
    system = "spiral"
    flags = [None, "-DNO_CRT"]#, "-mno-avx512f", "-mno-avx2"]
    for scenario in scenarios:
        results[str(scenario)] = {}
        for flag in flags:
            result = run_system_tr(system, scenario[0], scenario[1], False, addtl_flags=flag, trials=trials)
            results[str(scenario)][str(flag)] = result
    save_file(figurename, results)
    return results

def display_ablation(figurename, results):
    scenarios = scenarios_ablation
    scenario_names = ["{2^{20}\\times 256 \\text{B}}", "{2^{18}\\times 20\KB}"]
    system = "spiral"
    flags = [None, "-DNO_CRT"]#, "-mno-avx512f", "-mno-avx2"]
    headers = ["Database", "Baseline", "No CRT", "No AVX-512", "No AVX2"]
    table = [[scenario_names[i]] + [int(results[str(scenario)][str(flag)]["total_us"]) for flag in flags] for i,scenario in enumerate(scenarios)]
    print(tabulate(table, headers, tablefmt="plain"))

def scale_by(r, rounds, is_spiral=False):
    o = {
        "query_sz": r["query_sz"],
        "resp_sz": r["resp_sz"] * rounds,
        "total_us": r["total_us"] * rounds 
    }
    if is_spiral:
        o["total_us"] = (r["fdim_us"] + r["fold_us"])*rounds + r["conv_us"] + r["exp_us"]
    if 'params' in r:
        o['params'] = r['params']
    return o

def get_voice(rounds, trials):
    results = {}
    r = run_system_tr("nopriv", 20, 96, False, trials=trials)
    results["nopriv"] = scale_by(r, rounds)
    r = run_system_tr("fastpir", 20, 96, False, trials=trials)
    results["fastpir"] = scale_by(r, rounds)
    r = run_system_tr("spiralstream", 14, 6144, False, cmd_extras="--ignore-expansion --max-query-size 33000000", trials=trials)
    results["spiralstream"] = scale_by(r, rounds, is_spiral=True)
    r = run_system_tr("spiralstream-pack", 14, 6144, False, cmd_extras="--ignore-expansion --max-query-size 33000000", trials=trials)
    results["spiralstream-pack"] = scale_by(r, rounds, is_spiral=True)
    return results

systems_application_movie = ["nopriv", "spiralstream", "spiralstream-pack"]#["nopriv", "onionpir", "spiralstream"]
systems_application_wiki = ["nopriv", "spiral", "spiral-pack", "spiralstream", "spiralstream-pack"]#["nopriv", "onionpir", "spiral"]

def gen_application(figurename, trials=1):
    results = {"movie":{}, "wiki":{}, "voice":{}}
    # Movie scenario
    for i, system in enumerate(systems_application_movie):
        scenario = (14, 2000000000)
        results["movie"][system] = run_system_tr(system, scenario[0], scenario[1], False, trials=trials, cmd_extras="--max-query-size 33000000")
    # Wikipedia scenario
    for i, system in enumerate(systems_application_wiki):
        scenario = (20, 30000)
        results["wiki"][system] = run_system_tr(system, scenario[0], scenario[1], False, trials=trials)
    # Voice call scenario
    results["voice"] = get_voice(625, trials)
    save_file(figurename, results)
    return results

def get_e2e_time(res, cores=16):
    psz = res["param_sz"] if "param_sz" in res else 0
    k = 1 # number of queries of a given visitor
    return (res["total_us"] / (cores*1000000)) + (res["resp_sz"] / download_speed_B_s) + (res["query_sz"] / upload_speed_B_s)# + (psz / upload_speed_B_s)/k

def display_application(figurename, results):
    print(results)
    table = []
    metrics = ["Upload", "Download", "ServerComputation", "ParamSz"]
    metric_keys = ["query_sz", "resp_sz", "total_us", "param_sz"]
    headers = ["Scenario", "System"] + metrics + ["Cost"]
    # Movie
    for i, system in enumerate(systems_application_movie):
        assert check_corr(results["movie"][system])
        table.append(["2^14 x 2 GB movie", system] + [results["movie"][system][metric] if metric in results["movie"][system] else "-" for metric in metric_keys] + [get_cost_r(results["movie"][system])])
    # Wikipedia
    for i, system in enumerate(systems_application_wiki):
        assert check_corr(results["wiki"][system])
        table.append(["2^20 x 30 KB encyclopedia", system] + [results["wiki"][system][metric] if metric in results["wiki"][system] else "-" for metric in metric_keys] + [get_e2e_time(results["wiki"][system])])
    # Voice
    for i, system in enumerate(["nopriv", "fastpir", "spiralstream", "spiralstream-pack"]):
        assert check_corr(results["voice"][system])
        table.append(["5 min voice call, 2^20 users", system] + [results["voice"][system][metric] if metric in results["voice"][system] else "-" for metric in metric_keys] + [get_cost_r(results["voice"][system])])
    print(tabulate(table, headers, tablefmt="plain", disable_numparse=True))

def get_sys_pc(system, scenario):
    if system == "best":
        return best_sys_for_scenario_pc[scenario]
    else:
        return system

def get_systems_pc(scenario):
    return [get_sys_pc(system, scenario) for system in systems_packingcomp]

def gen_packingcomp(figurename, trials=1):
    results = {}
    scenarios = scenarios_packingcomp
    for scenario in scenarios:
        results[str(scenario)] = {}
        for system in get_systems_pc(scenario):
            result = run_system_tr(system, scenario[0], scenario[1], False, cmd_extras="--max-query-size 33000000", trials=trials)
            results[str(scenario)][system] = result
    save_file(figurename, results)
    return results

def postprocess_results_packingcomp(results):
    scenarios = scenarios_packingcomp
    for scenario in scenarios:
        num_items_log2 = scenario[0]
        item_size = scenario[1]
        for system in get_systems_pc(scenario):
            dat = results[str(scenario)][system]
            resp_sz = dat["resp_sz"]
            total_us = dat["total_us"]

            results[str(scenario)][system]["rate"] = item_size / resp_sz
            results[str(scenario)][system]["tput"] = int(round(get_total_size(scenario) / total_us)) if total_us > 0 else "-"
            results[str(scenario)][system]["param_sz"] = get_pp_size(system, dat)

def display_packingcomp(figurename, results):
    postprocess_results_packingcomp(results)
    scenarios = scenarios_packingcomp
    headers = ["Database", "Metric"] + headers_packingcomp
    table = []
    for scenario in scenarios:
        for key in ["param_sz", "query_sz", "resp_sz", "total_us", "rate", "tput"]:
            row = [get_row_start(figurename, scenario, key), get_key_tex(key)]
            assert all([check_corr(results[str(scenario)][system]) for system in get_systems_pc(scenario)])
            raw_data = [result[str(scenario)][system][key] for system in get_systems_pc(scenario)]
            data = emph_data(key, raw_data)
            table.append(row + data)
    print(finish_output_str(figurename, table, headers))

def gen_limits(figurename, trials=1):
    results = {}
    scenarios = scenarios_limits
    for scenario in scenarios:
        results[str(scenario)] = {}
        for system in systems_limits:
            results[str(scenario)][system] = {}
            result_rate = run_system_tr(system, scenario[0], scenario[1], False, cmd_extras="--optimize-for rate --max-query-size 33000000", trials=trials)
            result_tput = run_system_tr(system, scenario[0], scenario[1], False, cmd_extras="--optimize-for tput --max-query-size 33000000", trials=trials)
            results[str(scenario)][system]["optim-rate"] = result_rate
            results[str(scenario)][system]["optim-tput"] = result_tput
    save_file(figurename, results)
    return results

def postprocess_results_limits(results):
    scenarios = scenarios_limits
    for scenario in scenarios:
        num_items_log2 = scenario[0]
        item_size = scenario[1]
        for system in systems_limits:
            for criteria in ["rate", "tput"]:
                okey = "optim-"+criteria
                dat = results[str(scenario)][system][okey]
                resp_sz = dat["resp_sz"]
                total_us = dat["total_us"]

                results[str(scenario)][system][okey]["rate"] = item_size / resp_sz
                results[str(scenario)][system][okey]["tput"] = int(round(get_total_size(scenario) / total_us)) if total_us > 0 else "-"
                results[str(scenario)][system][okey]["param_sz"] = get_pp_size(system, dat)

def display_limits(figurename, results):
    postprocess_results_limits(results)
    scenarios = scenarios_limits
    headers = ["Database", "System"] + [get_key_tex(k) for k in keys_limits_disp]
    table = []
    for scenario in scenarios:
        for criteria in ["rate", "tput"]:
            best_system = max(systems_limits, key=lambda s: result[str(scenario)][s]["optim-"+criteria][criteria])
            cur_res = result[str(scenario)][best_system]["optim-"+criteria]
            assert check_corr(cur_res)

            row = [get_db_tex(scenario), "\\"+best_system.replace("-","")]
            data = [metric_str_output(k, cur_res[k]) for k in keys_limits_disp]
            idx_target = keys_limits_disp.index(criteria)
            data[idx_target] = emph_str(data[idx_target])
            table.append(row + data)
    print(finish_output_str(figurename, table, headers))

def gen_maxtotalquery(figurename, trials=1):
    results = {}
    scenarios = scenarios_maxtotalquery
    for scenario in scenarios:
        results[str(scenario)] = {}
        for system in systems_maxtotalquery:
            results[str(scenario)][system] = {}
            result_rate = None
            try:
                result_rate = run_system_tr(system, scenario[0], scenario[1], False, cmd_extras="--optimize-for rate --max-"+str(scenario[2])+"-size "+str(scenario[3]), trials=trials)
            except KeyboardInterrupt:
                sys.exit(0)
            except:
                pass
            result_tput = None
            try:
                result_tput = run_system_tr(system, scenario[0], scenario[1], False, cmd_extras="--optimize-for tput --max-"+str(scenario[2])+"-size "+str(scenario[3]), trials=trials)
            except KeyboardInterrupt:
                sys.exit(0)
            except:
                pass
            results[str(scenario)][system]["optim-rate"] = result_rate
            results[str(scenario)][system]["optim-tput"] = result_tput
    save_file(figurename, results)
    return results

def display_maxtotalquery_breakdown(figurename, results):
    for criteria in ["rate", "tput"]:
        print(criteria.capitalize()+":")
        headers = [criteria.capitalize(), "totalquery", "param", "query"] #["Max"] + [x.replace("-","") for x in max_predicates]
        table = []
        for max_val in sorted(list(set(map(lambda s: s[3], scenarios_maxtotalquery)))):
            row = []
            scenarios = list(filter(lambda s: s[3] == max_val, scenarios_maxtotalquery))
            for max_predicate in ["total-query"]: # assumes only 1 setting of nu_1/nu_2
                scenario = list(filter(lambda s: s[2] == max_predicate, scenarios))[0]
                res = lambda s: results[str(scenario)][s]["optim-"+criteria]
                best_system = max(systems_maxtotalquery, key=lambda s: res(s)[criteria] if res(s) is not None else -1)
                run = res(best_system)
                if run is None:
                    continue
                val = "{:.4f}".format(float(run[criteria]))
                if best_system == "spiralstream":
                    val = val+"(*)"                
                # dat = results[str(scenario)][best_system]["optim-"+criteria]
                # print(dat["param_sz"]+dat["query_sz"] < scenario[2], dat["param_sz"]+dat["query_sz"], scenario[2])
                row.append(run[criteria]) # tput/rate
                row.append(run["param_sz"] + run["query_sz"]) # total query size
                row.append(run["param_sz"]) # public params size
                row.append(run["query_sz"]) # online query size
            if len(row) > 0:
                table.append(row)
        print(tabulate(table, headers, tablefmt="plain", disable_numparse=True))

def display_maxtotalquery(figurename, results):
    max_vals = sorted(list(set(map(lambda s: s[3], scenarios_maxtotalquery))))
    runs = {}
    for criteria in ["rate", "tput"]:
        runs[criteria] = {}
        for predicate in ['query', 'total-query']:
            runs[criteria][predicate] = {"all": {}}
            for system in systems_maxtotalquery:
                runs[criteria][predicate][system] = {}

            table = []
            for max_val in max_vals:
                row = []
                scenario = list(filter(lambda s: s[2] == predicate and s[3] == max_val, scenarios_maxtotalquery))[0]
                res = lambda s: results[str(scenario)][s]["optim-"+criteria]
                best_system = max(systems_maxtotalquery, key=lambda s: res(s)[criteria] if res(s) is not None else -1)
                run = res(best_system)
                if run is None:
                    continue
                assert check_corr(run)
                val = "{:.4f}".format(float(run[criteria]))
                runs[criteria][predicate][best_system][max_val] = run
                runs[criteria][predicate]["all"][max_val] = run

    for system in ["all"] + systems_maxtotalquery[1:]:
        for predicate in ["query", "total-query"]:
            for criteria in ["rate", "tput"]:
                crit = criteria.capitalize()
                predicate_name = "OnlineQuery" if predicate == 'query' else "TotalQuery"
                pred = predicate_name.lower()

                print("% "+system+", "+predicate+", "+crit+":")
                print("\pgfplotstableread{")
                table = []
                headers = [predicate_name, crit]
                for max_val in max_vals:
                    objs = runs[criteria][predicate][system]
                    if max_val not in objs:
                        continue
                    run = objs[max_val]
                    table.append([max_val, run[criteria]])
                print(tabulate(table, headers, tablefmt="plain", disable_numparse=True))
                print("}\maxupload"+system.replace("-","")+pred+criteria)
    print()
    for system in ["onionpir", "fastpir"]:
        scenario = scenarios_maxtotalquery[0]
        scenario = (scenario[0], scenario[1])
        result_table = load_file("table")[str(scenario)]
        res = result_table[system]
        for predicate in ["query", "total-query"]:
            for criteria in ["rate", "tput"]:
                max_val = res["query_sz"]
                if predicate == "total-query":
                    max_val += get_pp_size(system, None)
                obj_val = scenario[1] / res["resp_sz"]
                if criteria == "tput":
                    obj_val = int(round(get_total_size(scenario) / res["total_us"]))
                
                crit = criteria.capitalize()
                predicate_name = "OnlineQuery" if predicate == 'query' else "TotalQuery"
                pred = predicate_name.lower()

                # print("% "+system+", "+predicate+", "+crit+":")
                print("\pgfplotstableread{")
                table = [[max_val, obj_val]]
                headers = [predicate_name, crit]
                print(tabulate(table, headers, tablefmt="plain", disable_numparse=True))
                print("}\maxupload"+system.replace("-","")+pred+criteria)

def parse_all_args():
    parser = argparse.ArgumentParser(description='Generate figures.')
    parser.add_argument('figures', nargs='+', default=[],
                        help='figures to generate (asympcomp, asympcomplarge, streaming, table, ubench)')
    parser.add_argument('--load', help='just load the results from pkl files', action='store_true')
    parser.add_argument('--spiral-only', help='only run spiral and variants', action='store_true')
    parser.add_argument('--path', help='path to load from', type=str)
    parser.add_argument('--graph', help='prefer graphs instead of tables', action='store_true')
    parser.add_argument('--trials', metavar='trials', type=int, nargs='?', const=1, default=1,
                        help='number of trials to run over')
    return parser.parse_args()

all_figures = ["asympcomp", "asympcomplarge", "streaming", "table", "ubench", "ablation", "application", "packingcomp", "limits", "maxtotalquery"]
fns_table = {
    "asympcomp": (gen_asympcomp, display_asympcomp),
    "asympcomplarge": (gen_asympcomp, display_asympcomp),
    "streaming": (gen_streaming, display_streaming),
    "table": (gen_table, display_table),
    "ubench": (gen_ubench, display_ubench),
    "ablation": (gen_ablation, display_ablation),
    "application": (gen_application, display_application),
    "packingcomp": (gen_packingcomp, display_packingcomp),
    "limits": (gen_limits, display_limits),
    "maxtotalquery": (gen_maxtotalquery, display_maxtotalquery)
}

if __name__ == "__main__":
    args = parse_all_args()

    if args.spiral_only:
        systems_packingcomp = systems_packingcomp[1:]
        headers_packingcomp = headers_packingcomp[1:]

    figures = args.figures
    if "all" in figures:
        figures = all_figures
    if args.graph:
        fns_table["streaming"] = (gen_streaming, display_streaming_graph)

    for figure in figures:
        start = timer()
        if len(figures) > 1:
            print(f"Figure \"{figure}\":")

        gen_fn = fns_table[figure][0]
        display_fn = fns_table[figure][1]
        result = None
        if args.load:
            result = load_file(figure, args.path)
        else:
            result = gen_fn(figure, trials=args.trials)
        display_fn(figure, result)

        if len(figures) > 1:
            print()
        end = timer()
        print(f"Took {end - start} s")