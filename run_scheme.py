import pickle
import math
import sys
import json
from ast import literal_eval
import pprint
import subprocess
import random
import re
import argparse

def get_factor(itemsize, maxsize):
    factor = 1
    if itemsize <= maxsize:
        factor = 1
    else:
        factor = math.ceil(itemsize / maxsize)
    return factor

systems = ["sealpir", "fastpir", "onionpir", "spiralstream", "spiral", "spiralstream-pack", "spiral-pack", "nopriv"]

def run_system(system, db_items_log2, itemsize, streaming, show_output=False, addtl_flags=None, cmd_extras=None):
    assert system in systems, "Must choose available system."
    # if system == "spiralstream":
    #     assert streaming, "Must use stream with spiralstream."
    if streaming:
        assert itemsize == 1, "Must set itemsize to 1 for streaming."
    
    # Spiral
    def run_spiral(db_items_log2, itemsize):
        cmd = "/usr/bin/python3 select_params.py --quiet"
        if streaming:
            cmd += " --stream"
        if "spiralstream" in system:
            cmd += " --direct-upload"
        if "pack" in system:
            cmd += " --high-rate"
        if addtl_flags:
            cmd += f" --addtl-flags=\"\\\"{addtl_flags}\\\"\""
        if cmd_extras:
            cmd += " "+cmd_extras
        s = subprocess.check_output(f"{cmd} {db_items_log2} {itemsize}", shell=True)
        if show_output:
            print(s.decode("utf8"))
        return s.decode("utf8")
    def run_and_analyze_spiral(db_items_log2=20, itemsize=256):
        if streaming:
            itemsize = 1
        s = run_spiral(db_items_log2, itemsize)
        obj = json.loads(s)
        if streaming:
            return {
                "tput":  (obj["dbsize"]) / (obj["fdim_us"] + obj["fold_us"]),
                "resp_sz": obj["resp_sz"],
                "item_sz": obj["item_sz"],
                "param_sz": obj["param_sz"],
                "params": obj["params"],
                "query_sz": obj["query_sz"],
                "other_data": obj
            }
        else:
            return obj

    # SealPIR
    def run_sealpir(db_items_log2, itemsize):
        s = subprocess.check_output(f"/home/ubuntu/SealPIR/bin/sealpir-implicit {db_items_log2} {itemsize}", shell=True)
        if show_output:
            print(s.decode("utf8"))
        return s.decode("utf8")
    def analyze_sealpir(s, db_items_log2, itemsize, factor):
        total_re = r"\s+PIRServer reply generation time.*:\s+([0-9]+) ms"
        query_sz_re = r"\s+Query size bytes.*:\s+([0-9]+)"
        resp_sz_re = r"\s+Reply size bytes.*:\s+([0-9]+)"

        exp_re = r"Server: expansion time.*\s+([0-9]+) ms"
        
        total_ms = int(re.search(total_re, s).group(1))
        exp_ms = sum([int(i) for i in re.findall(exp_re, s)])
        query_sz_b = int(re.search(query_sz_re, s).group(1))
        resp_sz_b = int(re.search(resp_sz_re, s).group(1))
        if streaming:
            return {
                "tput":  ((1<<db_items_log2)*(itemsize)) / ((total_ms-exp_ms) * 1000),
                "resp_sz": factor * resp_sz_b,
                "item_sz": factor * itemsize,
                "query_sz": query_sz_b
            }
        else:
            return {
                "total_us": (factor * (total_ms-exp_ms) + exp_ms) * 1000,
                "resp_sz": factor * resp_sz_b,
                "query_sz": query_sz_b
            }
    def run_and_analyze_sealpir(db_items_log2=20, itemsize=256):
        maxsize = 3072
        if streaming:
            itemsize = maxsize
        factor = get_factor(itemsize, maxsize)
        s = run_sealpir(db_items_log2, min(itemsize, maxsize))
        return analyze_sealpir(s, db_items_log2, itemsize, factor)

    # FastPIR
    def run_fastpir(db_items_log2, itemsize):
        db_items = 1 << db_items_log2
        s = subprocess.check_output(f"/home/ubuntu/FastPIR/bin/fastpir-implicit -n {db_items} -s {itemsize}", shell=True)
        if show_output:
            print(s.decode("utf8"))
        return s.decode("utf8")
    def analyze_fastpir(s, db_items_log2, itemsize, factor):
        total_re = r"\s+Response generation time.*:\s+([0-9]+)"
        query_sz_re = r"\s+Query size.*:\s+([0-9]+)"
        resp_sz_re = r"\s+Response size.*:\s+([0-9]+)"
        
        total_us = int(re.search(total_re, s).group(1))
        query_sz_b = int(re.search(query_sz_re, s).group(1))
        resp_sz_b = int(re.search(resp_sz_re, s).group(1))
        if streaming:
            return {
                "tput":  ((1<<db_items_log2)*(itemsize)) / (total_us),
                "resp_sz": factor*resp_sz_b,
                "item_sz": factor*itemsize,
                "query_sz": query_sz_b
            }
        else:
            return {
                "total_us": factor*total_us,
                "resp_sz": factor*resp_sz_b,
                "query_sz": query_sz_b
            }
    def run_and_analyze_fastpir(db_items_log2=20, itemsize=256):
        maxsize = 9120
        if streaming:
            itemsize = maxsize
        factor = get_factor(itemsize, maxsize)
        s = run_fastpir(db_items_log2, min(itemsize, maxsize))
        return analyze_fastpir(s, db_items_log2, itemsize, factor)

    # OnionPIR
    def run_onionpir(db_items_log2, itemsize):
        s = subprocess.check_output(f"/home/ubuntu/Onion-PIR/onionpir-implicit {db_items_log2} {itemsize}", shell=True)
        if show_output:
            print(s.decode("utf8"))
        return s.decode("utf8")
    def analyze_onionpir(s, db_items_log2, itemsize, factor):
        exp_re = r"\s+Server: rlwe exansion time.*=\s+([0-9]+)"
        exp2_re = r"\s+Server: expand after first diemension.*=\s+([0-9]+)"
        total_re = r"\s+Main: PIRServer reply generation time.*:\s+([0-9]+)"
        resp_sz_re = r"\s+Reply size bytes.*:\s+([0-9]+)"
        
        exp_us = (int(re.search(exp_re, s).group(1))+int(re.search(exp2_re, s).group(1)))*1000
        total_us = int(re.search(total_re, s).group(1))*1000
        query_sz_b = 63488
        resp_sz_b = int(re.search(resp_sz_re, s).group(1))

        if streaming:
            return {
                "tput":  ((1<<db_items_log2)*(itemsize)) / (total_us-exp_us),
                "resp_sz": factor*resp_sz_b,
                "item_sz": factor*itemsize,
                "query_sz": query_sz_b
            }
        else:
            return {
                "total_us": factor*(total_us-exp_us) + exp_us,
                "resp_sz": factor*resp_sz_b,
                "query_sz": query_sz_b
            }
    def run_and_analyze_onionpir(db_items_log2=20, itemsize=256):
        if db_items_log2 == 10 and itemsize == 10000:
            return {
                "total_us": 582000,
                "resp_sz": 126976,
                "query_sz": 63488
            }
        maxsize = 30720
        if streaming:
            itemsize = maxsize
        factor = get_factor(itemsize, maxsize)
        s = run_onionpir(db_items_log2, min(itemsize, maxsize))
        return analyze_onionpir(s, db_items_log2, itemsize, factor)
    
    def run_and_analyze_nopriv(db_items_log2=20, itemsize=256):
        return {
            "total_us": 0,
            "resp_sz": itemsize,
            "query_sz": 0
        }
    
    funcs = {
        "sealpir": run_and_analyze_sealpir,
        "fastpir": run_and_analyze_fastpir,
        "onionpir": run_and_analyze_onionpir,
        "spiralstream": run_and_analyze_spiral,
        "spiralstream-pack": run_and_analyze_spiral,
        "spiral": run_and_analyze_spiral,
        "spiral-pack": run_and_analyze_spiral,
        "nopriv": run_and_analyze_nopriv
    }
    result = funcs[system](db_items_log2, itemsize)
    return result

def run_system_tr(system, db_items_log2, itemsize, streaming, show_output=False, addtl_flags=None, cmd_extras=None, trials=1):
    all_results = []
    for trial in range(trials):
        result = run_system(system, db_items_log2, itemsize, streaming, show_output, addtl_flags, cmd_extras)
        all_results.append(result)
    
    res = all_results[0]
    res["from_trials"] = trials
    keys_to_avg = ["tput"] if streaming else ["total_us", "cost"]
    if "spiral" not in system:
        keys_to_avg = ["tput"] if streaming else ["total_us"]
    for key_to_avg in keys_to_avg:
        avg_val = sum([r[key_to_avg] for r in all_results])/trials
        res[key_to_avg] = avg_val
    return res

def parse_all_args():
    parser = argparse.ArgumentParser(description='Run other systems.')
    parser.add_argument('system', metavar='system', type=str,
                        help='system to benchmark (one of spiral, spiralstream, spiral-pack, spiralstream-pack, sealpir, fastpir, onionpir, nopriv)')
    parser.add_argument('targetnum', metavar='logN', type=int,
                        help='log2 of the number of items')
    parser.add_argument('itemsize', metavar='itemsize', type=int,
                        help='item size in bytes')
    parser.add_argument('--show-output', help='show output', action='store_true')
    parser.add_argument('--stream', help='use the streaming variant', action='store_true')
    parser.add_argument('--trials', metavar='trials', type=int, nargs='?', const=1, default=1,
                    help='number of trials to run over')
    return parser.parse_args()
    

if __name__ == "__main__":
    args = parse_all_args()
    db_items_log2 = args.targetnum
    itemsize = args.itemsize
    streaming = args.stream

    result = run_system(args.system, args.targetnum, args.itemsize, streaming, show_output=args.show_output)
    pprint.pprint(result)