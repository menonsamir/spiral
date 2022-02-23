import math

other_pp_sz = {
    "onionpir": 4600000,
    "fastpir": 1400000,
    "sealpir": 3400000
}

def get_pp_size(system, r=None):
    if 'spiral' in system:
        if "param_sz" in r:
            return r["param_sz"]
        else:
            return r["other_data"]["param_sz"]
        # params = r["params"]
        # if 'pack' in system:
        #     packing_sz = params["n"]*params["t_conv"]*(params["n"]+1)
        #     exp_mats = params["nu_1"]*2*(params["t_exp"]) + int(math.ceil(math.log2(params["nu_1"]*params["t_GSW"])))*2*(params["t_exp_right"])
        #     v_sz = 2 * 2 * params["t_conv"]
        #     total_elems = packing_sz
        #     if 'stream' not in system:
        #         total_elems += exp_mats + v_sz
        #     return int(total_elems * 2048 * 56 / 8)
        # else:
        #     exp_mats = params["nu_1"]*2*(params["t_exp"]) + int(math.ceil(math.log2(params["nu_1"]*params["t_GSW"])))*2*(params["t_exp_right"])
        #     w_sz = 3 * 2 * params["t_conv"]
        #     v_sz = 3 * 2 * params["t_conv"]
        #     total_elems = w_sz + v_sz
        #     if 'stream' not in system:
        #         total_elems += exp_mats
        #     return int(total_elems * 2048 * 56 / 8)
    else:
        return other_pp_sz[system]
