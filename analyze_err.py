import numpy as np
import math
import tabulate
import sys

def modulus_cutoff(errs, bins, p):
    error_rates = []
    for q in bins:
        count = 0
        for e in errs:
            if abs(e) * (p/q) > 0.5:
                count += 1
        error_rates.append(count / len(errs))
    return error_rates

def extend_subg(error_rate, modulus, p):
    logq = math.log(modulus,2)
    logp = math.log(p,2)
    logpi = math.log(math.pi,2)
    s_e_log2 = 2*(logq - (logp+1)) + logpi - math.log(math.log(2) - math.log(error_rate),2)
    # s_e_log2 = 2*logp - math.log((-4) * (modulus**2) * math.pi * math.log(error_rate/2), 2)
    # s_e_log2 = - math.log(4 * math.pi) - math.log(-math.log((1-error_rate)/2),2) - 2*math.log(p,2) + 2*math.log(modulus, 2)
    return s_e_log2

p = int(sys.argv[1])
data = open(sys.argv[2], 'r').read()
errs = [int(i) for i in data.strip().split(" ")]
print(len(errs))

bins = [2**i for i in np.arange(40,60,0.1)]
error_rates = modulus_cutoff(errs, bins, p)

# remove values close to zero
min_observations = 5
num_zeros = 0
for i in reversed(error_rates):
    if i > min_observations/len(errs):
        print(i)
        break
    num_zeros += 1
corr_error_rates = error_rates[:-num_zeros]
bins_for_errs = bins[:len(corr_error_rates)]

print(tabulate.tabulate(zip([math.log(i,2) for i in bins_for_errs], corr_error_rates), ["modulus", "err_rate"], tablefmt="plain"))

last_err = corr_error_rates[-1]
last_mod = bins_for_errs[-1]
print(last_err, last_mod)
extended_s_e_log2 = extend_subg(last_err, last_mod, p)

print("Extended subg. width:", extended_s_e_log2)