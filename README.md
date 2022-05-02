> **⚠️ IMPORTANT ⚠️**
> 
> A more complete and documented Rust implementation of Spiral is available at https://github.com/menonsamir/spiral-rs. Future work is being done on this new implementation; we preserve this old repository for reference purposes only, since our original benchmarks were collected using this code. We will not be updating this code.

# Spiral: Fast, High-Rate Single-Server PIR via FHE Composition

This is an implementation of our paper "Spiral: Fast, High-Rate Single-Server PIR via FHE Composition", available [here](https://eprint.iacr.org/2022/368.pdf).

> **WARNING**: This is research-quality code; it has not been checked for side-channel leakage or basic logical or memory safety issues. Do not use this in production.

## Building

To build Spiral, you will need:
- Intel HEXL 1.2.1, which should be installed using vcpkg
- Python 3
- A modern C++ compiler (Clang 12+ is ideal)

From a fresh install of Ubuntu 20.04, the steps are:

```
sudo apt-get update
sudo apt-get install -y build-essential clang-12 git-lfs
git lfs install
cd ~
git clone https://github.com/Microsoft/vcpkg.git
./vcpkg/bootstrap-vcpkg.sh -disableMetrics
./vcpkg/vcpkg install hexl
git clone https://github.com/menonsamir/spiral.git
cd spiral
python3 select_params.py 20 256
```

The `select_params.py` performs the actual Spiral build using CMake automatically.

To replicate the key table from the paper, just run:

```
python3 run_all.py packingcomp --spiral-only
```

All other figures from the paper can also be generated using this script. We used AWS EC2 `c5n.2xlarge` instances to produce our results.

To perform a manual build, you can run something like:

```
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=~/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build -j4 -- PARAMSET=PARAMS_DYNAMIC \
  TEXP=8 TEXPRIGHT=56 TCONV=4 TGSW=8 QPBITS=20 PVALUE=256 \
  QNUMFIRST=1 QNUMREST=0 OUTN=2 
./spiral 8 7 1234
```

The build steps are roughly similar on other platforms. If vcpkg is not installed at `~/vcpkg`, you will need to specify its location in the `--vcpkg-root` argument. Building on M1-based macOS devices is quite tricky - support for this is a TODO for now.

## Running

To run Spiral on a database of `2^20` (~1 million) items of 256 bytes each, run:
```
python3 select_params.py 20 256 --show-output
```

The script will perform automatic parameter selection, compile `./spiral`, and then perform a full end-to-end PIR test, outputting a JSON string summarizing the results.

To run the Spiral variants:

| Variant          | Options                  |
| ---------------- | ------------------------ |
| Spiral           | none                     |
| SpiralStream     | `--direct-upload`        |
| SpiralPack       | `--pack`                 |
| SpiralStreamPack | `--direct-upload --pack` |

By default, we use a implicit database representation, to allow measurement of large databases without concern for memory usage. To use an explicit database representation, pass `--explicit-db`. The details of other options are availible through `python3 select_params.py --help`.

## Parameters

This repository contains the parameters collected using `generate_all_schemes.py` cached as `.pkl` files. These contain all parameter sets in the search space that satisfy our correctness threshold of `2^-40`. If you'd like to change the search space or correctness threshold, you will need to regenerate these. To regenerate all parameters, simply run:

```
python3 generate_all_schemes.py && \
python3 generate_all_schemes.py --stream && \
python3 generate_all_schemes.py --high-rate && \
python3 generate_all_schemes.py --stream --high-rate
```

The cost model used to estimate the running time of a given parameter set has been constructed manually through regression and lookup tables on an AWS EC2 `c5n.2xlarge` instance. The accuracy of this cost model could be significantly degraded on systems without AVX2 support.
