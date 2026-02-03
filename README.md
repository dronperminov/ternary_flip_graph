# ternary_flip_graph

[![arXiv:2511.20317](https://img.shields.io/badge/arXiv-2511.20317-b31b1b.svg)](https://arxiv.org/abs/2511.20317)
[![arXiv:2512.13365](https://img.shields.io/badge/arXiv-2512.21980-b31b1b.svg)](https://arxiv.org/abs/2512.21980)

A comprehensive toolkit for discovering fast matrix multiplication algorithms using flip graph exploration with support for Z₂, Z₃ and ternary coefficient sets.


## Key Features
* Multi-ring support: works with Z₂ (`{0, 1}`), Z₃ (`{0, 1, 2}`) and ternary coefficient sets (`{-1, 0, 1}`);
* No lifting required: ternary coefficients (`{-1, 0, 1}`) produce schemes valid over any ring automatically;
* Parallel exploration: multi-threaded search with configurable parallel runners;
* Multiple optimization strategies: rank minimization, naive complexity minimization and alternative scheme discovery;
* Flexible operations: support for `flip`, `plus`, `split`, `reduce`, `extend`, `merge`, `product`, and `project` operations.


## Tools Overview
### `flip_graph`: single dimension discovery

Discover optimal matrix multiplication schemes for specific dimensions. This tool performs random walks through the flip graph space to find low-rank decompositions.

**Use when**: you need to find the best scheme for a specific matrix multiplication problem (e.g., 3×3×3, 2x4x5).

### `meta_flip_graph`: multi-dimensional meta search
Discover schemes using higher-level operations that transform schemes between dimensions. Supports extension, merging, product, and projection operations.

**Use when**: you want to discover schemes for various matrix dimensions or explore how optimal decompositions can be transformed between different problem sizes.

### `complexity_minimizer`: naive complexity optimization
Given a fixed-rank scheme, find variants with minimal "naive complexity" (number of non-zero coefficients in the decomposition).

**Use when**: you have a rank-optimal scheme and want to minimize its practical implementation cost without using common subexpression elimination (CSE). For optimizing real addition complexity with CSE, please use the companion repository [ternary_addition_reducer](https://github.com/dronperminov/ternary_addition_reducer).

### `find_alternative_schemes`: alternative solution discovery
Find alternative schemes with the same rank but different coefficient patterns, providing implementation choices.

**Use when**: you need multiple valid schemes for the same problem, possibly for hardware-specific optimizations.

### `lift`: Hensel Lifting from Z₂/Z₃ to general ring
Lifts a scheme from the rings Z₂ or Z₃ to higher-order rings (Z₂ → Z₄ → Z₈ → Z₁₆ ... or Z₃ → Z₉ → Z₂₇ → Z₈₁) using Hensel lifting. The method solves the system `Jx = T - UxVxW` and applies rational reconstruction with subsequent validation to ensure correct lifting.

**Use when**: you need to convert a scheme from Z₂/Z₃ to a general ring (ZT, Z, or Q).

## Quick Start

### Installation
```bash
git clone https://github.com/dronperminov/ternary_flip_graph
cd ternary_flip_graph
make -j$(nproc)
```

### Basic Discovery

#### `4x4x4` matrices with ternary coefficients
```bash
./flip_graph -n1 4 -n2 4 -n3 4 --ring ZT --count 64 --threads 16
```

#### single binary scheme initialization from input.txt
```bash
./flip_graph -i input.txt --ring Z2 --count 128
```

#### multiple Z3 scheme initialization from input.txt
```bash
./flip_graph -m -i input.txt --ring Z3 --count 128
```

### Meta Search
```bash
./meta_flip_graph -i input.txt --ring ZT --count 64 --resize-probability 0.1
```

### Complexity Minimization
```bash
./complexity_minimizer -i input_scheme.txt --ring Z2 --count 32 --threads 12
```

## Understanding the Output

### Flip Graph Progress Table
```text
+-----------------------------------------------------------------------------------+
| dimension: 4x4x4            seed: 1766925187                        best rank: 52 |
| threads: 24                 flip iters: 1.0M                        iteration: 64 |
| count: 256                  reset iters: 100.0M                 elapsed: 00:18:34 |
| ring: ZT                    plus diff: 4                     improvements: 6 / 10 |
+===================================================================================+
| runner | scheme rank |   naive    |            |        flips        |    plus    |
|   id   | best | curr | complexity | iterations |  count  | available | iterations |
+--------+------+------+------------+------------+---------+-----------+------------+
|    108 |   52 |   52 |        685 |          0 |   6.20K |         4 |     12.75K |
|    159 |   54 |   52 |        685 |          0 |  16.87K |         4 |     25.07K |
|     93 |   54 |   54 |        650 |       3.0M |   4.71K |         6 |      6.77K |
+--------+------+------+------------+------------+---------+-----------+------------+
- iteration time (last / min / max / mean): 17.97 / 16.30 / 19.07 / 17.25
```

#### Key metrics:

* `best rank`: lowest rank discovered so far;
* `naive complexity`: number of non-zero coefficients in current scheme;
* `iterations`: number of iterations since last rank improvement for this runner;
* `flips count`: number of successful flip operations since last rank improvement;
* `flips available`: number of currently possible flips for the current scheme;
* `plus iterations`: period (in iterations) between calls to the plus operator for this runner.


### Complexity Minimizer Table
```text
+----------------------------------+
| dimension       rank        ring |
|     5x5x5         93          Z2 |
+----------------------------------+
| count: 32 (12 threads)           |
| seed: 1766927866                 |
| best complexity: 843             |
| iteration: 3                     |
| elapsed: 10.11                   |
+==================================+
| runner | naive scheme complexity |
|   id   |    best    |    curr    |
+--------+------------+------------+
| 4      | 843        | 1013       |
| 1      | 843        | 1023       |
| 3      | 843        | 1013       |
+--------+------------+------------+
- iteration time (last / min / max / mean): 1.38 / 0.86 / 1.95 / 1.40
```


## Key Parameters Guide

### Common Parameters
* `-o`: output directory for discovered schemes (default: `schemes`);
* `--format`: output format for saved schemes: `txt` or `json` (default: `json`);
* `--ring`: coefficient ring (`Z2`, `Z3`, or `ZT` for `{-1, 0, 1}`);
* `--count`: number of parallel search processes;
* `--threads`: OpenMP threads per process;
* `--seed`: random seed (`0` for time-based).

### Flip Graph Specific
* Initial scheme specification: two mutually exclusive approaches:
    * Start from naive scheme: use `-n1 -n2 -n3` to begin exploration from the naive matrix multiplication scheme of dimension n₁×n₂×n₃;
    * Start from existing schemes: use `-i` to load schemes from file for further optimization;
    * `-n1`: number of rows in first matrix (A);
    * `-n2`: number of columns in A / rows in second matrix (B);
    * `-n3`: number of columns in second matrix (B);
    * `-i`: path to input file with initial scheme(s):
    * `-m`: boolean flag to read multiple schemes from file (first line contains number of schemes).

**Note**: these options are mutually exclusive. You cannot specify both dimension flags and an input file simultaneously.

* `--flip-iterations`: flip operations between progress reports;
* `--min-plus-iterations`: minimum period for plus operator calls (default: `5K`)
* `--max-plus-iterations`: maximum period for plus operator calls (default: `100K`)
* `--plus-diff`: maximum rank difference for plus operations (default: `4`);
* `--reset-iterations`: total operations before resetting search (default: `100M`);
* `--target-rank`: stop when this rank is found (`0` = find minimum).

**Note about plus iterations**: the plus operator is called every `plus_iterations` iterations, where `plus_iterations` is randomly chosen between `--min-plus-iterations` and `--max-plus-iterations`.

### Meta Flip Graph Additional Parameters
* `--resize-probability`: probability of resize operation (`0.0` to `1.0`, default: `0`). Controls how often the algorithm attempts to change scheme dimensions.
* `--int-width`: integer bit width controlling maximum matrix dimensions (`16|32|64|128`, default: `64`). Smaller values improve performance but limit matrix size.

### Complexity Minimizer
* `-i`: path to input file with initial scheme(s) (required);
* `-m`: boolean flag to read multiple schemes from file (first line contains number of schemes);
* `--max-no-improvements`: stop after `N` iterations without improvement;
* `--plus-probability`: chance to attempt "plus" operation for explore alternative schemes.

### Find Alternative Schemes
* `-i`: path to input file with initial scheme (required);
* `--max-count`: number of alternative schemes to find (default: `10000`).


## Input files format
The tool supports two modes of initialization from files:

* **Single scheme**: the file contains only one scheme;
* **Multiple schemes**: the file starts with the number of schemes, followed by each scheme sequentially.

### Single scheme format
The file contains exactly one scheme without any preceding count.
```txt
2 2 3 11
1 0 0 1 1 0 0 0 0 0 0 1 -1 1 0 0 0 1 0 1 1 0 1 0 0 0 1 -1 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1
1 0 0 0 1 0 0 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 -1 0 -1 1 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 1
1 0 0 1 0 0 0 0 1 -1 0 0 -1 1 0 0 0 0 1 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1
```

### Multiple schemes format
The file begins with an integer specifying the number of schemes, followed by each scheme.
```txt
3
2 2 3 11
1 0 0 1 1 0 0 0 0 0 0 1 -1 1 0 0 0 1 0 1 1 0 1 0 0 0 1 -1 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1
1 0 0 0 1 0 0 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 -1 0 -1 1 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 1
1 0 0 1 0 0 0 0 1 -1 0 0 -1 1 0 0 0 0 1 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1
2 2 3 11
0 0 0 1 0 0 0 1 0 0 1 0 0 0 1 1 0 1 0 0 0 1 0 -1 1 0 0 0 1 0 0 0 1 0 0 1 1 0 -1 0 1 1 0 0
0 0 0 0 0 1 0 -1 0 0 1 0 0 0 1 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 1 0 0 0 1 0 0 -1 0 0 0 1 0 1 0 0 -1 -1 0 0 0 0 0 0 0 1 0 0
0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 1 0 1 0 -1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 1 0 -1 0 0 0
2 2 3 11
0 0 1 0 0 0 1 1 0 1 1 0 0 1 0 0 1 1 0 0 1 0 -1 0 1 0 0 0 0 1 0 -1 0 0 1 -1 0 1 0 0 0 0 0 1
0 1 0 0 -1 0 0 0 0 0 1 0 0 0 1 0 1 0 0 0 1 0 0 -1 0 0 1 0 0 0 0 1 1 0 0 0 1 0 -1 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 1 -1 0 0 1 0 0 1 0 0
0 0 1 1 0 0 0 0 0 1 0 -1 0 0 1 0 0 1 0 0 0 0 -1 -1 1 0 -1 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 -1 0 1 0 0 0 0 -1 0 0 0 0 0 0 1 0 0 0 0
```

For correct run flip graph with this file you need use `-m` flag:
```bash
./flip_graph -i input.txt -m --count 128
```

## Understanding the Approach

### Why Ternary Coefficients?
Traditional approaches using Z₂ or Z₃ require "lifting" - converting schemes to work over arbitrary fields. Ternary coefficient set `{-1, 0, 1}` produces schemes that are automatically valid over any ring, eliminating the lifting step entirely.

### What are Flip Graphs?
Flip graphs represent the space of possible matrix multiplication schemes, where nodes are schemes and edges are transformations ("flips") between them. Introduced by Kauers and Moosbauer in ["Flip Graphs for Matrix Multiplication"](https://arxiv.org/abs/2212.01175), this graph-based approach models matrix multiplication algorithms as tensor decompositions and uses flips to swap components while preserving correctness. By exploring this graph via random walks, optimal decompositions can be discovered.

### What is "Naive Complexity"?
The number of non-zero coefficients in a scheme. Lower naive complexity means fewer operations in practical implementations, even if the theoretical rank is the same. Do not confuse with CSE-optimized complexity - this metric counts all non-zero coefficients without common subexpression elimination. For optimizing real addition complexity with CSE, please see the companion repository [ternary_addition_reducer](https://github.com/dronperminov/ternary_addition_reducer).

### How Does Naive Complexity Minimizer Work?
The complexity_minimizer tool performs a specialized random walk on the flip graph where all schemes have the same fixed rank, but differ in their naive complexity. It explores the neighborhood of the input scheme by applying flips and plus operations, always accepting moves that reduce the number of non-zero coefficients. This process continues until no improvements are found for a specified number of iterations, yielding a scheme with minimal naive complexity for that particular rank.

### How Does Alternative Scheme Discovery Work?
The find_alternative_schemes tool explores the space of schemes with the same rank as the input scheme. Using random flips, it generates structurally different schemes that maintain the same theoretical rank. The tool identifies "alternative" schemes by comparing sorted rows of coefficients - two schemes are considered different if their sorted coefficient patterns don't match exactly. This allows discovery of multiple valid decompositions with the same efficiency, providing options for hardware-specific optimizations or further refinement.


## Citation
If you use this software in your research, please cite:

```bibtex
@article{perminov2025fast,
    title={Fast Matrix Multiplication via Ternary Meta Flip Graphs},
    author={Perminov, Andrew I},
    journal={arXiv preprint arXiv:2511.20317},
    url={https://arxiv.org/abs/2511.20317},
    year={2025}
}
```

```bibtex
@article{perminov202558,
    title={A 58-Addition, Rank-23 Scheme for General 3x3 Matrix Multiplication},
    author={Perminov, Andrew I},
    journal={arXiv preprint arXiv:2512.21980},
    url={https://arxiv.org/abs/2512.21980},
    year={2025}
}
```