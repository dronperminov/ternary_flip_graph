# ternary_flip_graph

[![arXiv:2511.20317](https://img.shields.io/badge/arXiv-2603.02398-b31b1b.svg)](https://arxiv.org/abs/2603.02398)
[![arXiv:2511.20317](https://img.shields.io/badge/arXiv-2511.20317-b31b1b.svg)](https://arxiv.org/abs/2511.20317)
[![arXiv:2512.13365](https://img.shields.io/badge/arXiv-2512.21980-b31b1b.svg)](https://arxiv.org/abs/2512.21980)

A collection of C++ tools for discovering and transforming fast matrix multiplication schemes using the [flip graph](https://arxiv.org/abs/2212.01175)
approach. The repository supports searches over multiple coefficient rings, with a particular focus on integer ternary coefficients `ZT`, which allows
schemes to be valid over the general ring automatically, without any lifting step.

All tools are implemented in pure C++, require only a standard `g++` compiler, have no external dependencies, and support parallel execution via OpenMP.
The implementation handles matrices with up to 128 elements, supporting a wide range of matrix dimensions.

For a comprehensive collection of discovered schemes and research results on fast matrix multiplication schemes, visit the companion repository:
[FastMatrixMultiplication](https://github.com/dronperminov/FastMatrixMultiplication)


## Key features

- **Multiple coefficient rings**: `Z2` binary coefficients `{0, 1}`, `Z3` ternary modular coefficients `{0, 1, 2}`,
  and `ZT` integer ternary coefficients `{-1, 0, 1}`;
- **No lifting required for `ZT`**: schemes over `ZT` are immediately valid over the general ring, avoiding Hensel lifting and rational reconstruction entirely;
- **Flip graph-based search**: local transformations of schemes using `flip`, `plus`, `split`, `reduce`, and related operations, supporting both fixed-dimension
  and variable-dimension searches;
- **Meta-operations**: `extend`, `project`, `merge` and `product` operations allow exploration across different matrix sizes;
- **Parallel and reproducible**: multi-runner architecture combined with OpenMP threading and fully controllable random seeds;
- **Validation and post-processing**: automatic verification of Brent equations, lifting from modular rings when needed, and naive additive complexity
  minimization;
- **Large matrix support**: handles matrices with up to 128 elements, accommodating a wide variety of dimensions.


## Installation

```bash
git clone https://github.com/dronperminov/ternary_flip_graph
cd ternary_flip_graph
make -j$(nproc)
```

The build produces standalone binaries for all tools.


## Tools overview
The repository provides the following command-line tools:

- `flip_graph` — random walk search in the flip graph for fixed dimensions.
- `meta_flip_graph` — flip graph search augmented with dimension-changing meta operations.
- `find_alternative_schemes` — generation of distinct schemes of the same size.
- `lift` — Hensel lifting and rational reconstruction from modular rings (`Z2` / `Z3`).
- `validate_schemes` — verification of Brent equations.
- `optimize_scheme` — optimization of naive additive complexity or potential flips count.

Each tool is described in detail below.


## Core search tools

### flip_graph
`flip_graph` performs random-walk search in the flip graph for schemes of a fixed dimension `(n1, n2, n3)`. The dimension never changes during execution.
Initialization can be done either from a naive matrix multiplication scheme (via dimensions), or from one or more schemes loaded from file. These modes
are mutually exclusive.

#### Main parameters
- `--ring {ZT, Z2, Z3}` — coefficient ring (default: `ZT`);
- `--count INT` — number of parallel runners (default: `8`);
- `--threads INT` — number of OpenMP threads;
- `--format {txt, json}` — output format (default: `txt`);
- `--output-path PATH` — output directory for discovered schemes (default: `schemes`);
- `--seed INT` — random seed, 0 uses time-based seed (default: `0`);
- `--top-count INT` — number of best schemes displayed (default: `10`).

#### Initialization parameters
- `-n1 INT` — rows in matrix A;
- `-n2 INT` — columns in A / rows in B;
- `-n3 INT` — columns in matrix B;
- `--input-path PATH` — path to input file with initial scheme(s);
- `--multiple` — input file contains multiple schemes (read multiple schemes from file, with total count on first line);
- `--no-verify` — skip checking Brent equations for correctness.

#### Random walk parameters
* `--flip-iterations INT` — flips performed before reporting (default: `1M`);
* `--min-plus-iterations INT` — minimal period for `plus` / `split` operations (default: `5K`);
* `--max-plus-iterations INT` — maximal period for `plus` / `split` operations (default: `100K`);
* `--reset-iterations INT` — iterations before reset (default: `10B`);
* `--plus-diff INT` — allowed rank difference for `plus` operation (default: `4`);
* `--sandwiching-probability REAL` — probability of `sandwiching` operation (default: `0`);
* `--reduce-probability REAL` — probability of `reduce` operation (default: `0`).

#### Pool parameters

The pool strategy enables systematic rank reduction by iteratively discovering schemes of decreasing rank. Initially, a pool is populated with schemes of
rank `r`. Each runner randomly selects a scheme from this pool and performs flip operations until a reduction to rank `r-1` is found. Once such a scheme is
discovered, it is added to a new pool, and the runner selects another random scheme from the current pool. When the new pool reaches at least `--pool-size`
schemes, it replaces the current pool, and the search continues for rank `r-2`, then `r-3`, and so on.

- `--use-pool` — enable pool strategy;
- `--pool-size INT` — size of the next-rank pool before switching (default: `1000`).
- `--pool-min-size INT` — minimum number of schemes required in the pool to consider it valid after max iterations (default: `5`).
- `--pool-max-iterations INT` — maximum number of random walk iterations allowed to reach the minimum pool size (default: `1000`).
- `--pool-select-strategy {uniform, flips}` — strategy for selecting schemes from the pool:
  * `uniform` — each scheme is selected with equal probability;
  * `flips` — schemes are selected with weights proportional to the number of potential flip operations.

If `--pool-max-iterations` is reached and the pool contains at least `--pool-min-size` schemes, the search proceeds to the next rank. However, if the pool reaches `--pool-size` schemes earlier, the transition occurs immediately.

#### Metrics parameters

The metrics system allows tracking the progress of the search by recording statistical information about scheme ranks at each iteration. When enabled, a JSON
Lines file is created containing run metadata followed by per-iteration metrics. This feature works with both standard flip graph search and the pool strategy
(`--use-pool`), though the recorded metrics differ between the two modes.

- `--save-metrics` — enable metrics collection and saving;
- `--metrics-path PATH` — path to the metrics output file (default: `schemes/metrics.jsonl`).

##### Standard mode metrics

When metrics are enabled in standard mode, the first line of the file contains the run configuration as a JSON object:

```json
{"dimension": [4, 4, 4], "count": 128, "ring": "ZT", "copy_best_probability": 0, "seed": 1772651458, "random_walk_parameters": {"iterations": 1000000, "min_plus_iterations": 5000, "max_plus_iterations": 100000, "reset_iterations": 10000000000, "plus_diff": 4, "sandwiching_probability": 0, "reduce_probability": 0}}
```

Subsequent lines record metrics at each iteration (including iteration 0 before any flips):
```json
{"iteration": 0, "ranks": {"mean": 64, "std": 0, "min": 64, "max": 64, "min_count": 128, "max_count": 128}, "ranks_best": {"mean": 64, "std": 0, "min": 64, "max": 64, "min_count": 128, "max_count": 128}}
{"iteration": 1, "ranks": {"mean": 63.6094, "std": 0.589417, "min": 62, "max": 64, "min_count": 7, "max_count": 85}, "ranks_best": {"mean": 63.5938, "std": 0.592103, "min": 62, "max": 64, "min_count": 7, "max_count": 83}}
{"iteration": 2, "ranks": {"mean": 63.3047, "std": 0.755796, "min": 61, "max": 64, "min_count": 2, "max_count": 60}, "ranks_best": {"mean": 63.2578, "std": 0.76287, "min": 61, "max": 64, "min_count": 2, "max_count": 56}}
{"iteration": 3, "ranks": {"mean": 63.0156, "std": 0.892542, "min": 61, "max": 64, "min_count": 9, "max_count": 43}, "ranks_best": {"mean": 62.9453, "std": 0.886607, "min": 61, "max": 64, "min_count": 10, "max_count": 37}}
```

For each iteration, the following statistics are provided:
- `iteration`: iteration number;
- `rank`: best rank achieved so far across all runners;
- `ranks`: statistics for the current ranks of all runners;
- `ranks_best`: statistics for the best rank achieved by each runner;

Each statistics object contains `mean`, `std` (standard deviation), `min` and `max` values, along with `min_count` and `max_count` indicating how many schemes currently have ranks equal to the minimum and maximum values respectively.

##### Pool mode metrics

When metrics are enabled with the pool strategy (`--use-pool`), the first line contains run configuration including pool parameters:
```json
{"dimension": [4, 4, 4], "count": 128, "ring": "ZT", "seed": 1772651458, "random_walk_parameters": {"iterations": 1000000, "min_plus_iterations": 5000, "max_plus_iterations": 10000, "reset_iterations": 10000000000, "plus_diff": 4, "sandwiching_probability": 0, "reduce_probability": 0}, "pool_parameters": {"size": 256, "min_size": 32, "max_iterations": 10, "select_strategy": "uniform"}}
```

Subsequent lines track the pool progression:
```json
{"iteration": 0, "rank": 64, "pool_size": 1, "step": 0}
{"iteration": 1, "rank": 63, "pool_size": 31, "step": 1}
{"iteration": 2, "rank": 63, "pool_size": 32, "step": 2}
{"iteration": 3, "rank": 62, "pool_size": 47, "step": 1}
```

Fields in pool mode:
- `iteration`: iteration number;
- `rank`: current target rank being searched for;
- `pool_size`: current size of the pool for the target rank;
- `step`: number of iterations spent searching for the current rank.

#### Other parameters
- `--target-rank INT` — stop search when this rank is found, 0 searches for minimum (default: `0`);
- `--copy-best-probability REAL` — probability to replace scheme with best scheme after improvement (default: `0.5`);
- `--max-improvements INT` — maximum saved recent improvements for reset sampling (default: `10`).

#### Examples
Start search from naive `4x4x4` scheme with ternary coefficients, using `128` parallel runners and `16` threads:
```bash
./flip_graph -n1 4 -n2 4 -n3 4 --ring ZT --count 128 --threads 16
```

Search binary schemes loaded from file:
```bash
./flip_graph -i input.txt --ring Z2 --count 128
```

Search ternary modular schemes from multiple-scheme file:
```bash
/flip_graph -i input.txt -m --ring Z3 --count 128
```

Use pool strategy to systematically reduce rank:
```bash
./flip_graph -i input.txt --ring ZT --use-pool --pool-size 5000
```

Track search progress with metrics in standard mode:
```bash
./flip_graph -n1 4 -n2 4 -n3 4 --ring ZT --count 128 --save-metrics --metrics-path metrics_4x4x4.jsonl
```

Track pool strategy progress with metrics:

```bash
./flip_graph -n1 4 -n2 4 -n3 4 --ring ZT --use-pool --pool-size 1000 --save-metrics --metrics-path metrics_pool_4x4x4.jsonl
```

#### Example output
```text
+-----------------------------------------------------------------------------------+
| dimension: 4x4x4            seed: 1771326467                        best rank: 52 |
| threads: 12                 flip iters: 1.0M                        iteration: 33 |
| count: 256                  reset iters: 10.0B                  elapsed: 00:10:30 |
| ring: ZT                    plus diff: 4                     improvements: 9 / 10 |
+===================================================================================+
| runner | scheme rank |   naive    |            |        flips        |    plus    |
|   id   | best | curr | complexity | iterations |  count  | available | iterations |
+--------+------+------+------------+------------+---------+-----------+------------+
|    155 |   52 |   52 |        362 |   1000.00K |   3.92K |        13 |      6.14K |
|     27 |   52 |   52 |        363 |   1000.00K |  89.09K |        11 |     91.07K |
|     31 |   52 |   52 |        363 |   1000.00K |     703 |        12 |     90.82K |
+--------+------+------+------------+------------+---------+-----------+------------+
- iteration time (last / min / max / mean): 19.68 / 12.85 / 27.88 / 19.09
```


### meta_flip_graph

`meta_flip_graph` extends flip_graph with meta-operations that can modify scheme dimensions during search. Supports all flip_graph operations plus `project`,
`extend`, `merge` and `product` operations.

After every random walk phase, each runner may invoke a meta operation with probability `--meta-probability`. Before executing a meta operation, the current
rank is compared with known best ranks. If the gap exceeds `--meta-max-rank-diff`, the runner resets to an initial scheme (not necessarily of the same size).

#### Meta operation strategies
Three strategies are available:
- `default`: with probability `0.5`, permute dimensions and attempt merge with random best scheme. If merge fails, with probability `0.5` attempt `projection`
  to random dimension, otherwise perform `extension`.
- `proj`: performs only random projection operations.
- `ext`: performs only random extension operations.

#### Main parameters
Same as `flip_graph`.

#### Initialization parameters
Same as `flip_graph`.

#### Random walk parameters
Same as `flip_graph`.

#### Meta parameters
- `--meta-probability REAL` — probability of call meta operations (default: `0`);
- `--meta-strategy {default, proj, ext}` — strategy of meta operations (default: `default`);
- `--meta-min-dimension INT` — minimum dimension for projection (default: `2`);
- `--meta-max-dimension INT` — maximum dimension for merge/extend (default: `16`);
- `--meta-max-rank INT` — maximum rank for merge/extend (default: `350`);
- `--meta-max-rank-diff INT` — reset threshold relative to known rank (default: `10`).

#### Additional parameters
- `--improve-ring {Z2, ZT, Q}` — save only schemes improving known rank (saves all by default);
- `--int-width {16, 32, 64, 128}` — integer bit width, determines maximum matrix elements (default: `64`).

#### Example
Search binary schemes with meta-operations, dimensions varying from 4 to 10, rank limit 400:
```bash
./meta_flip_graph -i input.txt -m --ring Z2 --meta-probability 0.5 --meta-strategy default \
  --meta-min-dimension 4 --meta-max-dimension 10 --meta-max-rank 400
```


## Supporting tools

### find_alternative_schemes
Generates alternative schemes of the same dimensions as the input scheme using `flip` / `plus` / `split` / `sandwiching` operations. Useful for expanding
existing schemes for independent analysis or as broader initialization for new searches.

#### Parameters
- `--ring {ZT, Z2, Z3}` — coefficient ring (default: `ZT`);
- `--max-count INT` — number of alternatives to generate (default: `10K`);
- `--input-path PATH` — path to file with input scheme (required);
- `--output-path PATH` — output directory for alternative schemes (default: `schemes`);
- `--target-rank INT` — keep only schemes with this rank, if not specified, the rank of the input scheme is used as the target;
- `--plus-probability REAL` — probability of plus (default: `0.001`);
- `--plus-diff INT` — rank difference allowed for plus (default: `2`);
- `--sandwiching-probability REAL` — sandwiching probability (default: `0`);
- `--seed INT` — random seed (default: `0`).

The tool expects a single scheme in the input file. If `--target-rank` is specified, only schemes with that rank are saved (can be higher or lower than
initial rank).

#### Example
Find 128 alternative ternary schemes for `3x3x8` with rank `56`:
```bash
./find_alternative_schemes -i inputs/3x3x8_m56_ZT.txt --max-count 128 --ring ZT
```

#### Example output
```text
+-----------+-------------+------------+-----------------+
| iteration | alternative | complexity | available flips |
+-----------+-------------+------------+-----------------+
|         1 |           1 |        510 |              15 |
|         2 |           2 |        503 |              15 |
|         3 |           3 |        517 |              15 |
...
```


### lift
Attempts to lift schemes from modular rings (`Z₂` / `Z₃`) to general rings (`ZT` / `Z` / `Q`) using Hensel lifting followed by rational reconstruction. For `Z₂`
schemes, lifts through `Z₄ → Z₈ → Z₁₆ → ...`; for `Z₃` schemes, through `Z₉ → Z₂₇ → Z₈₁ → ...`. After each lifting step, attempts rational reconstruction to
obtain rational/integer coefficients.

#### Parameters
- `--ring {Z2, Z3}` — source ring (required);
- `--input-path PATH` — path to file with input scheme(s);
- `--output-path PATH` — output directory for lifted schemes (default: `schemes`);
- `--multiple` — input file contains multiple schemes;
- `--steps INT` — number of lifting steps; (default: `10`)
- `--canonize` — canonize reconstructed schemes;
- `--threads INT` — OpenMP threads;
- `--int-width {16, 32, 64, 128}` — integer width (default: `64`).

#### Example
Lift multiple binary schemes:
```bash
./lift -i input.txt -m --ring Z2 --threads 1 --steps 10
```

Example output lifting 10 binary schemes:
```text
+--------+-----------+------+----------------------------+-------+--------------+
| scheme | dimension | rank |           status           | steps | elapsed time |
+--------+-----------+------+----------------------------+-------+--------------+
|      1 |    4x5x10 |  148 |        reconstructed in ZT |     1 |         3.88 |
|      2 |     6x6x6 |  153 |        reconstructed in ZT |     1 |         4.77 |
|      3 |     4x7x7 |  144 |         reconstructed in Z |     6 |        27.02 |
|      4 |     5x6x7 |  150 |        reconstructed in ZT |     1 |         4.07 |
|      5 |     6x7x8 |  239 |         reconstructed in Q |     8 |     00:06:46 |
|      7 |     5x5x7 |  127 |        reconstructed in ZT |     1 |         2.32 |
|      8 |     4x6x7 |  123 |        reconstructed in ZT |     1 |         2.10 |
|      9 |     4x5x7 |  104 |        reconstructed in ZT |     1 |         0.86 |
|     10 |     6x6x7 |  183 |         reconstructed in Q |     6 |     00:02:54 |
+--------+-----------+------+----------------------------+-------+--------------+
- elapsed time (total / mean): 00:07:25 / 00:01.02
```


### validate_schemes
Verifies that schemes satisfy Brent equations, determining whether they are valid matrix multiplication schemes.
Supports both integer and fractional coefficients.

#### Parameters
- `--input-path PATH` — path to input file with scheme(s) (required)
- `--multiple` — input file contains multiple schemes;
- `--format {int, frac}` — integer or rational input format (required);
- `--show-ring` — show detected ring;
- `--show-coefficients` — show unique coefficient values.

For fractional coefficients (`--format frac`), each value must be written as a pair of integers (numerator and denominator), even for integer values
(e.g., `7 1` for 7).

#### Example

```bash
./validate_schemes -i input.txt -m --format frac --shor-ring --show-coefficients
```

Example output:
```text
Start checking 10 schemes in "input.txt"
- correct scheme 1 / 10 (4, 4, 7: 85), ring: ZT, values: {-1, 0, 1}
- correct scheme 2 / 10 (2, 2, 12: 42), ring: ZT, values: {-1, 0, 1}
- correct scheme 3 / 10 (2, 2, 4: 14), ring: ZT, values: {-1, 0, 1}
- correct scheme 4 / 10 (2, 4, 7: 45), ring: ZT, values: {-1, 0, 1}
- correct scheme 5 / 10 (3, 3, 8: 56), ring: ZT, values: {-1, 0, 1}
- correct scheme 6 / 10 (3, 4, 7: 64), ring: ZT, values: {-1, 0, 1}
- correct scheme 7 / 10 (3, 3, 9: 63), ring: Z, values: {-1, -2, -3, -4, 0, 1, 2, 3}
- correct scheme 8 / 10 (3, 4, 6: 54), ring: Q, values: {-1, -1/2, 0, 1, 1/2, 2}
- correct scheme 9 / 10 (2, 2, 8: 28), ring: Z, values: {-1, -2, 0, 1, 2}
- correct scheme 10 / 10 (4, 5, 6: 90), ring: Z, values: {-1, -2, 0, 1, 2}

All 10 schemes are correct
```

### optimize_scheme
Uses flip operations to optimize schemes according to a selected metric. By default, it minimizes naive additive complexity (maximizes the number of zero
coefficients), but can also optimize for the number of available flips.

#### Parameters
- `--ring {ZT, Z2, Z3}` — coefficient ring (default: `ZT`);
- `--input-path PATH` — path to input file with scheme(s) (required);
- `--output-path PATH` — output directory for optimized schemes (default: `schemes`);
- `--multiple` — input file contains multiple schemes;
- `--count INT` — number of runners (default: `8`);
- `--threads INT` — OpenMP threads;
- `--metric {complexity, flips}` — metric to optimize (default: `complexity`);
- `--format {txt, json}` — output format (default: `txt`);
- `--flip-iterations INT` — flips per report (default: `100K`);
- `--plus-probability REAL` — probability of `plus` operation (default: `0.01`);
* `--plus-diff INT` — allowed rank difference for `plus` operation (default: `2`);
- `--maximize` — maximize instead of minimize;
- `--copy-best-probability REAL` — probability to replace scheme with best scheme after improvement (default: `0.5`);
- `--max-no-improvements INT` — termination threshold (default: `3`).

#### Example
Minimizing a `3x3x8` rank `56` ternary scheme:
```bash
./optimize_scheme -c 16 -i inputs/3x3x8_m56_ZT.txt --top-count 3
```

Maximize the number of available flips for the same scheme:
```bash
./optimize_scheme -c 16 -i inputs/3x3x8_m56_ZT.txt --metric flips --maximize
```

Example output complexity minimizing a `3x3x8` rank `56` ternary scheme:
```text
+----------------------------------+
| dimension       rank        ring |
|     3x3x8         56          ZT |
+----------------------------------+
| count: 16 (12 threads)           |
| seed: 1771334776                 |
| best complexity: 439             |
| iteration: 4                     |
| elapsed: 1.81                    |
+==================================+
| runner |    scheme complexity    |
|   id   |    best    |    curr    |
+--------+------------+------------+
| 5      | 439        | 527        |
| 2      | 440        | 439        |
| 1      | 440        | 439        |
+--------+------------+------------+
- iteration time (last / min / max / mean): 0.52 / 0.42 / 0.52 / 0.45
```


## File Formats
### Single Scheme Format

- First line: dimensions and rank: `n1 n2 n3 rank`;
- Second line: `U` coefficients;
- Third line: `V` coefficients;
- Fourth line: `W` coefficients.

Example file with `2x2x3` scheme with `11` multiplications:
```text
2 2 3 11
1 0 0 1 1 0 0 0 0 0 0 1 -1 1 0 0 0 1 0 1 1 0 1 0 0 0 1 -1 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1
1 0 0 0 1 0 0 1 0 0 1 0 1 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1 -1 0 -1 1 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 1
1 0 0 1 0 0 0 0 1 -1 0 0 -1 1 0 0 0 0 1 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 1
```

### Multiple Schemes Format
First line contains the number of schemes, followed by that many schemes in single scheme format.

Example with three schemes:
```text
3
2 2 2 7
1 0 0 1 1 0 0 0 0 1 0 -1 -1 0 1 0 0 0 0 1 1 1 0 0 0 0 1 1
1 0 0 1 0 1 0 -1 0 0 1 1 1 1 0 0 -1 0 1 0 0 0 0 1 1 0 0 0
1 0 0 1 0 0 1 1 1 0 0 0 0 0 0 1 1 1 0 0 -1 0 1 0 0 1 0 -1
2 2 3 11
0 0 0 1 0 0 0 1 0 0 1 0 0 0 1 1 0 1 0 0 0 1 0 -1 1 0 0 0 1 0 0 0 1 0 0 1 1 0 -1 0 1 1 0 0
0 0 0 0 0 1 0 -1 0 0 1 0 0 0 1 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 1 0 0 0 1 0 0 -1 0 0 0 1 0 1 0 0 -1 -1 0 0 0 0 0 0 0 1 0 0
0 0 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 1 0 1 0 -1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 1 1 0 0 0 0 1 0 0 0 0 1 0 -1 0 0 0
2 2 3 11
0 0 1 0 0 0 1 1 0 1 1 0 0 1 0 0 1 1 0 0 1 0 -1 0 1 0 0 0 0 1 0 -1 0 0 1 -1 0 1 0 0 0 0 0 1
0 1 0 0 -1 0 0 0 0 0 1 0 0 0 1 0 1 0 0 0 1 0 0 -1 0 0 1 0 0 0 0 1 1 0 0 0 1 0 -1 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 1 -1 0 0 1 0 0 1 0 0
0 0 1 1 0 0 0 0 0 1 0 -1 0 0 1 0 0 1 0 0 0 0 -1 -1 1 0 -1 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 -1 0 1 0 0 0 0 -1 0 0 0 0 0 0 1 0 0 0 0
```

When reading multiple schemes, use the `-m` flag with tools that support it.

## Important Notes

- When reading `Z₂` / `Z₃` schemes from files, coefficients are automatically reduced modulo 2 or 3.
- For ternary schemes, input coefficients must already be in `{-1, 0, 1}` — automatic conversion is not possible.
- Cannot specify both naive dimensions (`-n1`, `-n2`, `-n3`) and input file (`--input-path`) simultaneously.
- All tools verify scheme correctness by default unless `--no-verify` is specified.


## Citation
If you use this software in your research, please cite:

```bibtex
@article{perminov2026fast,
    title={Fast Matrix Multiplication in Small Formats: Discovering New Schemes with an Open-Source Flip Graph Framework},
    author={Perminov, Andrew I},
    journal={arXiv preprint arXiv:2603.02398},
    url={https://arxiv.org/abs/2603.02398},
    year={2026}
}
```

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
