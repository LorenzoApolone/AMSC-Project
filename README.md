# AMSC Particle Swarm Optimization Project

This repository contains a set of Particle Swarm Optimization implementations
used to compare serial, MPI-parallel, topology-based, DMS-PSO-HS and CPSO
approaches on a common benchmark suite.

The code is organized around shared benchmark interfaces in `src/interfaces/`
and a shared function library in `src/interfaces/functions.hpp`. Each solver returns
standard output objects or run artifacts so that convergence, execution time and
communication costs can be compared across methods.

## Repository Layout

- `src/`: common interfaces, benchmark functions and shared build files.
- `src/pso/`: standard serial and MPI PSO implementation and drivers.
- `src/topology/`: topology-based PSO, with classic, small-world, scale-free and
  random communication graphs.
- `src/dpso/`: DMS-PSO-HS serial and MPI implementations.
- `src/cpso/`: Cooperative PSO implementation with serial and MPI solvers.
- `src/benchmarks/`: benchmark launchers, raw results, analysis scripts and plots.
- `Bibliography/`: reference papers and source material.

## Quick Start: Build and Run

The commands below assume a Linux, WSL or HPC login-node environment. Each
section states the directory where the command should be launched.

### 1. Install or load prerequisites

Required tools:

- C++17 compiler;
- `make`;
- MPI compiler and launcher (`mpic++`, `mpirun`);
- Python 3 for benchmark orchestration scripts.

Optional documentation tools:

- `doxygen`;
- Graphviz `dot`, used by Doxygen for class and call graphs.

On Ubuntu/WSL, the usual package set is:

```bash
sudo apt update
sudo apt install build-essential make openmpi-bin libopenmpi-dev python3 doxygen graphviz
```

On an HPC cluster, load the equivalent compiler, MPI and Python modules before
building. For example:

```bash
module load gcc
module load openmpi
module load python
```

### 2. Build all local executables

```bash
cd AMSC-Project
make -C src
make -C src/topology
make -C src/cpso
```

This creates:

- `src/main_serial` and `src/main_parallel` for standard PSO;
- `src/main_dpso` and `src/main_dpso_serial` for DMS-PSO-HS;
- `src/topology/topology_parallel` and `src/topology/topology_serial`;
- `src/cpso/test_cpso` and `src/cpso/test_cpso_parallel`.

To remove object files and dependency files:

```bash
cd AMSC-Project
make -C src clean
make -C src/topology clean
make -C src/cpso clean
```

To also remove the top-level `src` executables:

```bash
cd AMSC-Project
make -C src distclean
```

### 3. Run quick smoke tests

These examples use small dimensions and iteration counts so that each executable
can be checked quickly. Increase the parameters for full experiments.

Standard PSO:

```bash
cd src

./main_serial 8 64 1000 1e-4 123
mpirun -np 4 ./main_parallel 8 64 1000 1e-4
```

DMS-PSO-HS:

```bash
cd src

./main_dpso_serial 8 64 1000 1e-4 dpso/params.txt 123
mpirun -np 4 ./main_dpso 8 64 1000 1e-4 dpso/params.txt 123
```

Topology-based PSO:

```bash
cd src/topology

./topology_serial 8 64 1000 1e-4 123
mpirun -np 4 ./topology_parallel 8 64 1000 1e-4 123
```

Cooperative PSO:

```bash
cd src/cpso

./test_cpso 8 4 16 1000 123
mpirun -np 4 ./test_cpso_parallel 8 4 16 1000 50 50 123
```

For CPSO communication ablation runs, disable the greedy merge fallback with:

```bash
CPSO_MPI_DISABLE_GREEDY_MERGE=1 mpirun -np 4 ./test_cpso_parallel 8 4 16 1000 50 50 123
```

### 4. Run benchmark batteries

The benchmark scripts build the needed executables by default, generate a case
matrix and store raw results under ignored `results/` folders. Use `--dry-run`
first to inspect the commands.

### CPSO Benchmark Batteries

Available batteries via the `--battery <NAME>` argument:
- **`validation`**: Short correctness check to verify the pipeline. Runs `dim=64` and `dim=128` with `k=28` subswarms on serial, np=1, np=4, and np=28.
- **`comparable`**: Inter-method comparison. Strong scaling with a fixed total problem (`dim=32, 64, 128`), `k=32`, and 8 particles per swarm. Varies MPI processes.
- **`cpso`**: CPSO specific runs. Includes weak scaling tests (constant particles per process, constant dimension per process, constant local load) and a dimension sweep up to `dim=256` for `np=8`.
- **`communication`**: Ablation study on MPI communication policies. Compares baseline greedy merge versus no-greedy fallback for `dim=64, 128`.
- **`appendix`**: Transition observation runs. Fixed `k=28` varying MPI processes to observe the transition towards `np=28`.
- **`all`**: Runs all the above batteries.

```bash
cd AMSC-Project
bash src/benchmarks/cpso/run_cpso_benchmarks.sh --battery all --dry-run
bash src/benchmarks/cpso/run_cpso_benchmarks.sh --battery all --seeds 123
```

### DMS-PSO-HS (DPSO) Benchmark Batteries

Available batteries via the `--battery <NAME>` argument:
- **`strong`**: Strong scaling benchmark. Fixed total problem (dim=32, 64, 128) with 256 particles. Varies MPI processes.
- **`weak_ppc`**: Weak scaling with constant particles per process. `dim=64` fixed, total particles scale with MPI processes.
- **`weak_local`**: Weak scaling with constant local load. Each process gets a fixed portion of dimensions and particles. Total dimensions and particles scale with MPI processes.
- **`dim_sweep`**: Dimension sweep with a fixed number of MPI processes (`np=8`) and fixed total particles (`128`).
- **`serial`**: Serial baseline. Runs `main_dpso_serial` with the same configurations as the strong scaling and dimension sweep for direct comparison.
- **`all`**: Runs all the above specific DPSO tests.

```bash
cd AMSC-Project
bash src/dpso/run_benchmarks.sh --battery all
```

### Topology Benchmark Batteries

The topology benchmarks are launched sequentially from a predefined input list. It executes the following test campaigns:
- **Strong Scaling**: Fixed total problem size for `dim=32`, `dim=64`, and `dim=128` with 256 particles. Runs serial and MPI from np=1 up to np=28.
- **Weak Scaling (Constant Particles)**: `dim=64`, particles scale with MPI processes (`32 * np`).
- **Weak Scaling (Constant Dimension)**: Total particles fixed at 224, dimension scales with MPI processes (`8 * np`).
- **Weak Scaling (Constant Local Load)**: Both dimensions (`4 * np`) and particles (`32 * np`) scale with MPI processes.
- **Dimension Sweep**: Fixed `np=8`, dimension varies from 16 to 256.

```bash
cd AMSC-Project
bash src/benchmarks/topology/run_topology_benchmarks.sh
```

### 5. Generate and view the Doxygen documentation

The first command generates
`docs/doxygen/html/`; the second command serves that generated HTML locally.

```bash
cd AMSC-Project
doxygen Doxyfile
python3 -m http.server 8765 --bind 127.0.0.1 --directory docs/doxygen/html
```

Then open:

```text
http://127.0.0.1:8765/index.html
```

The generated `docs/doxygen/` directory is ignored by Git. Commit the Doxygen
configuration and theme files, but not the generated HTML output.

## Standard PSO

The base implementation is the reference Particle Swarm Optimization solver.
Particles move in the full `D`-dimensional search space and update their
velocity from three terms: inertia, personal best attraction and swarm-level
social attraction.

The serial version is implemented in `src/pso/pso_serial.cpp`, while the MPI version
is implemented in `src/pso/pso_mpi.cpp`. Both rely on the shared `Particle`,
`TestFunction`, `OutputObject` and `StoppingCriteriaManager` abstractions.

Typical build and run:

```bash
cd src
make

./main_serial 32 256 10000 0.0001 123
mpirun -np 4 ./main_parallel 32 256 10000 0.0001
```

Serial command-line arguments are:

```text
<dimension> <number_of_particles> <max_iterations> <target_error> [seed]
```

Parallel command-line arguments are:

```text
<dimension> <number_of_particles> <max_iterations> <target_error>
```

## Variant 1: Topology-Based PSO

The topology-based variant changes the way particles exchange information.
Instead of always using the global best particle as the social reference, each
particle communicates through a graph neighborhood. This makes the information
flow slower and more structured, which can help preserve diversity.

The supported topologies are:

- classic global PSO;
- small-world network;
- scale-free network;
- random Erdos-Renyi network.

The main files are in `src/topology/`:

- `create_network.hpp/.cpp`: topology generators;
- `pso_topology.hpp`, `pso_topology.cpp`: MPI topology-based solver;
- `pso_serial_topology.cpp`: serial topology-based solver;
- `main_topology.cpp`: parallel benchmark entry point;
- `main_topology_serial.cpp`: serial benchmark entry point.

Build and run:

```bash
cd src/topology
make

mpirun -np 4 ./topology_parallel 64 512 10000 0.0001 456
./topology_serial 64 512 10000 0.0001 456
```

## Variant 2: DMS-PSO-HS

The DMS-PSO-HS implementation combines Dynamic Multi-Swarm PSO with a Harmony
Search inspired refinement phase. The swarm is periodically regrouped into
smaller dynamic sub-swarms, so particles do not keep the same neighborhood for
the whole run.

The Harmony Search phase introduces additional exploration through harmony
memory parameters such as HMCR and PAR. These parameters are collected in
`DPSOParameters` inside `src/dpso/methods_dpso.hpp`.

Main files:

- `src/dpso/dpso.cpp`: MPI implementation;
- `src/dpso/dpso_serial.cpp`: serial implementation;
- `src/dpso/main_dpso.cpp`: MPI benchmark driver;
- `src/dpso/main_dpso_serial.cpp`: serial benchmark driver;
- `src/dpso/methods_dpso.hpp`: public interface and parameters.
- `src/dpso/params.txt`: optional per-function parameter configuration.

Build and run:

```bash
cd src
make dpso

mpirun -np 4 ./main_dpso 32 256 10000 0.0001
./main_dpso_serial 32 256 10000 0.0001
```

To use the provided DPSO parameter file from the `src` directory:

```bash
mpirun -np 4 ./main_dpso 32 256 10000 1e-6 dpso/params.txt 123
./main_dpso_serial 32 256 10000 1e-6 dpso/params.txt 123
```

## Variant 3: CPSO

The CPSO implementation decomposes the original `D`-dimensional search space
into `k` sub-swarms. Each sub-swarm owns only a subset of coordinates and
optimizes those coordinates through a local topology-aware PSO step. Partial
solutions are evaluated through a shared `ContextVector`, so every local
candidate is still tested as part of a full objective vector.

The serial solver accepts improvements greedily: each sub-swarm is processed in
sequence and can immediately update the shared context. The MPI solver uses a
batch-based protocol instead. Each rank owns a contiguous range of sub-swarms,
proposes sparse coordinate deltas, exchanges them with `MPI_Allgather`, and then
uses full, salvaged and greedy merge policies to keep all ranks synchronized on
the same accepted context.

Main files:

- `CPSOBase.hpp/.cpp`: common setup, decomposition, topologies and context init;
- `CPSOSerial.hpp/.cpp`: serial CPSO loop;
- `CPSOParallel.hpp/.cpp`: MPI CPSO loop and merge policy;
- `SubSwarm.hpp/.cpp`: local particles, memories, active dimensions and PSO step;
- `ContextVector.hpp/.cpp`: full-vector context used for partial evaluations;
- `SubSwarmOwnershipUtils.hpp`: sub-swarm to rank ownership ranges;
- `SubSwarmTopologyFactory.hpp/.cpp`: local graph construction;
- `SwarmMetrics.hpp/.cpp`: diversity and distance metrics;
- `CPSOBenchmarkUtils.hpp`: CPSO benchmark factory and output helpers.

Build and run:

```bash
cd src/cpso
make

./test_cpso 32 32 8 10000 789
mpirun -np 4 ./test_cpso_parallel 32 32 8 10000 50 50 789
```

Parallel arguments are:

```text
<dimension> <k_subswarms> <particles_per_swarm> [max_iters]
[shuffle_freq] [stagnation_patience] [seed]
```

The greedy merge fallback is enabled by default. For communication ablation
experiments it can be disabled with:

```bash
CPSO_MPI_DISABLE_GREEDY_MERGE=1 mpirun -np 4 ./test_cpso_parallel ...
```

## Benchmark Functions

All benchmark functions are defined in `src/interfaces/functions.hpp` and derive from
`TestFunction`. The current suite contains:

- `Sphere`
- `Ellipsoid`
- `SumOfDiffPowers`
- `QuinticFunction`
- `DropWave`
- `Weierstrass`
- `Alpine1`
- `Ackley`
- `Griewank`
- `Rastrigin`
- `HappyCat`
- `HGBat`
- `Rosenbrock`
- `HighCondElliptic`
- `Discus`
- `BentCigar`
- `Schafferf7Func`
- `ExpSchafferF6`
- `RotatedHyper`
- `Schwefel`
- `SumOfDifferentPowers2`
- `XinSheYang1`
- `Schwefel221`
- `Schwefel222`
- `Salomon`
- `ModifiedRidge`
- `Zakharov`
- `ModifiedXinSheYang3`
- `ModifiedXinSheYang5`
- `Levy`
- `Michalewicz`
- `Bohachevsky`
- `Powell`
- `DixonPrice`
- `StyblinskiTang`
- `Step`
- `Qing`
- `Trid`
- `Shubert`
- `Alpine2`
- `Eggholder`
- `Easom`
- `Brown`
- `Csendes`
- `Vincent`

Each function stores its dimension, domain, known reference solution and
typology tags, such as unimodal/multimodal, separable/non-separable and
differentiable/non-differentiable.
