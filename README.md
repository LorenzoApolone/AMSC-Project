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

Command-line arguments are:

```text
<dimension> <number_of_particles> <max_iterations> <target_error> [seed]
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
