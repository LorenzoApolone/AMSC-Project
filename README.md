# Particle Swarm Optimization Framework

This project implements two variants of the **Particle Swarm Optimization (PSO)** algorithm in C++ for solving global optimization problems.  
It includes **serial and parallel versions** (MPI and OpenMP), an **automated testing system**, and tools for **performance analysis** on standard benchmark functions.

## Introduction

This project was developed as part of the "Advanced Mathod for Scientific Computing" course (2025/2026) teached by Luca Formaggia in "Politecnico di Milano ".  
Its goal is to explore the Particle Swarm Optimization algorithm.  
The framework allows easy testing of multiple mathematical benchmark functions and comparison between different PSO variants on both local and distributed systems.

This project is made by the following students: 

- Lorenzo Apolone
- Francesco Calzona
- Alessandro Masini
- Matteo Parimbelli
- Giovanni Vaccaro


## Objectives

- Implement the **serial PSO algorithm** and   the **parallel PSO algorithm**
- Test the algorithm on **benchmark functions** such as Sphere, Ellipsoid, Sum of Different Powers, DropWave, Weierstrass, Alpine1, and others.  
- Develop **parallel implementations** using **OpenMP** and **MPI**, and analyze **speed-up**, **efficiency**, and **convergence**.  


## Requirements

- **Compiler:** GCC or Clang with C++17 support  
- **Parallelism:** OpenMPI (for distributed version) or OpenMP (for shared-memory version)  
- **Tools:**  
  - GNU Make or CMake  
  - Python 3 + Matplotlib (optional, for plots)  


## Installation

Clone the repository:

```bash
git clone https://github.com/AMSC-25-26/swarm-search-2-swarmsearch

```


## Execution

Two type of executions are supported: the serial one and the parallel one. Both can be compiled simultaneously.

### Compilation

From the folder ./swarm-search-2-swarmsearch/:

To compile both files:

```cd src
make 
```
Or compiling individually

```make serial
make parallel 
```

### Parallel Execution
To execute it starting from: 

```
mpirun -np <Num of prcessor>  ./main_parallel <Dim> <Number of particel> <Max iterations> <Error supported>
```
Example: 

```
mpirun -np 4 ./main_parallel 10 100 1000 0.001
```

### Serial Execution

If you want to use the serial:  


```
./main_serial <Dim> <Number_of_particles> <Max_iterations> <Error_tolerance>
```

Example: 

```
mpirun ./main_serial 10 100 1000 0.001
```


## Project Structure 
```bash
.
|-- Bibliography
|   |-- data-07-00046-v2.pdf
|   |-- s11831-021-09694-4.pdf
|   |-- s12065-019-00210-z.pdf
|   `-- s13369-018-03713-6.pdf
|-- Description_SwarmOptimization.pdf
|-- README.md
|-- guides
|   `-- massive_testing.md
|-- poetry.lock
|-- pyproject.toml
`-- src
    |-- Makefile
    |-- __pycache__
    |-- functions.cpp
    |-- functions.o
    |-- interfaces.cpp
    |-- interfaces.hpp
    |-- interfaces.o
    |-- main_parallel
    |-- main_parallel.cpp
    |-- main_parallel.o
    |-- main_serial.cpp
    |-- make.dep
    |-- methods.hpp
    |-- particle.hpp
    |-- plots.py
    |-- pso_mpi.cpp
    |-- pso_mpi.o
    |-- pso_serial.cpp
    |-- run_plots.py
    |-- run_tests.sh
    |-- tests
    `-- utilities.hpp
```
## functions.cpp

This module implements a **suite of continuous optimization benchmarks** as concrete classes derived from  Bibliography/data-07-00046-v2.pdf
Each class defines:
- a **name** (for logging/output),
- **domain bounds** `(lower, upper)`,
- the **known global minimizer** `x*`, and
- the scalar **objective** `double value(const std::vector<double>& x) const`.

The set covers convex, ill-conditioned, separable/non-separable, and strongly multimodal landscapes to stress different PSO behaviors.

### Implemented functions
#### 1. Sphere
$$
f(x) = \sum_{i=1}^{D} x_i^2
$$

#### 2. Ellipsoid
$$
f(x) = \sum_{i=1}^{D} i \, x_i^2
$$

#### 3. Sum of Different Powers
$$
f(x) = \sum_{i=1}^{D} |x_i|^{i+1}
$$

#### 4. Quintic Function
$$
f(x) = \sum_{i=1}^{D} \big(x_i^5 - 3x_i^4 + 4x_i^3 + 2x_i^2 - 10x_i - 4\big)
$$

#### 5. Drop-Wave
$$
f(x) = 1 - \frac{1 + \cos\left(12\sqrt{\sum_{i=1}^{D} x_i^2}\right)}{0.5 \sum_{i=1}^{D} x_i^2 + 2}
$$

#### 6. Weierstrass
$$
f(x) = \sum_{i=1}^{D} \sum_{k=0}^{k_{max}} \big[a^k \cos(2\pi b^k (x_i + 0.5))\big]- D \sum_{k=0}^{k_{max}} \big[a^k \cos(\pi b^k)\big]
$$
where $$ a = 0.5, \, b = 3, \, k_{max} = 20 $$

#### 7. Alpine1
$$
f(x) = \sum_{i=1}^{D} |x_i \sin(x_i) + 0.1x_i|
$$

#### 8. Ackley
$$
f(x) = -20 \exp\left(-0.2 \sqrt{\frac{1}{D}\sum x_i^2}\right)-\exp\left(\frac{1}{D}\sum \cos(2\pi x_i)\right) + 20 + e
$$

#### 9. Griewank
$$
f(x) = \frac{1}{4000}\sum_{i=1}^{D} x_i^2 - \prod_{i=1}^{D} \cos\left(\frac{x_i}{\sqrt{i}}\right) + 1
$$

#### 10. Rastrigin
$$
f(x) = \sum_{i=1}^{D} [x_i^2 - 10 \cos(2\pi x_i) + 10]
$$

#### 11. HappyCat
$$
f(x) = \left(\frac{\|x\|^2 - D}{4} \right)^2 + \frac{1}{D}\sum x_i + 0.5
$$

#### 12. HGBat
$$
f(x) = \sqrt{|\|x\|^2 - (\sum x_i)^2|} + \frac{0.5(\|x\|^2 + (\sum x_i)^2)}{D} + 0.5
$$

#### 13. Rosenbrock
$$
f(x) = \sum_{i=1}^{D-1} [100(x_{i+1} - x_i^2)^2 + (x_i - 1)^2]
$$

#### 14. High Conditioned Elliptic
$$
f(x) = \sum_{i=1}^{D} (10^6)^{\frac{i-1}{D-1}} x_i^2
$$

#### 15. Discus
$$
f(x) = 10^6 x_1^2 + \sum_{i=2}^{D} x_i^2
$$

#### 16. Bent Cigar
$$
f(x) = x_1^2 + 10^6 \sum_{i=2}^{D} x_i^2
$$

#### 17. Perm D, Beta (PermdbFunc)
$$
f(x) = \sum_{k=1}^{D} \left[\sum_{i=1}^{D} \left(\frac{i^k + \beta}{D}(x_i^k - \frac{1}{i^k})\right)\right]^2
$$
where $$ \beta = 0.5 $$

#### 18. Schaffer F7
$$
f(x) = \left(\frac{1}{D-1}\sum_{i=1}^{D-1} \sqrt{s_i}(\sin(50 s_i^{0.2}) + 1)\right)^2
$$
where $$ s_i = x_i^2 + x_{i+1}^2 $$

#### 19. Expanded Schaffer F6
$$
f(x) = \sum_{i=1}^{D-1} \left[0.5 + \frac{\sin^2(\sqrt{x_i^2 + x_{i+1}^2}) - 0.5}{(1 + 0.001(x_i^2 + x_{i+1}^2))^2}\right]
$$

#### 20. Rotated Hyper-Ellipsoid
$$
f(x) = \sum_{i=1}^{D} \sum_{j=1}^{i} x_j^2
$$

#### 21. Schwefel
$$
f(x) = 418.9829D - \sum_{i=1}^{D} x_i \sin(\sqrt{|x_i|})
$$

#### 22. Sum of Different Powers 2
$$
f(x) = \sum_{i=1}^{D} |x_i|^{i+2}
$$

#### 23. Xin-She Yang 1
$$
f(x) = \sum_{i=1}^{D} \epsilon_i |x_i|^i, \quad \epsilon_i \sim U(0,1)
$$

#### 24. Schwefel 2.21
$$
f(x) = \max_i |x_i|
$$

#### 25. Schwefel 2.22
$$
f(x) = \sum_{i=1}^{D} |x_i| + \prod_{i=1}^{D} |x_i|
$$

#### 26. Salomon
$$
f(x) = 1 - \cos(2\pi \sqrt{\sum x_i^2}) + 0.1\sqrt{\sum x_i^2}
$$

#### 27. Modified Ridge
$$
f(x) = x_1^2 + 100\sqrt{\sum_{i=2}^{D} x_i^2}
$$

#### 28. Zakharov
$$
f(x) = \sum_{i=1}^{D} x_i^2 + \left(\sum_{i=1}^{D} 0.5 i x_i\right)^2 + \left(\sum_{i=1}^{D} 0.5 i x_i\right)^4
$$

#### 29. Modified Xin-She Yang 3
$$
f(x) = \exp\left(-\sum_{i=1}^{D} \sin(x_i^2)\right) - 2\exp\left(-\sum_{i=1}^{D} (x_i - \pi)^2\right)
$$

#### 30. Modified Xin-She Yang 5
$$
f(x) = \sum_{i=1}^{D} \sin^6(x_i) + 0.1 \sum_{i=1}^{D} (x_i - \pi)^2
$$

## Output and Results Management (`interfaces.cpp`)

This component handles **result reporting and data storage** for each PSO optimization run.  
It is responsible for both **terminal output** (human-readable summaries) and **structured file output** (for later analysis or plotting).

---

### **Main Responsibilities**

1. **Display concise run summaries**
   - Function: `void OutputObject::terminal_info()`
   - Prints to the terminal key information about a completed optimization:
     - Function name  
     - Problem dimension  
     - Number of particles  
     - Number of cores  
     - Final convergence value (last Δx)  
     - Execution time  
     - Total number of iterations  

   Used for **quick inspection** of single runs.  
   For large-scale experiments, it’s preferable to use file export instead.

---

2. **Store results in structured directories**
   - Function: `void OutputObject::output_to_file()`
   - Automatically generates a hierarchical directory tree:

     ```
     tests/<function_name>/<dimension>/<n_points>/<n_cores>/
     ```

   - Each test run is saved as:
     ```
     test_<N>.txt
     ```
     where `<N>` is an auto-incremented index determined by scanning existing files (via `get_max_test_number()`).

---

3. **File contents**

   Each output file stores one line per iteration of the PSO run with columns:



## Serial PSO (`pso_serial.cpp`)

Implements the **sequential Particle Swarm Optimization (PSO)** algorithm.

- Function: `pso_serial(const TestFunction& f, int d, const StopCriterion& stop, int n_points)`
- Uses time-varying inertia (`0.9 → 0.4`), personal/global best updates, and boundary clamping.
- Stops when reaching the **maximum iterations** or **convergence tolerance**.

**Key steps:**
1. Random initialization of particles within `[LB, UB]`.
2. Iterative update of velocity and position:

## Parallel PSO (MPI) (`pso_mpi.cpp`)

Implements a **parallel PSO** using MPI with global-best synchronization.

- Function: `pso_mpi(const TestFunction& f, int d, const StopCriterion& stop, int n_points)`
- Work split: particles are partitioned across ranks (`local_n` balanced with remainder).
- Initialization: per-rank RNG (`seed = rank + 42`), positions in `[LB,UB]`, small random velocities.
- Iteration:
  - Velocity/position update with linearly decreasing inertia (`0.9 → 0.4`) and boundary clamping.
  - Local update of `pbest` and tentative `gbest`.
  - **Global reduction** of best value via `MPI_Allreduce` with `MPI_MINLOC`, then broadcast of the corresponding `gbest_pos`.
- Stopping:
  - Rank 0 computes error `f.error(gbest_pos)`, appends to `history`, evaluates `stop.should_stop(...)`, broadcasts stop signal.
- Timing & output:
  - Wall time measured with `MPI_Wtime()`.
  - Returns `OutputObject` with best solution, convergence history (from rank 0), total time, iterations, and `n_cores = size`.

## Main Parallel (`main_parallel.cpp`)

This file launches the **parallel execution of the PSO algorithm** on all implemented benchmark functions.  
It initializes MPI, sets problem parameters from command-line arguments, and runs `pso_mpi` for each test function in sequence.

## Main Serial (`main_serial.cpp`)

This file launches the **serial execution of the PSO algorithm** on all implemented benchmark functions.  
It initializes and sets problem parameters from command-line arguments, and runs `pso_serial` for each test function in sequence.
