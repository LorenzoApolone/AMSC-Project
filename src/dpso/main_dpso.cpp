#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include "dpso/methods_dpso.hpp"
#include "functions.cpp"

void run_benchmark(const TestFunction& func, unsigned int dim, 
                   int max_iter, unsigned int particles_per_rank,
                   int rank, int size) {
    StopCriterion stop(max_iter, 1e-3);
    OutputObject res = pso_mpi(func, dim, stop, particles_per_rank);

    if (rank == 0) {
        double fval = func.value(res.x_best);
        double err  = func.error(res.x_best);

        std::cout << std::left << std::setw(22) << func.get_name()
                  << std::right << std::setw(6)  << dim
                  << std::setw(8)  << particles_per_rank * size
                  << std::setw(8)  << max_iter
                  << "   " << std::scientific << std::setprecision(4) << std::setw(13) << fval
                  << "   " << std::setw(13) << err
                  << "   " << std::fixed << std::setprecision(2) << std::setw(8) << res.execution_time << "s"
                  << std::endl;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double start_time = MPI_Wtime();

    std::vector<unsigned int> dims = {20};

    const unsigned int particles_per_rank = 25;
    const int base_iter = 10000;

    if (rank == 0) {
        std::cout << "\n+" << std::string(90, '-') << "+" << std::endl;
        std::cout << "|  DMS-PSO + Harmony Search  --  Benchmark Suite" << std::setw(59) << " " << "|" << std::endl;
        std::cout << "|  MPI ranks: " << std::left << std::setw(4) << size << std::setw(74) << " " << "|" << std::endl;
        std::cout << "+" << std::string(90, '-') << "+\n" << std::endl;

        std::cout << std::left  << std::setw(22) << "Function"
                  << std::right << std::setw(6)  << "Dim"
                  << std::setw(8)  << "Swarm"
                  << std::setw(8)  << "Iters"
                  << "   " << std::setw(13) << "f(x_best)"
                  << "   " << std::setw(13) << "Error"
                  << "   " << std::setw(8) << "Time"
                  << std::endl;
        std::cout << std::string(90, '-') << std::endl;
    }

    for (unsigned int dim : dims) {
        unsigned int ppr = particles_per_rank;
        int iters = base_iter;

        { Sphere f(dim);           run_benchmark(f, dim, iters, ppr, rank, size); }
        { Ellipsoid f(dim);        run_benchmark(f, dim, iters, ppr, rank, size); }
        { SumOfDiffPowers f(dim);  run_benchmark(f, dim, iters, ppr, rank, size); }
        { QuinticFunction f(dim);  run_benchmark(f, dim, iters, ppr, rank, size); }
        { DropWave f(dim);         run_benchmark(f, dim, iters, ppr, rank, size); }
        { Weierstrass f(dim);      run_benchmark(f, dim, iters, ppr, rank, size); }
        { Alpine1 f(dim);          run_benchmark(f, dim, iters, ppr, rank, size); }
        { Ackley f(dim);           run_benchmark(f, dim, iters, ppr, rank, size); }
        { Griewank f(dim);         run_benchmark(f, dim, iters, ppr, rank, size); }
        { Rastrigin f(dim);        run_benchmark(f, dim, iters, ppr, rank, size); }
        { HappyCat f(dim);         run_benchmark(f, dim, iters, ppr, rank, size); }
        { HGBat f(dim);            run_benchmark(f, dim, iters, ppr, rank, size); }
        { Rosenbrock f(dim);       run_benchmark(f, dim, iters, ppr, rank, size); }
        { HighCondElliptic f(dim); run_benchmark(f, dim, iters, ppr, rank, size); }
        { Discus f(dim);           run_benchmark(f, dim, iters, ppr, rank, size); }
        { BentCigar f(dim);        run_benchmark(f, dim, iters, ppr, rank, size); }
        { PermdbFunc f(dim);       run_benchmark(f, dim, iters, ppr, rank, size); }
        { Schafferf7Func f(dim);   run_benchmark(f, dim, iters, ppr, rank, size); }
        { ExpSchafferF6 f(dim);    run_benchmark(f, dim, iters, ppr, rank, size); }
        { RotatedHyper f(dim);     run_benchmark(f, dim, iters, ppr, rank, size); }
        { Schwefel f(dim);         run_benchmark(f, dim, iters, ppr, rank, size); }
        { SumOfDifferentPowers2 f(dim); run_benchmark(f, dim, iters, ppr, rank, size); }
        { XinSheYang1 f(dim);      run_benchmark(f, dim, iters, ppr, rank, size); }
        { Schwefel221 f(dim);      run_benchmark(f, dim, iters, ppr, rank, size); }
        { Schwefel222 f(dim);      run_benchmark(f, dim, iters, ppr, rank, size); }
        { Salomon f(dim);          run_benchmark(f, dim, iters, ppr, rank, size); }
        { ModifiedRidge f(dim);    run_benchmark(f, dim, iters, ppr, rank, size); }
        { Zakharov f(dim);         run_benchmark(f, dim, iters, ppr, rank, size); }
        { ModifiedXinSheYang3 f(dim); run_benchmark(f, dim, iters, ppr, rank, size); }
        { ModifiedXinSheYang5 f(dim); run_benchmark(f, dim, iters, ppr, rank, size); }

        if (rank == 0) std::cout << std::string(90, '-') << std::endl;
    }

    if (rank == 0) {
        std::cout << "\n=== ALL BENCHMARKS COMPLETED ===" << std::endl;
        double end_time = MPI_Wtime();
        std::cout << "Tempo totale esecuzione: " << (end_time - start_time) << " s" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
