/*
//creazione Dummy object

StopCriterion stopcriterion(100, 1e-6);
int dim = 5;
std::vector<double> x_best(dim, 0.0);
std::vector<double> x_star(dim, 0.0);
std::vector<double> conv_history = {1.0, 0.5, 0.25, 0.1};

OutputObject out(
    "DummyFunction",   // function_name
    dim,                 // dimension
    100,               // n_points
    x_best,            // x_best
    x_star,            // x_star
    0.0,               // f_val
    conv_history,      // convergence history
    4,                 // n_cores
    0.01,              // execution time (s)
    10,                // iterations
    stopcriterion      // stop criterion
);
*/






  /*  
    // --- stampa ID globali delle particelle (ordinato per rank) ---
for (int r = 0; r < size; ++r) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == r) {
        std::cout << "Rank " << rank
                  << " gestisce local_n=" << local_n
                  << " particelle: ";
        for (int i = 0; i < local_n; ++i) {
            int gid = displs[rank] + i; // ID globale
            std::cout << gid << " ";
        }
        std::cout << "\n" << std::flush;
    }
}
MPI_Barrier(MPI_COMM_WORLD);
                  
*/   