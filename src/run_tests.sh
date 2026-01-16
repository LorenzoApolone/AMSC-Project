#!/bin/bash

# Usage example:
# ./run_tests.sh --n_tests 5 --dim 3 4 --n_points 50 100 --max_iter 100 --error_tol 0.001 --n_cores 2 4

param_file="./run_parameters.json"

# ------------------------------
# Parse inputs
# ------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --n_tests) n_tests="$2"; shift 2;;
        --dim) shift; dims=(); while [[ $# -gt 0 && $1 != --* ]]; do dims+=("$1"); shift; done;;
        --n_points) shift; npoints=(); while [[ $# -gt 0 && $1 != --* ]]; do npoints+=("$1"); shift; done;;
        --max_iter) shift; maxiters=(); while [[ $# -gt 0 && $1 != --* ]]; do maxiters+=("$1"); shift; done;;
        --error_tol) shift; errtols=(); while [[ $# -gt 0 && $1 != --* ]]; do errtols+=("$1"); shift; done;;
        --n_cores) shift; ncores=(); while [[ $# -gt 0 && $1 != --* ]]; do ncores+=("$1"); shift; done;;
        *) echo "Unknown parameter $1"; exit 1;;
    esac
done

# ------------------------------
# Check mandatory parameters
# ------------------------------
if [[ -z "$n_tests" || ${#dims[@]} -eq 0 || ${#npoints[@]} -eq 0 || ${#maxiters[@]} -eq 0 || ${#errtols[@]} -eq 0 || ${#ncores[@]} -eq 0 ]]; then
    echo "Missing mandatory parameters. Please provide all inputs."
    exit 1
fi

# ------------------------------
# Save parameter combinations in JSON
# ------------------------------
echo "[" > "$param_file"
first_entry=true
for dim_val in "${dims[@]}"; do
    for npts in "${npoints[@]}"; do
        for maxit in "${maxiters[@]}"; do
            for etol in "${errtols[@]}"; do
                for ncore in "${ncores[@]}"; do
                    if ! $first_entry; then
                        echo "," >> "$param_file"
                    fi
                    first_entry=false
                    echo "  {\"dim\": $dim_val, \"n_points\": $npts, \"max_iter\": $maxit, \"error_tol\": $etol,  \"n_cores\": $ncore}" >> "$param_file"
                done
            done
        done
    done
done
echo "]" >> "$param_file"
echo "Saved all parameter combinations in $param_file"

# ------------------------------
# Run C++ program for all configurations
# ------------------------------
echo "Running C++ program for all parameter combinations..."
for dim_val in "${dims[@]}"; do
    for npts in "${npoints[@]}"; do
        for maxit in "${maxiters[@]}"; do
            for etol in "${errtols[@]}"; do
                for ncore in "${ncores[@]}"; do
                    for ((i=0;i<n_tests;i++)); do
                        echo "=== RUN $i: dim=$dim_val n_points=$npts max_iter=$maxit error_tol=$etol n_cores=$ncore ==="
                        
                        echo "--- SERIAL VERSION (1 core) ---"
                        ./main_serial "$dim_val" "$npts" "$maxit" "$etol"

                        echo "--- PARALLEL VERSION ($ncore cores) ---"
                        echo "Executing: mpiexec -np $ncore ./main_parallel $dim_val $npts $maxit $etol"
                        mpiexec -np "$ncore" ./main_parallel "$dim_val" "$npts" "$maxit" "$etol"
                        echo "MPI execution completed with exit code: $?"
                        
                        echo ""
                    done
                done
            done
        done
    done
done

echo "All C++ runs completed."
