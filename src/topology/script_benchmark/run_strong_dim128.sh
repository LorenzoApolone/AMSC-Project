#!/bin/bash

cd "$(dirname "$0")/.." || exit 1

mkdir -p script_benchmark/risultati_run

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUT_FILE="script_benchmark/risultati_run/strong_dim128_${TIMESTAMP}.txt"

echo "========================================================" | tee -a "$OUT_FILE"
echo "ESECUZIONE: STRONG SCALING - PROBLEMA TOTALE FISSO - DIM 128, K 32, PPS 8" | tee -a "$OUT_FILE"
echo "Avviata il: $(date)" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
echo "" | tee -a "$OUT_FILE"

run_cmd() {
    local cmd="$*"
    echo "--------------------------------------------------------" | tee -a "$OUT_FILE"
    echo "=> ESECUZIONE COMANDO: $cmd" | tee -a "$OUT_FILE"
    eval "$cmd" 2>&1 | tee -a "$OUT_FILE"
    echo "" | tee -a "$OUT_FILE"
}


echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 789" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"

run_cmd "mpirun --map-by core --bind-to core -np 1 ./topology_parallel 128 256 20000 0.0001 789"
run_cmd "mpirun --map-by core --bind-to core -np 2 ./topology_parallel 128 256 20000 0.0001 789"
run_cmd "mpirun --map-by core --bind-to core -np 4 ./topology_parallel 128 256 20000 0.0001 789"
run_cmd "mpirun --map-by core --bind-to core -np 8 ./topology_parallel 128 256 20000 0.0001 789"
run_cmd "mpirun --map-by core --bind-to core -np 16 ./topology_parallel 128 256 20000 0.0001 789"
run_cmd "mpirun --map-by core --bind-to core -np 28 ./topology_parallel 128 256 20000 0.0001 789"

# RUN PER SEED 2024
echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 2024" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 128 256 20000 0.0001 2024"
run_cmd "mpirun --map-by core --bind-to core -np 1 ./topology_parallel 128 256 20000 0.0001 2024"
run_cmd "mpirun --map-by core --bind-to core -np 2 ./topology_parallel 128 256 20000 0.0001 2024"
run_cmd "mpirun --map-by core --bind-to core -np 4 ./topology_parallel 128 256 20000 0.0001 2024"
run_cmd "mpirun --map-by core --bind-to core -np 8 ./topology_parallel 128 256 20000 0.0001 2024"
run_cmd "mpirun --map-by core --bind-to core -np 16 ./topology_parallel 128 256 20000 0.0001 2024"
run_cmd "mpirun --map-by core --bind-to core -np 28 ./topology_parallel 128 256 20000 0.0001 2024"

# RUN PER SEED 4242
echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 4242" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 128 256 20000 0.0001 4242"
run_cmd "mpirun --map-by core --bind-to core -np 1 ./topology_parallel 128 256 20000 0.0001 4242"
run_cmd "mpirun --map-by core --bind-to core -np 2 ./topology_parallel 128 256 20000 0.0001 4242"
run_cmd "mpirun --map-by core --bind-to core -np 4 ./topology_parallel 128 256 20000 0.0001 4242"
run_cmd "mpirun --map-by core --bind-to core -np 8 ./topology_parallel 128 256 20000 0.0001 4242"
run_cmd "mpirun --map-by core --bind-to core -np 16 ./topology_parallel 128 256 20000 0.0001 4242"
run_cmd "mpirun --map-by core --bind-to core -np 28 ./topology_parallel 128 256 20000 0.0001 4242"

echo "========================================================" | tee -a "$OUT_FILE"
echo "ESECUZIONE COMPLETATA" | tee -a "$OUT_FILE"
