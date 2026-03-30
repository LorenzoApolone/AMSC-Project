#!/bin/bash

cd "$(dirname "$0")/.." || exit 1

mkdir -p script_benchmark/risultati_run

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUT_FILE="script_benchmark/risultati_run/weak_dim64_${TIMESTAMP}.txt"

echo "========================================================" | tee -a "$OUT_FILE"
echo "ESECUZIONE: WEAK SCALING - PARTICELLE PER PROCESSO COSTANTI - DIM 64, PPS 32, K = NP" | tee -a "$OUT_FILE"
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

# RUN PER SEED 123
echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 123" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 64 32 10000 0.0001 123"
run_cmd "mpirun -np 1 ./topology_parallel 64 32 10000 0.0001 123"
run_cmd "mpirun -np 2 ./topology_parallel 64 64 10000 0.0001 123"
run_cmd "mpirun -np 4 ./topology_parallel 64 128 10000 0.0001 123"
run_cmd "mpirun -np 8 ./topology_parallel 64 256 10000 0.0001 123"
run_cmd "mpirun -np 16 ./topology_parallel 64 512 10000 0.0001 123"
run_cmd "mpirun -np 28 ./topology_parallel 64 896 10000 0.0001 123"

# RUN PER SEED 456
echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 456" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 64 32 10000 0.0001 456"
run_cmd "mpirun -np 1 ./topology_parallel 64 32 10000 0.0001 456"
run_cmd "mpirun -np 2 ./topology_parallel 64 64 10000 0.0001 456"
run_cmd "mpirun -np 4 ./topology_parallel 64 128 10000 0.0001 456"
run_cmd "mpirun -np 8 ./topology_parallel 64 256 10000 0.0001 456"
run_cmd "mpirun -np 16 ./topology_parallel 64 512 10000 0.0001 456"
run_cmd "mpirun -np 28 ./topology_parallel 64 896 10000 0.0001 456"

# RUN PER SEED 789
echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 789" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 64 32 10000 0.0001 789"
run_cmd "mpirun -np 1 ./topology_parallel 64 32 10000 0.0001 789"
run_cmd "mpirun -np 2 ./topology_parallel 64 64 10000 0.0001 789"
run_cmd "mpirun -np 4 ./topology_parallel 64 128 10000 0.0001 789"
run_cmd "mpirun -np 8 ./topology_parallel 64 256 10000 0.0001 789"
run_cmd "mpirun -np 16 ./topology_parallel 64 512 10000 0.0001 789"
run_cmd "mpirun -np 28 ./topology_parallel 64 896 10000 0.0001 789"

# RUN PER SEED 2024
echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 2024" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 64 32 10000 0.0001 2024"
run_cmd "mpirun -np 1 ./topology_parallel 64 32 10000 0.0001 2024"
run_cmd "mpirun -np 2 ./topology_parallel 64 64 10000 0.0001 2024"
run_cmd "mpirun -np 4 ./topology_parallel 64 128 10000 0.0001 2024"
run_cmd "mpirun -np 8 ./topology_parallel 64 256 10000 0.0001 2024"
run_cmd "mpirun -np 16 ./topology_parallel 64 512 10000 0.0001 2024"
run_cmd "mpirun -np 28 ./topology_parallel 64 896 10000 0.0001 2024"

# RUN PER SEED 4242
echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 4242" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 64 32 10000 0.0001 4242"
run_cmd "mpirun -np 1 ./topology_parallel 64 32 10000 0.0001 4242"
run_cmd "mpirun -np 2 ./topology_parallel 64 64 10000 0.0001 4242"
run_cmd "mpirun -np 4 ./topology_parallel 64 128 10000 0.0001 4242"
run_cmd "mpirun -np 8 ./topology_parallel 64 256 10000 0.0001 4242"
run_cmd "mpirun -np 16 ./topology_parallel 64 512 10000 0.0001 4242"
run_cmd "mpirun -np 28 ./topology_parallel 64 896 10000 0.0001 4242"

echo "========================================================" | tee -a "$OUT_FILE"
echo "ESECUZIONE COMPLETATA" | tee -a "$OUT_FILE"
