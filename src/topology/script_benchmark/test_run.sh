#!/bin/bash

cd "$(dirname "$0")/.." || exit 1

mkdir -p script_benchmark/risultati_run
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUT_FILE="script_benchmark/risultati_run/strong_dim64_${TIMESTAMP}.txt"

echo "========================================================" | tee -a "$OUT_FILE"
echo "Script test" | tee -a "$OUT_FILE"
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
echo ">>> RUN PER SEED 123" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 4 60 1000 0.0001 123"
run_cmd "mpirun -np 1 ./topology_parallel 4 60 1000 0.0001 123"
run_cmd "mpirun -np 2 ./topology_parallel 4 60 1000 0.0001 123"
run_cmd "mpirun -np 4 ./topology_parallel 4 60 1000 0.0001 123"

echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 456" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 4 60 1000 0.0001 456"
run_cmd "mpirun -np 1 ./topology_parallel 4 60 1000 0.0001 456"
run_cmd "mpirun -np 2 ./topology_parallel 4 60 1000 0.0001 456"
run_cmd "mpirun -np 4 ./topology_parallel 4 60 1000 0.0001 456"

echo "========================================================" | tee -a "$OUT_FILE"
echo "ESECUZIONE COMPLETATA" | tee -a "$OUT_FILE"
