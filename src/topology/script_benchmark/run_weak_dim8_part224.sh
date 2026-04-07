#!/bin/bash

cd "$(dirname "$0")/.." || exit 1

mkdir -p script_benchmark/risultati_run

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUT_FILE="script_benchmark/risultati_run/weak_dim8_part224_${TIMESTAMP}.txt"

echo "========================================================" | tee -a "$OUT_FILE"
echo "ESECUZIONE: WEAK SCALING - DIMENSIONE PER PROCESSO COSTANTE - LOCAL DIM 8, PARTICELLE TOTALI 224, K = NP" | tee -a "$OUT_FILE"
echo "Avviata il: $(date)" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
echo "" | tee -a "$OUT_FILE"

# Funzione per eseguire e loggare il comando
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
run_cmd "./topology_serial 8 224 20000 0.0001 123"
run_cmd "mpirun --map-by core --bind-to core -np 1 ./topology_parallel 8 224 20000 0.0001 123"
run_cmd "mpirun --map-by core --bind-to core -np 2 ./topology_parallel 16 224 20000 0.0001 123"
run_cmd "mpirun --map-by core --bind-to core -np 4 ./topology_parallel 32 224 20000 0.0001 123"
run_cmd "mpirun --map-by core --bind-to core -np 8 ./topology_parallel 64 224 20000 0.0001 123"
run_cmd "mpirun --map-by core --bind-to core -np 16 ./topology_parallel 128 224 20000 0.0001 123"
run_cmd "mpirun --map-by core --bind-to core -np 28 ./topology_parallel 224 224 20000 0.0001 123"

# RUN PER SEED 456
echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 456" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 8 224 20000 0.0001 456"
run_cmd "mpirun --map-by core --bind-to core -np 1 ./topology_parallel 8 224 20000 0.0001 456"
run_cmd "mpirun --map-by core --bind-to core -np 2 ./topology_parallel 16 224 20000 0.0001 456"
run_cmd "mpirun --map-by core --bind-to core -np 4 ./topology_parallel 32 224 20000 0.0001 456"
run_cmd "mpirun --map-by core --bind-to core -np 8 ./topology_parallel 64 224 20000 0.0001 456"
run_cmd "mpirun --map-by core --bind-to core -np 16 ./topology_parallel 128 224 20000 0.0001 456"
run_cmd "mpirun --map-by core --bind-to core -np 28 ./topology_parallel 224 224 20000 0.0001 456"

# RUN PER SEED 789
echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 789" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 8 224 20000 0.0001 789"
run_cmd "mpirun --map-by core --bind-to core -np 1 ./topology_parallel 8 224 20000 0.0001 789"
run_cmd "mpirun --map-by core --bind-to core -np 2 ./topology_parallel 16 224 20000 0.0001 789"
run_cmd "mpirun --map-by core --bind-to core -np 4 ./topology_parallel 32 224 20000 0.0001 789"
run_cmd "mpirun --map-by core --bind-to core -np 8 ./topology_parallel 64 224 20000 0.0001 789"
run_cmd "mpirun --map-by core --bind-to core -np 16 ./topology_parallel 128 224 20000 0.0001 789"
run_cmd "mpirun --map-by core --bind-to core -np 28 ./topology_parallel 224 224 20000 0.0001 789"

# RUN PER SEED 2024
echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 2024" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 8 224 20000 0.0001 2024"
run_cmd "mpirun --map-by core --bind-to core -np 1 ./topology_parallel 8 224 20000 0.0001 2024"
run_cmd "mpirun --map-by core --bind-to core -np 2 ./topology_parallel 16 224 20000 0.0001 2024"
run_cmd "mpirun --map-by core --bind-to core -np 4 ./topology_parallel 32 224 20000 0.0001 2024"
run_cmd "mpirun --map-by core --bind-to core -np 8 ./topology_parallel 64 224 20000 0.0001 2024"
run_cmd "mpirun --map-by core --bind-to core -np 16 ./topology_parallel 128 224 20000 0.0001 2024"
run_cmd "mpirun --map-by core --bind-to core -np 28 ./topology_parallel 224 224 20000 0.0001 2024"

# RUN PER SEED 4242
echo "========================================================" | tee -a "$OUT_FILE"
echo ">>> RUN PER SEED 4242" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
run_cmd "./topology_serial 8 224 20000 0.0001 4242"
run_cmd "mpirun --map-by core --bind-to core -np 1 ./topology_parallel 8 224 20000 0.0001 4242"
run_cmd "mpirun --map-by core --bind-to core -np 2 ./topology_parallel 16 224 20000 0.0001 4242"
run_cmd "mpirun --map-by core --bind-to core -np 4 ./topology_parallel 32 224 20000 0.0001 4242"
run_cmd "mpirun --map-by core --bind-to core -np 8 ./topology_parallel 64 224 20000 0.0001 4242"
run_cmd "mpirun --map-by core --bind-to core -np 16 ./topology_parallel 128 224 20000 0.0001 4242"
run_cmd "mpirun --map-by core --bind-to core -np 28 ./topology_parallel 224 224 20000 0.0001 4242"

echo "========================================================" | tee -a "$OUT_FILE"
echo "ESECUZIONE COMPLETATA" | tee -a "$OUT_FILE"
