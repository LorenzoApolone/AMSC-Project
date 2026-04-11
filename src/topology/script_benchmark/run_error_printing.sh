#!/bin/bash

# Spostati nella cartella root degli eseguibili (src/topology)
cd "$(dirname "$0")/.." || exit 1

# Crea la cartella per i risultati se non esiste
mkdir -p script_benchmark/risultati_run

# Imposta il nome del file di output basato sul timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUT_FILE="script_benchmark/risultati_run/user_benchmark_${TIMESTAMP}.txt"

echo "========================================================" | tee -a "$OUT_FILE"
echo "  ESECUZIONE BENCHMARK RICHIESTO (Tutti gli output nel file)" | tee -a "$OUT_FILE"
echo "  Avviata il: $(date)" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
echo "" | tee -a "$OUT_FILE"

# Funzione per eseguire e loggare il comando
run_cmd() {
    local cmd="$*"
    echo "--------------------------------------------------------" | tee -a "$OUT_FILE"
    echo "=> ESECUZIONE COMANDO: $cmd" | tee -a "$OUT_FILE"
    # Esegue il comando, redirigendo sia stdout che stderr (2>&1)
    # in console (tee) e nel file
    eval "$cmd" 2>&1 | tee -a "$OUT_FILE"
    echo "" | tee -a "$OUT_FILE"
}

# Dimension 5
run_cmd "mpirun -np 8 ./topology_parallel 5 256 10000 0.0001 123"
run_cmd "mpirun -np 8 ./topology_parallel 5 256 10000 0.0001 456"
run_cmd "mpirun -np 8 ./topology_parallel 5 256 10000 0.0001 789"
run_cmd "mpirun -np 8 ./topology_parallel 5 256 10000 0.0001 2024"
run_cmd "mpirun -np 8 ./topology_parallel 5 256 10000 0.0001 4242"

# Dimension 32
run_cmd "mpirun -np 8 ./topology_parallel 32 256 10000 0.0001 123"
run_cmd "mpirun -np 8 ./topology_parallel 32 256 10000 0.0001 456"
run_cmd "mpirun -np 8 ./topology_parallel 32 256 10000 0.0001 789"
run_cmd "mpirun -np 8 ./topology_parallel 32 256 10000 0.0001 2024"
run_cmd "mpirun -np 8 ./topology_parallel 32 256 10000 0.0001 4242"

# Dimension 64
run_cmd "mpirun -np 8 ./topology_parallel 64 256 10000 0.0001 123"
run_cmd "mpirun -np 8 ./topology_parallel 64 256 10000 0.0001 456"
run_cmd "mpirun -np 8 ./topology_parallel 64 256 10000 0.0001 789"
run_cmd "mpirun -np 8 ./topology_parallel 64 256 10000 0.0001 2024"
run_cmd "mpirun -np 8 ./topology_parallel 64 256 10000 0.0001 4242"

# Dimension 128
run_cmd "mpirun -np 8 ./topology_parallel 128 256 20000 0.0001 123"
run_cmd "mpirun -np 8 ./topology_parallel 128 256 20000 0.0001 456"
run_cmd "mpirun -np 8 ./topology_parallel 128 256 20000 0.0001 789"
run_cmd "mpirun -np 8 ./topology_parallel 128 256 20000 0.0001 2024"
run_cmd "mpirun -np 8 ./topology_parallel 128 256 20000 0.0001 4242"

echo "========================================================" | tee -a "$OUT_FILE"
echo "  ESECUZIONE COMPLETATA" | tee -a "$OUT_FILE"
echo "  I risultati completi sono salvati in: $OUT_FILE" | tee -a "$OUT_FILE"
echo "========================================================" | tee -a "$OUT_FILE"
