#!/bin/bash
#PBS -N benchmark_16cores
#PBS -q cpu
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -l walltime=10:00:00
#PBS -j oe

# Caricamento ambiente
source /software/spack/share/spack/setup-env.sh
spack load gcc@15.2.0
spack load openmpi@5.0.8

# Spostamento nella cartella
cd "$PBS_O_WORKDIR"

# Definizione file di output
mkdir -p risultati_run
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUT_FILE="print_here.txt"

echo "========================================================" > "$OUT_FILE"
echo "ESECUZIONE: BENCHMARK 16 CORES" >> "$OUT_FILE"
echo "Avviata il: $(date)" >> "$OUT_FILE"
echo "========================================================" >> "$OUT_FILE"
echo "" >> "$OUT_FILE"

# Funzione interna per loggare
log_cmd() {
    echo "--------------------------------------------------------" >> "$OUT_FILE"
    echo "=> ESECUZIONE COMANDO: $1" >> "$OUT_FILE"
    eval "$1" >> "$OUT_FILE" 2>&1
    echo "" >> "$OUT_FILE"
    echo "" >> "$OUT_FILE"
}

# Esecuzione dei calcoli richiesti
log_cmd "mpirun --map-by core --bind-to core -np 16 ./topology_parallel 64 512 20000 0.0001 789"
log_cmd "mpirun --map-by core --bind-to core -np 16 ./topology_parallel 64 512 10000 0.0001 789"
log_cmd "mpirun -np 16 ./topology_parallel 64 512 10000 0.0001 456"
log_cmd "mpirun -np 16 ./topology_parallel 64 512 20000 0.0001 456"

echo "========================================================" >> "$OUT_FILE"
echo "ESECUZIONE COMPLETATA" >> "$OUT_FILE"

