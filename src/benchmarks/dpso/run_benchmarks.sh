#!/usr/bin/env bash
# =============================================================================
# run_benchmarks.sh
# Script di benchmark per DPSO (parallelo MPI) e DPSO seriale.
# Analogo a run_cpso_benchmarks.sh ma adattato alla firma di main_dpso:
#   mpirun -np <p> main_dpso <dim> <total_particles> <max_iter> [conv_tol]
#   main_dpso_serial     <dim> <total_particles> <max_iter> [conv_tol]
#
# Utilizzo:
#   bash run_benchmarks.sh [--battery <battery>] [--tag <tag>]
#
# --battery:  strong | weak_ppc | weak_local | dim_sweep | serial | all
#             (default: all)
# --tag:      etichetta per la cartella output (default: test)
# =============================================================================
set -euo pipefail
# ---------------------------------------------------------------------------
# Parametri di default (adattabili)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(cd "${SCRIPT_DIR}/../../" && pwd)"
DPSO_EXE="${SRC_DIR}/main_dpso"
DPSO_SERIAL_EXE="${SRC_DIR}/main_dpso_serial"
BATTERY="all"
TAG="test"
SEEDS=(123 456 789 2024 4242)
NPS=(1 2 4 8 16 32)
CONV_TOL="1e-4"
# Cartella output
OUTPUT_BASE="${SCRIPT_DIR}/results"
# ---------------------------------------------------------------------------
# Parsing argomenti
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --battery) BATTERY="$2"; shift 2 ;;
    --tag)     TAG="$2";     shift 2 ;;
    *) echo "Argomento non riconosciuto: $1"; exit 1 ;;
  esac
done
OUTPUT_DIR="${OUTPUT_BASE}/${TAG}"
mkdir -p "${OUTPUT_DIR}"
# ---------------------------------------------------------------------------
# Utilità: lancio singola run e salvataggio output
# ---------------------------------------------------------------------------
run_dpso() {
  local label="$1"   # nome famiglia/tag della run
  local np="$2"      # numero processi MPI (0 = seriale)
  local dim="$3"
  local total_particles="$4"
  local max_iter="$5"
  local seed="$6"    # usato solo per il naming del file output (DPSO non ha seed CLI)
  local conv_tol="${7:-${CONV_TOL}}"
  local out_dir="${OUTPUT_DIR}/${label}"
  mkdir -p "${out_dir}"
  local out_file="${out_dir}/np${np}_dim${dim}_p${total_particles}_iter${max_iter}_seed${seed}.txt"
  if [[ "${np}" -eq 0 ]]; then
    # --- esecuzione seriale ---
    echo "[SERIAL] ${label} | dim=${dim} p=${total_particles} iter=${max_iter} seed(label)=${seed}"
    "${DPSO_SERIAL_EXE}" "${dim}" "${total_particles}" "${max_iter}" "${conv_tol}" \
      > "${out_file}" 2>&1
  else
    # --- esecuzione MPI ---
    echo "[MPI np=${np}] ${label} | dim=${dim} p=${total_particles} iter=${max_iter} seed(label)=${seed}"
    mpirun -np "${np}" "${DPSO_EXE}" "${dim}" "${total_particles}" "${max_iter}" "${conv_tol}" \
      > "${out_file}" 2>&1
  fi
}
# ---------------------------------------------------------------------------
# BATTERY 1 – STRONG SCALING
# Problema totale fisso: dim e total_particles costanti, np varia.
# Equivalente al benchmark CPSO strong_scaling_fixed_algorithm.
#
# Configurazioni:
#   dim=32,  total_particles=256, max_iter=10000
#   dim=64,  total_particles=256, max_iter=10000
#   dim=128, total_particles=256, max_iter=20000
# ---------------------------------------------------------------------------
battery_strong() {
  echo "=== STRONG SCALING ==="
  for seed in "${SEEDS[@]}"; do
    for np in "${NPS[@]}"; do
      run_dpso "strong_dim32_p256"  "${np}"  32  256  10000 "${seed}"
      run_dpso "strong_dim64_p256"  "${np}"  64  256  10000 "${seed}"
      run_dpso "strong_dim128_p256" "${np}" 128  256  20000 "${seed}"
    done
  done
}
# ---------------------------------------------------------------------------
# BATTERY 2 – WEAK SCALING (particelle per processo costanti)
# total_particles = PPS * np, dove PPS=32.
# dim fissa = 64, np varia.
# Equivalente al CPSO weak_particles_per_process_constant.
# ---------------------------------------------------------------------------
battery_weak_ppc() {
  echo "=== WEAK SCALING – PARTICELLE PER PROCESSO COSTANTI (PPS=32) ==="
  local PPS=32
  for seed in "${SEEDS[@]}"; do
    for np in "${NPS[@]}"; do
      local total=$(( PPS * np ))
      run_dpso "weak_ppc_dim64_pps32" "${np}" 64 "${total}" 10000 "${seed}"
    done
  done
}
# ---------------------------------------------------------------------------
# BATTERY 3 – WEAK SCALING (carico locale costante, dim/np cresce)
# Ogni processo gestisce una porzione fissa: local_dim=4, local_p=32.
# dim = 4*np, total_particles = 32*np.
# Equivalente al CPSO weak_local_load_constant.
# ---------------------------------------------------------------------------
battery_weak_local() {
  echo "=== WEAK SCALING – CARICO LOCALE COSTANTE (local_dim=4, local_p=32) ==="
  local LOCAL_DIM=4
  local LOCAL_P=32
  for seed in "${SEEDS[@]}"; do
    for np in "${NPS[@]}"; do
      local dim=$(( LOCAL_DIM * np ))
      local total=$(( LOCAL_P * np ))
      local iters=10000
      [[ "${np}" -ge 16 ]] && iters=20000
      run_dpso "weak_local_dim${LOCAL_DIM}_p${LOCAL_P}" "${np}" "${dim}" "${total}" "${iters}" "${seed}"
    done
  done
}
# ---------------------------------------------------------------------------
# BATTERY 4 – DIMENSION SWEEP
# np=8, total_particles=128, dim varia: 16, 32, 64, 128, 256.
# Equivalente al CPSO dimension_sweep_fixed_np8.
# ---------------------------------------------------------------------------
battery_dim_sweep() {
  echo "=== DIMENSION SWEEP (np=8, total_particles=128) ==="
  local NP_FIXED=8
  local P_FIXED=128
  declare -A ITERS_MAP=([16]=10000 [32]=10000 [64]=15000 [128]=20000 [256]=40000)
  for seed in "${SEEDS[@]}"; do
    for dim in 16 32 64 128 256; do
      local iters="${ITERS_MAP[${dim}]}"
      run_dpso "dim_sweep_np8_p128" "${NP_FIXED}" "${dim}" "${P_FIXED}" "${iters}" "${seed}"
    done
  done
}
# ---------------------------------------------------------------------------
# BATTERY 5 – SERIAL BASELINE
# Esegue main_dpso_serial con le stesse configurazioni dello strong scaling
# e del dimension sweep per confronto diretto.
# ---------------------------------------------------------------------------
battery_serial() {
  echo "=== SERIAL BASELINE ==="
  for seed in "${SEEDS[@]}"; do
    # Strong scaling configs (np=0 → serial)
    run_dpso "serial_strong_dim32_p256"  0  32  256  10000 "${seed}"
    run_dpso "serial_strong_dim64_p256"  0  64  256  10000 "${seed}"
    run_dpso "serial_strong_dim128_p256" 0 128  256  20000 "${seed}"
    # Dim sweep configs
    for dim_item in "16:10000" "32:10000" "64:15000" "128:20000" "256:40000"; do
      local dim="${dim_item%%:*}"
      local iters="${dim_item##*:}"
      run_dpso "serial_dim_sweep_p128" 0 "${dim}" 128 "${iters}" "${seed}"
    done
  done
}
# ---------------------------------------------------------------------------
# Dispatching battery
# ---------------------------------------------------------------------------
echo "========================================"
echo " DPSO Benchmark Runner"
echo " Battery : ${BATTERY}"
echo " Tag     : ${TAG}"
echo " Output  : ${OUTPUT_DIR}"
echo "========================================"
case "${BATTERY}" in
  strong)    battery_strong    ;;
  weak_ppc)  battery_weak_ppc  ;;
  weak_local) battery_weak_local ;;
  dim_sweep) battery_dim_sweep ;;
  serial)    battery_serial    ;;
  all)
    battery_strong
    battery_weak_ppc
    battery_weak_local
    battery_dim_sweep
    battery_serial
    ;;
  *)
    echo "Battery non riconosciuta: ${BATTERY}"
    echo "Valori validi: strong | weak_ppc | weak_local | dim_sweep | serial | all"
    exit 1
    ;;
esac
echo ""
echo "=== TUTTE LE RUN COMPLETATE ==="
echo "Output salvato in: ${OUTPUT_DIR}"
