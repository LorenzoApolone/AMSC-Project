#!/usr/bin/env bash
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
MATRIX_SCRIPT="$SCRIPT_DIR/benchmark_matrix.py"
RESULTS_ROOT="$SCRIPT_DIR/results"
BATTERY="all"
SEEDS_ARG=""
RUN_TAG="$(date +%Y%m%d_%H%M%S)"
DRY_RUN=0
SKIP_BUILD=0
FAILURES=0

MPIEXEC_BIN="${DPSO_MPIEXEC_BIN:-mpirun}"
MPIEXEC_ARGS="${DPSO_MPIEXEC_ARGS:-}"
TEST_PARALLEL_BINARY="${DPSO_TEST_BINARY:-$SRC_DIR/main_dpso}"
TEST_SERIAL_BINARY="${DPSO_SERIAL_BINARY:-$SRC_DIR/main_dpso_serial}"

timestamp() {
  date "+%Y-%m-%d %H:%M:%S"
}

usage() {
  cat <<'EOF'
Usage: run_dpso_benchmarks.sh [options]

Options:
  --battery NAME                  Select benchmark battery (all, comparable, dpso, dpso_weak_particles, dpso_weak_dimension, dpso_weak_local, dpso_sweep_dim).
  --seeds 123,456,789             Override default seed list.
  --results-root PATH             Directory where raw results are saved.
  --tag NAME                      Custom subdirectory name under results/root.
  --skip-build                    Do not rebuild dpso/main_dpso_serial and dpso/main_dpso.
  --dry-run                       Print commands without executing them.
  --help                          Show this message.

Environment variables:
  DPSO_MPIEXEC_BIN                MPI launcher binary (default: mpirun).
  DPSO_MPIEXEC_ARGS               Extra MPI launcher arguments.
  DPSO_TEST_BINARY                Override path to main_parallel.
  DPSO_SERIAL_BINARY              Override path to test_dpso.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --battery)
      BATTERY="$2"
      shift 2
      ;;
    --seeds)
      SEEDS_ARG="$2"
      shift 2
      ;;
    --results-root)
      RESULTS_ROOT="$2"
      shift 2
      ;;
    --tag)
      RUN_TAG="$2"
      shift 2
      ;;
    --skip-build)
      SKIP_BUILD=1
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

RAW_ROOT="$RESULTS_ROOT/raw/$RUN_TAG"
mkdir -p "$RAW_ROOT"

if [[ ! -x "$MATRIX_SCRIPT" ]]; then
  chmod +x "$MATRIX_SCRIPT"
fi

if [[ $SKIP_BUILD -eq 0 ]]; then
  echo "[build] Rebuilding DPSO serial and parallel test binaries..."
  make -C "$SRC_DIR" main_dpso main_dpso_serial -j2
fi

if [[ ! -x "$TEST_PARALLEL_BINARY" ]]; then
  echo "DPSO parallel test binary not found or not executable: $TEST_PARALLEL_BINARY" >&2
  exit 1
fi

if [[ ! -x "$TEST_SERIAL_BINARY" ]]; then
  echo "DPSO serial test binary not found or not executable: $TEST_SERIAL_BINARY" >&2
  exit 1
fi

MATRIX_CMD=(python3 "$MATRIX_SCRIPT" --battery "$BATTERY" --format tsv --no-header)
if [[ -n "$SEEDS_ARG" ]]; then
  MATRIX_CMD+=(--seeds "$SEEDS_ARG")
fi

CASE_LINES=()
while IFS= read -r line; do
  CASE_LINES+=("$line")
done < <("${MATRIX_CMD[@]}")
TOTAL_CASES="${#CASE_LINES[@]}"
if [[ $TOTAL_CASES -eq 0 ]]; then
  echo "No benchmark cases generated for battery '$BATTERY'." >&2
  exit 1
fi

RUN_MANIFEST="$RAW_ROOT/run_manifest.json"
python3 - <<PY > "$RUN_MANIFEST"
import json
manifest = {
    "battery": "$BATTERY",
    "results_root": "$RAW_ROOT",
    "mpiexec_bin": "$MPIEXEC_BIN",
    "mpiexec_args": "$MPIEXEC_ARGS",
    "parallel_test_binary": "$TEST_PARALLEL_BINARY",
    "serial_test_binary": "$TEST_SERIAL_BINARY",
    "seeds_override": "$SEEDS_ARG",
    "total_cases": int("$TOTAL_CASES"),
}
print(json.dumps(manifest, indent=2, ensure_ascii=True))
PY

if command -v stdbuf >/dev/null 2>&1; then
  RUN_PREFIX=(stdbuf -oL -eL)
else
  RUN_PREFIX=()
fi

echo "[run] Writing raw outputs to $RAW_ROOT"
echo "[run] Total cases: $TOTAL_CASES"

CASE_INDEX=0
for case_line in "${CASE_LINES[@]}"; do
  CASE_INDEX=$((CASE_INDEX + 1))
  IFS=$'\t' read -r battery suite family case_id execution_mode mpi_processes dim total_particles max_iters seed w c1 c2 regrouping_period sub_swarm_size hmcr par_min par_max description <<< "$case_line"

  if [[ "$execution_mode" == "serial" ]]; then
    CASE_LABEL="serial_dim${dim}"
    OUT_DIR="$RAW_ROOT/$battery/$suite/$family/$CASE_LABEL/seed_${seed}"
    CMD=(
      "$TEST_SERIAL_BINARY" "$dim" "$total_particles" "$max_iters" "1e-6" "$OUT_DIR/params.txt" "$seed"
    )
  else
    CASE_LABEL="np_${mpi_processes}_dim${dim}"
    OUT_DIR="$RAW_ROOT/$battery/$suite/$family/$CASE_LABEL/seed_${seed}"
    CMD=("$MPIEXEC_BIN")
    if [[ -n "$MPIEXEC_ARGS" ]]; then
      read -r -a EXTRA_MPI_ARGS <<< "$MPIEXEC_ARGS"
      CMD+=("${EXTRA_MPI_ARGS[@]}")
    fi
    CMD+=(
      -np "$mpi_processes"
      "$TEST_PARALLEL_BINARY"
      "$dim"
      "$total_particles"
      "$max_iters"
      "1e-6"
      "$OUT_DIR/params.txt"
      "$seed"
    )
  fi

  mkdir -p "$OUT_DIR"
  
  cat <<EOF > "$OUT_DIR/params.txt"
[GLOBAL]
w $w
c1 $c1
c2 $c2
regrouping_period $regrouping_period
sub_swarm_size $sub_swarm_size
hmcr $hmcr
par_min $par_min
par_max $par_max
EOF

  printf '%q ' "${CMD[@]}" > "$OUT_DIR/command.sh"
  printf '
' >> "$OUT_DIR/command.sh"

  if [[ $DRY_RUN -eq 1 ]]; then
    printf '[dry-run %d/%d] %s
' "$CASE_INDEX" "$TOTAL_CASES" "$case_id"
    cat "$OUT_DIR/command.sh"
    printf '
'
    continue
  fi

  printf '
[%s] [run %d/%d] %s
' "$(timestamp)" "$CASE_INDEX" "$TOTAL_CASES" "$case_id"
  printf '[%s]          mode=%s | suite=%s | family=%s | np=%s | dim=%s | tp=%s | max_iters=%s | seed=%s
'     "$(timestamp)" "$execution_mode" "$suite" "$family" "$mpi_processes" "$dim" "$total_particles" "$max_iters" "$seed"
  printf '[%s]          command: %s
' "$(timestamp)" "$(cat "$OUT_DIR/command.sh")"

  START_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  START_NS="$(date +%s%N)"
  set +e
  "${RUN_PREFIX[@]}" "${CMD[@]}"     > >(tee "$OUT_DIR/stdout.log")     2> >(tee "$OUT_DIR/stderr.log" >&2)
  RETURN_CODE=$?
  set -e
  END_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  END_NS="$(date +%s%N)"

  WALL_TIME_S="$(python3 - <<PY
start_ns = int("$START_NS")
end_ns = int("$END_NS")
print(f"{(end_ns - start_ns) / 1_000_000_000.0:.9f}")
PY
)"

  python3 - <<PY > "$OUT_DIR/metadata.json"
import json
metadata = {
    "battery": "$battery",
    "suite": "$suite",
    "family": "$family",
    "case_id": "$case_id",
    "execution_mode": "$execution_mode",
    "mpi_processes": int("$mpi_processes"),
    "dim": int("$dim"),
    "total_particles": int("$total_particles"),
    "max_iters": int("$max_iters"),
    "seed": int("$seed"),
    "w": float("$w"),
    "c1": float("$c1"),
    "c2": float("$c2"),
    "regrouping_period": int("$regrouping_period"),
    "sub_swarm_size": int("$sub_swarm_size"),
    "hmcr": float("$hmcr"),
    "par_min": float("$par_min"),
    "par_max": float("$par_max"),
    "description": "$description",
    "command": open("$OUT_DIR/command.sh", "r", encoding="utf-8").read().strip(),
    "stdout_log": "$OUT_DIR/stdout.log",
    "stderr_log": "$OUT_DIR/stderr.log",
    "start_utc": "$START_UTC",
    "end_utc": "$END_UTC",
    "suite_wall_time_s": float("$WALL_TIME_S"),
    "return_code": int("$RETURN_CODE"),
}
print(json.dumps(metadata, indent=2, ensure_ascii=True))
PY

  printf '[%s] [done %d/%d] %s | rc=%s | wall=%ss
'     "$(timestamp)" "$CASE_INDEX" "$TOTAL_CASES" "$case_id" "$RETURN_CODE" "$WALL_TIME_S"

  if [[ $RETURN_CODE -ne 0 ]]; then
    echo "[warn] Command failed for $case_id (seed=$seed, rc=$RETURN_CODE)" >&2
    FAILURES=$((FAILURES + 1))
  fi
done

if [[ $DRY_RUN -eq 1 ]]; then
  echo "[done] Dry run completed."
  exit 0
fi

if [[ $FAILURES -ne 0 ]]; then
  echo "[done] Completed with $FAILURES failing runs." >&2
  exit 1
fi

echo "[done] All runs completed successfully."
