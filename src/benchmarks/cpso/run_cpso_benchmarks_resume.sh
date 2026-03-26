#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
MATRIX_SCRIPT="$SCRIPT_DIR/benchmark_matrix.py"
RESULTS_ROOT="$SCRIPT_DIR/results"
BATTERY="all"
SEEDS_ARG=""
RUN_TAG="$(date +%Y%m%d_%H%M%S)_resume"
DRY_RUN=0
SKIP_BUILD=0
FAILURES=0
START_INDEX=""
END_INDEX=""

MPIEXEC_BIN="${CPSO_MPIEXEC_BIN:-mpirun}"
MPIEXEC_ARGS="${CPSO_MPIEXEC_ARGS:-}"
TEST_PARALLEL_BINARY="${CPSO_TEST_BINARY:-$SRC_DIR/cpso/test_cpso_parallel}"
TEST_SERIAL_BINARY="${CPSO_SERIAL_BINARY:-$SRC_DIR/cpso/test_cpso}"

timestamp() {
  date "+%Y-%m-%d %H:%M:%S"
}

usage() {
  cat <<'EOF'
Usage: run_cpso_benchmarks_resume.sh [options]

Options:
  --start-index N                  First global case index to execute (1-based, required).
  --end-index N                    Last global case index to execute (inclusive).
  --battery comparable|cpso|all    Select the benchmark battery to run.
  --seeds 123,456,789              Override default seed list.
  --results-root PATH              Directory where raw results are saved.
  --tag NAME                       Custom subdirectory name under results/root.
  --skip-build                     Do not rebuild cpso/test_cpso and cpso/test_cpso_parallel.
  --dry-run                        Print commands without executing them.
  --help                           Show this message.

Environment variables:
  CPSO_MPIEXEC_BIN                 MPI launcher binary (default: mpirun).
  CPSO_MPIEXEC_ARGS                Extra MPI launcher arguments.
  CPSO_TEST_BINARY                 Override path to test_cpso_parallel.
  CPSO_SERIAL_BINARY               Override path to test_cpso.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start-index)
      START_INDEX="$2"
      shift 2
      ;;
    --end-index)
      END_INDEX="$2"
      shift 2
      ;;
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

if [[ -z "$START_INDEX" ]]; then
  echo "--start-index is required for resume runs." >&2
  usage >&2
  exit 1
fi

if ! [[ "$START_INDEX" =~ ^[0-9]+$ ]]; then
  echo "--start-index must be a positive integer." >&2
  exit 1
fi

if (( START_INDEX < 1 )); then
  echo "--start-index must be >= 1." >&2
  exit 1
fi

if [[ -n "$END_INDEX" ]]; then
  if ! [[ "$END_INDEX" =~ ^[0-9]+$ ]]; then
    echo "--end-index must be a positive integer." >&2
    exit 1
  fi
  if (( END_INDEX < START_INDEX )); then
    echo "--end-index must be >= --start-index." >&2
    exit 1
  fi
fi

RAW_ROOT="$RESULTS_ROOT/raw/$RUN_TAG"
mkdir -p "$RAW_ROOT"

if [[ ! -x "$MATRIX_SCRIPT" ]]; then
  chmod +x "$MATRIX_SCRIPT"
fi

if [[ $SKIP_BUILD -eq 0 ]]; then
  echo "[build] Rebuilding CPSO serial and parallel test binaries..."
  make -C "$SRC_DIR/cpso" -j2
fi

if [[ ! -x "$TEST_PARALLEL_BINARY" ]]; then
  echo "CPSO parallel test binary not found or not executable: $TEST_PARALLEL_BINARY" >&2
  exit 1
fi

if [[ ! -x "$TEST_SERIAL_BINARY" ]]; then
  echo "CPSO serial test binary not found or not executable: $TEST_SERIAL_BINARY" >&2
  exit 1
fi

MATRIX_CMD=(python3 "$MATRIX_SCRIPT" --battery "$BATTERY" --format tsv --no-header)
if [[ -n "$SEEDS_ARG" ]]; then
  MATRIX_CMD+=(--seeds "$SEEDS_ARG")
fi

mapfile -t CASE_LINES < <("${MATRIX_CMD[@]}")
TOTAL_CASES="${#CASE_LINES[@]}"
if [[ $TOTAL_CASES -eq 0 ]]; then
  echo "No benchmark cases generated for battery '$BATTERY'." >&2
  exit 1
fi

if (( START_INDEX > TOTAL_CASES )); then
  echo "--start-index $START_INDEX is beyond the available cases ($TOTAL_CASES)." >&2
  exit 1
fi

if [[ -z "$END_INDEX" ]] || (( END_INDEX > TOTAL_CASES )); then
  END_INDEX="$TOTAL_CASES"
fi

SELECTED_CASES=$((END_INDEX - START_INDEX + 1))

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
    "total_generated_cases": int("$TOTAL_CASES"),
    "start_index": int("$START_INDEX"),
    "end_index": int("$END_INDEX"),
    "selected_cases": int("$SELECTED_CASES"),
    "resume_mode": True,
}
print(json.dumps(manifest, indent=2, ensure_ascii=True))
PY

if command -v stdbuf >/dev/null 2>&1; then
  RUN_PREFIX=(stdbuf -oL -eL)
else
  RUN_PREFIX=()
fi

echo "[run] Writing raw outputs to $RAW_ROOT"
echo "[run] Resuming cases $START_INDEX-$END_INDEX out of $TOTAL_CASES"

CASE_INDEX=0
EXECUTED_INDEX=0
for case_line in "${CASE_LINES[@]}"; do
  CASE_INDEX=$((CASE_INDEX + 1))

  if (( CASE_INDEX < START_INDEX || CASE_INDEX > END_INDEX )); then
    continue
  fi

  EXECUTED_INDEX=$((EXECUTED_INDEX + 1))
  IFS=$'\t' read -r battery suite family case_id execution_mode mpi_processes dim k_subswarms particles_per_swarm max_iters shuffle_freq stagnation_patience seed description <<< "$case_line"

  if [[ "$execution_mode" == "serial" ]]; then
    CASE_LABEL="serial"
    OUT_DIR="$RAW_ROOT/$battery/$suite/$family/$CASE_LABEL/seed_${seed}"
    CMD=(
      "$TEST_SERIAL_BINARY"
      "$dim"
      "$k_subswarms"
      "$particles_per_swarm"
      "$max_iters"
      "$seed"
    )
  else
    CASE_LABEL="np_${mpi_processes}"
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
      "$k_subswarms"
      "$particles_per_swarm"
      "$max_iters"
      "$shuffle_freq"
      "$stagnation_patience"
      "$seed"
    )
  fi

  mkdir -p "$OUT_DIR"
  printf '%q ' "${CMD[@]}" > "$OUT_DIR/command.sh"
  printf '\n' >> "$OUT_DIR/command.sh"

  if [[ $DRY_RUN -eq 1 ]]; then
    printf '[dry-run %d/%d | global %d/%d] %s\n' "$EXECUTED_INDEX" "$SELECTED_CASES" "$CASE_INDEX" "$TOTAL_CASES" "$case_id"
    cat "$OUT_DIR/command.sh"
    printf '\n'
    continue
  fi

  printf '\n[%s] [run %d/%d | global %d/%d] %s\n' "$(timestamp)" "$EXECUTED_INDEX" "$SELECTED_CASES" "$CASE_INDEX" "$TOTAL_CASES" "$case_id"
  printf '[%s]          mode=%s | suite=%s | family=%s | np=%s | dim=%s | k=%s | pps=%s | max_iters=%s | seed=%s\n' \
    "$(timestamp)" "$execution_mode" "$suite" "$family" "$mpi_processes" "$dim" "$k_subswarms" "$particles_per_swarm" "$max_iters" "$seed"
  printf '[%s]          command: %s\n' "$(timestamp)" "$(cat "$OUT_DIR/command.sh")"

  START_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  START_NS="$(date +%s%N)"
  set +e
  "${RUN_PREFIX[@]}" "${CMD[@]}" \
    > >(tee "$OUT_DIR/stdout.log") \
    2> >(tee "$OUT_DIR/stderr.log" >&2)
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
    "k_subswarms": int("$k_subswarms"),
    "particles_per_swarm": int("$particles_per_swarm"),
    "max_iters": int("$max_iters"),
    "shuffle_freq": int("$shuffle_freq"),
    "stagnation_patience": int("$stagnation_patience"),
    "seed": int("$seed"),
    "description": "$description",
    "command": open("$OUT_DIR/command.sh", "r", encoding="utf-8").read().strip(),
    "stdout_log": "$OUT_DIR/stdout.log",
    "stderr_log": "$OUT_DIR/stderr.log",
    "start_utc": "$START_UTC",
    "end_utc": "$END_UTC",
    "suite_wall_time_s": float("$WALL_TIME_S"),
    "return_code": int("$RETURN_CODE"),
    "global_case_index": int("$CASE_INDEX"),
}
print(json.dumps(metadata, indent=2, ensure_ascii=True))
PY

  printf '[%s] [done %d/%d | global %d/%d] %s | rc=%s | wall=%ss\n' \
    "$(timestamp)" "$EXECUTED_INDEX" "$SELECTED_CASES" "$CASE_INDEX" "$TOTAL_CASES" "$case_id" "$RETURN_CODE" "$WALL_TIME_S"

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
  echo "[done] Resume run completed with $FAILURES failing runs." >&2
  exit 1
fi

echo "[done] Resume run completed successfully."