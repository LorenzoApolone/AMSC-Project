#!/usr/bin/env bash
set -euo pipefail

HOST="${CPSO_CLUSTER_HOST:-}"
REMOTE_SRC="${CPSO_REMOTE_SRC:-}"
BATTERY="all"
SEEDS_ARG=""
TAG=""
REMOTE_RESULTS_ROOT=""
LOCAL_PULL_DIR=""
PULL_RESULTS=0

usage() {
  printf '%s\n' \
    "Usage: dispatch_cpso_ssh.sh --host HOST --remote-src REMOTE_SRC [options]" \
    "" \
    "Options:" \
    "  --host HOST                    SSH host for the cluster/login node." \
    "  --remote-src PATH              Absolute remote path of the src/ directory." \
    "  --battery comparable|cpso|appendix|all  Select the benchmark battery." \
    "  --seeds 123,456,789            Override default seed list." \
    "  --tag NAME                     Custom run tag forwarded to the remote runner." \
    "  --remote-results-root PATH     Remote results root (default: REMOTE_SRC/benchmarks/cpso/results)." \
    "  --pull-results PATH            Copy remote results back into PATH after completion." \
    "  --help                         Show this message." \
    "" \
    "Environment variables:" \
    "  CPSO_CLUSTER_HOST              Default value for --host." \
    "  CPSO_REMOTE_SRC                Default value for --remote-src."
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --host)
      HOST="$2"
      shift 2
      ;;
    --remote-src)
      REMOTE_SRC="$2"
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
    --tag)
      TAG="$2"
      shift 2
      ;;
    --remote-results-root)
      REMOTE_RESULTS_ROOT="$2"
      shift 2
      ;;
    --pull-results)
      LOCAL_PULL_DIR="$2"
      PULL_RESULTS=1
      shift 2
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

if [[ -z "$HOST" || -z "$REMOTE_SRC" ]]; then
  echo "Both --host and --remote-src are required." >&2
  usage >&2
  exit 1
fi

if [[ -z "$REMOTE_RESULTS_ROOT" ]]; then
  REMOTE_RESULTS_ROOT="$REMOTE_SRC/benchmarks/cpso/results"
fi

REMOTE_RUNNER="$REMOTE_SRC/benchmarks/cpso/run_cpso_benchmarks.sh"
REMOTE_CMD="cd $(printf '%q' "$REMOTE_SRC") && bash $(printf '%q' "$REMOTE_RUNNER") --battery $(printf '%q' "$BATTERY") --results-root $(printf '%q' "$REMOTE_RESULTS_ROOT")"
if [[ -n "$SEEDS_ARG" ]]; then
  REMOTE_CMD+=" --seeds $(printf '%q' "$SEEDS_ARG")"
fi
if [[ -n "$TAG" ]]; then
  REMOTE_CMD+=" --tag $(printf '%q' "$TAG")"
fi

echo "[ssh] Launching remote benchmark run on $HOST"
ssh "$HOST" "$REMOTE_CMD"

if [[ $PULL_RESULTS -eq 1 ]]; then
  mkdir -p "$LOCAL_PULL_DIR"
  echo "[ssh] Pulling results into $LOCAL_PULL_DIR"
  rsync -az "$HOST:$REMOTE_RESULTS_ROOT/" "$LOCAL_PULL_DIR/"
fi

