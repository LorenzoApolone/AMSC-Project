#!/usr/bin/env bash
set -euo pipefail

HOST="${DPSO_CLUSTER_HOST:-}"
REMOTE_SRC="${DPSO_REMOTE_SRC:-}"
BATTERY="all"
SEEDS_ARG=""
TAG=""
REMOTE_RESULTS_ROOT=""
LOCAL_PULL_DIR=""
PULL_RESULTS=0
SKIP_BUILD=0
DRY_RUN=0

usage() {
  printf '%s\n' \
    "Usage: dispatch_dpso_ssh.sh --host HOST --remote-src REMOTE_SRC [options]" \
    "" \
    "Options:" \
    "  --host HOST                    SSH host for the cluster/login node." \
    "  --remote-src PATH              Absolute remote path of the src/ directory." \
    "  --battery comparable|dpso|all  Select the benchmark battery." \
    "  --seeds 123,456,789            Override default seed list." \
    "  --tag NAME                     Custom run tag forwarded to the remote runner." \
    "  --remote-results-root PATH     Remote results root (default: REMOTE_SRC/benchmarks/dpso/results)." \
    "  --pull-results PATH            Copy remote results back into PATH after completion." \
    "  --skip-build                   Do not rebuild remotely before running." \
    "  --dry-run                      Print remote commands without executing them." \
    "  --help                         Show this message." \
    "" \
    "Environment variables:" \
    "  DPSO_CLUSTER_HOST              Default value for --host." \
    "  DPSO_REMOTE_SRC                Default value for --remote-src."
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

if [[ -z "$HOST" || -z "$REMOTE_SRC" ]]; then
  echo "Both --host and --remote-src are required." >&2
  usage >&2
  exit 1
fi

if [[ -z "$REMOTE_RESULTS_ROOT" ]]; then
  REMOTE_RESULTS_ROOT="$REMOTE_SRC/benchmarks/dpso/results"
fi

REMOTE_RUNNER="$REMOTE_SRC/benchmarks/dpso/run_benchmarks.sh"
REMOTE_CMD="cd $(printf '%q' "$REMOTE_SRC") && bash $(printf '%q' "$REMOTE_RUNNER") --battery $(printf '%q' "$BATTERY")"
if [[ -n "$SEEDS_ARG" ]]; then
  REMOTE_CMD+=" --seeds $(printf '%q' "$SEEDS_ARG")"
fi
if [[ -n "$TAG" ]]; then
  REMOTE_CMD+=" --tag $(printf '%q' "$TAG")"
fi
if [[ $SKIP_BUILD -eq 1 ]]; then
  REMOTE_CMD+=" --skip-build"
fi
if [[ $DRY_RUN -eq 1 ]]; then
  REMOTE_CMD+=" --dry-run"
fi

echo "[ssh] Launching remote benchmark run on $HOST"
ssh "$HOST" "$REMOTE_CMD"

if [[ $PULL_RESULTS -eq 1 ]]; then
  mkdir -p "$LOCAL_PULL_DIR"
  echo "[ssh] Pulling results into $LOCAL_PULL_DIR"
  rsync -az "$HOST:$REMOTE_RESULTS_ROOT/" "$LOCAL_PULL_DIR/"
fi
