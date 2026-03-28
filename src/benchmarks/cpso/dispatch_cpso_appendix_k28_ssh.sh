#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HOST="${CPSO_CLUSTER_HOST:-}"
REMOTE_SRC="${CPSO_REMOTE_SRC:-}"
SEEDS_ARG=""
TAG="cpso_appendix_k28"
REMOTE_RESULTS_ROOT=""
LOCAL_PULL_DIR=""
PULL_RESULTS=0

usage() {
  printf '%s\n' \
    "Usage: dispatch_cpso_appendix_k28_ssh.sh --host HOST --remote-src REMOTE_SRC [options]" \
    "" \
    "Options:" \
    "  --host HOST                SSH host for the cluster/login node." \
    "  --remote-src PATH          Absolute remote path of the src/ directory." \
    "  --seeds 123,456,789        Override default seed list." \
    "  --tag NAME                 Custom run tag (default: cpso_appendix_k28)." \
    "  --remote-results-root PATH Remote results root (default: REMOTE_SRC/benchmarks/cpso/results)." \
    "  --pull-results PATH        Copy remote results back into PATH after completion." \
    "  --help                     Show this message." \
    "" \
    "This launcher always runs the CPSO appendix battery with fixed k=28."
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

CMD=(bash "$SCRIPT_DIR/dispatch_cpso_ssh.sh" --host "$HOST" --remote-src "$REMOTE_SRC" --battery appendix --tag "$TAG")
if [[ -n "$SEEDS_ARG" ]]; then
  CMD+=(--seeds "$SEEDS_ARG")
fi
if [[ -n "$REMOTE_RESULTS_ROOT" ]]; then
  CMD+=(--remote-results-root "$REMOTE_RESULTS_ROOT")
fi
if [[ $PULL_RESULTS -eq 1 ]]; then
  CMD+=(--pull-results "$LOCAL_PULL_DIR")
fi

"${CMD[@]}"
