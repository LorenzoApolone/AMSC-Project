#!/usr/bin/env bash

# this files dele the data created by run_sweep.sh
set -euo pipefail

OUTDIR="results"



if [[ ! -d "$OUTDIR" ]]; then
  echo "Nessuna cartella '$OUTDIR' trovata. Nulla da pulire."
  exit 0
fi

echo "Cartella trovata: $OUTDIR"
echo "Contenuto attuale:"
ls -lah "$OUTDIR" || true
echo

read -r -p "Vuoi cancellare TUTTO dentro '$OUTDIR'? (yes/no): " ans
if [[ "$ans" != "yes" ]]; then
  echo "Annullato."
  exit 1
fi

rm -rf "$OUTDIR"/runs.csv "$OUTDIR"/logs
mkdir -p "$OUTDIR"/logs

echo "Pulizia completata."
echo "Ricreato: $OUTDIR/logs"