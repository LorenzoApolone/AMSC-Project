#!/usr/bin/env bash
set -euo pipefail

# =========================
# Configurazione (modifica qui)
# =========================
EXE="./a"
OUTDIR="results"
LOGDIR="$OUTDIR/logs"
CSV="$OUTDIR/runs.csv"

# Sweep parameters (modifica facilmente queste liste)
NPS=(1 2 4)
DIMS=(3 10)
NPOINTS=(100)
MAXITERS=(100 1000)
DELTAS=(0.001)

# Ripetizioni per ogni combinazione
REPS=1

# =========================
# Setup cartelle + header CSV
# =========================
mkdir -p "$OUTDIR" "$LOGDIR"

# Header nuovo (coerente): time_total = tempo totale; time_comm = tempo comunicazione (allgatherv, se presente)
if [[ ! -f "$CSV" ]]; then
  echo "np,dim,n_points,max_iter,delta,rep,version,time_total,time_comm,conv,total" > "$CSV"
fi

# Estrae il valore di key=... da una riga "RESULT,..."
get_kv() {
  local line="$1"
  local key="$2"
  echo "$line" | sed -n "s/.*${key}=\\([^,]*\\).*/\\1/p"
}

run_idx=0

# =========================
# Sweep
# =========================
for np in "${NPS[@]}"; do
  for dim in "${DIMS[@]}"; do
    for npt in "${NPOINTS[@]}"; do
      for mi in "${MAXITERS[@]}"; do
        for dx in "${DELTAS[@]}"; do
          for rep in $(seq 1 "$REPS"); do

            run_idx=$((run_idx+1))
            run_id="run${run_idx}_np${np}_d${dim}_n${npt}_mi${mi}_dx${dx}_rep${rep}"
            logfile="$LOGDIR/${run_id}.out"

            echo ">>> $run_id"
            mpirun -np "$np" "$EXE" "$dim" "$npt" "$mi" "$dx" | tee "$logfile" >/dev/null

            # Legge SOLO le righe strutturate
            while IFS= read -r line; do
              version="$(get_kv "$line" "version")"
              t_time="$(get_kv "$line" "time")"
              t1="$(get_kv "$line" "time1")"
              t2="$(get_kv "$line" "time2")"
              conv="$(get_kv "$line" "conv")"
              total="$(get_kv "$line" "total")"

              # Normalizzazione:
              # - classic: ha "time" (tempo totale)
              # - altre topologie: hanno "time1" (comm) e "time2" (totale)
              if [[ -n "${t_time:-}" ]]; then
                time_total="$t_time"
                time_comm=""
              else
                time_total="${t2:-}"
                time_comm="${t1:-}"
              fi

              echo "$np,$dim,$npt,$mi,$dx,$rep,$version,${time_total:-},${time_comm:-},${conv:-},${total:-}" >> "$CSV"
            done < <(grep '^RESULT,version=' "$logfile" || true)

          done
        done
      done
    done
  done
done

echo "Fatto. CSV in: $CSV"