#!/bin/bash

# Si assicura di essere nella cartella dove risiede lo script (src/benchmarks/topology)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
cd "$SCRIPT_DIR" || exit 1

INPUT_FILE="topology_benchmark_runs.txt"
OUTPUT_FILE="result_pso.txt"

# Usiamo un percorso assoluto per il file di output così non lo perdiamo quando cambiamo cartella
ABS_OUTPUT_FILE="$SCRIPT_DIR/$OUTPUT_FILE"

# Pulisce/crea il file di output
> "$ABS_OUTPUT_FILE"

# Verifica che il file dei test esista in questa cartella
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Errore: file $INPUT_FILE non trovato in $SCRIPT_DIR."
    echo "Assicurati di avere il file di testo nella stessa cartella di questo script."
    exit 1
fi

echo "Inizio esecuzione dei comandi da $INPUT_FILE..."
echo "I risultati saranno salvati in $ABS_OUTPUT_FILE"
echo ""

# Poiché gli eseguibili vengono creati con `make` in src/topology/, 
# ci spostiamo temporaneamente in quella cartella per lanciare i file
cd ../../topology || { echo "Errore: cartella src/topology/ non trovata. Assicurati che lo script sia nella cartella giusta (src/benchmarks/topology/)"; exit 1; }

# Legge dal file txt originale
while IFS= read -r line || [[ -n "$line" ]]; do
    # Rimuovi eventuali spazi iniziali e finali
    cmd=$(echo "$line" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
    
    # Controlla se la riga è un comando (inizia con ./ o mpirun)
    if [[ "$cmd" == ./* ]] || [[ "$cmd" == mpirun* ]]; then
        
        # Corregge gli errori di battitura nel file txt
        # Rimuove l'eventuale "/topology" in eccesso dopo il ./ per lanciare gli eseguibili correttamente
        cmd=$(echo "$cmd" | sed 's|\./topology/|\./|g')
        
        # Stampa a schermo
        echo "start executing $cmd"
        
        # Stampa su file
        echo "start executing $cmd" >> "$ABS_OUTPUT_FILE"
        
        # Esegue il comando e salva l'output
        eval "$cmd" >> "$ABS_OUTPUT_FILE" 2>&1
    fi
done < "$SCRIPT_DIR/$INPUT_FILE"

echo ""
echo "Esecuzione completata! Controlla il file result_pso.txt nella cartella dei benchmark."
