#!/bin/bash

#Definisci il percorso della sottocartella
base_dir="./photon_dist/pinhole"
mkdir -p "$base_dir"  # Crea la cartella base se non esiste

# Crea una cartella con data e ora all'interno della sottocartella
timestamp=$(date +"%Y%m%d_%H%M%S")
output_dir="$base_dir/simulations_$timestamp"
mkdir "$output_dir"

# Ciclo che esegue la simulazione 100 volte
for i in {0..499}
do

  # Imposta il seed uguale all'indice 'i', o un altro calcolo se preferisci
  seed=$((4000 + i))

  # Definisce il nome del file di output
  output_file="$output_dir/hole500um_${seed}.raw" 
  log_file="$output_dir/log_hole500um_${seed}.log"
  
  # Comando per lanciare la simulazione
  ./flash init_vis.mac $seed $output_file | tee "$log_file"

  # Facoltativo: stampa per vedere l'iterazione in corso
    echo "Eseguito run $i con seed $seed, file di output $output_file e log $log_file"

done