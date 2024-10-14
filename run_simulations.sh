#!/bin/bash

#Definisci il percorso della sottocartella
base_dir="./photon_dist/telecentric"
mkdir -p "$base_dir"  # Crea la cartella base se non esiste

# Crea una cartella con data e ora all'interno della sottocartella
timestamp=$(date +"%Y%m%d_%H%M%S")
output_dir="$base_dir/simulations_$timestamp"
mkdir "$output_dir"

# Ciclo che esegue la simulazione 100 volte
for i in {0..30}
do

  # Imposta il seed uguale all'indice 'i', o un altro calcolo se preferisci
  seed=$((4000 + i))

  # Definisce il nome del file di output
  output_file="$output_dir/my_outputfile_${i}.raw" 
  
  # Comando per lanciare la simulazione
  #./flash init_vis.mac $seed $output_file >> "$output_dir/log_simulation.log" 2>&1
  ./flash init_vis.mac $seed $output_file | tee -a "$output_dir/log_simulation.log"

  # Facoltativo: stampa per vedere l'iterazione in corso
  echo "Eseguito run $i con seed $seed e file di output $output_file"
done