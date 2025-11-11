#!/bin/bash

q_ranges=("q1to10" "q10to100") 
bools=("false" "true")  # boolean values

# Set permissions once
mkdir -p log results
chmod 777 log results

for qrange in "${q_ranges[@]}"; do
  mkdir -p "./log/${qrange}"   
  mkdir -p "./results/${qrange}"
  
  for i in $(seq -w 0 9); do
    for b1 in "${bools[@]}"; do
      for b2 in "${bools[@]}"; do
        input_file="inputfiles/ep_10x100_${qrange}/split/part_${i}.list"
        output_file="./results/${qrange}/JetOut_ep_${qrange}_${i}_b1_${b1}_b2_${b2}.histos.root"
        log_file="./log/${qrange}/out_ep_${qrange}_${i}_b1_${b1}_b2_${b2}.log"

        echo "Running ROOT job for $qrange part $i with booleans $b1, $b2 ..."
        root -l -q "src/JetTreesRecluster.C(\"$input_file\",\"$output_file\",$b1,$b2)" &> "$log_file"
      done
    done
  done
done

