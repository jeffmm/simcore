#!/bin/bash
[ -z "$1" ] && echo "No run_name supplied" && exit
[ ! -f "sjob_$1_j00.sh" ] && echo "No slurm job file with run_name $1 found" && exit
for i in sjob_"$1"_j*.sh; do
    echo "$i"
done
