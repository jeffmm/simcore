#!/bin/bash

params=""
for var in "$@"
do
  param="-v $PWD/$var:/app/run/$var"
  params="$params $param"
done

run_folder=`date +"%Y%m%d_%H%M%S"`
mkdir $run_folder

docker run -it --rm -v $PWD/$run_folder:/app/run/outputs $params jeffmm/simcore:latest
