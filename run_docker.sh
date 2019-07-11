#!/bin/bash

params=""
flags=""
file=true
for var in "$@"
do
  case "$var" in
    *.yaml)
      if [ "$file" = true ]
      then
        param="-v $PWD/$var:/app/run/params.yaml"
        params="$params $param"
        let "file=!file"
      fi
      ;;
    -*)
      flags="$flags $var"
      ;;
    *)
      varname="${var##*/}"
      param="-v $PWD/$var:/app/run/$varname"
      params="$params $param"
      ;;
  esac
done
#echo $params
#echo $flags

run_folder=`date +"%Y%m%d_%H%M%S"`
mkdir $run_folder

docker run -it --rm -v $PWD/$run_folder:/app/run/outputs $params jeffmm/simcore:latest /bin/bash -c "./simcore $flags params.yaml; sleep 2; cp *.y* outputs; if ls *.[^y]* 1> /dev/null 2>&1; then mv *.[^y]* outputs; fi"
