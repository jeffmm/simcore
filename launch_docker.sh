#! /bin/bash
#
# Run to launch a docker container named "simcore_dev".
# Provide the argument 'build' to force rebuilding the project image:
# ./run_docker.sh build

# If 'build' argument was provided, rebuild image
if [ "$1" == "build" ]; then
  echo "Building docker image"
  docker build -t jeffmm/simcore docker
# Otherwise, just start up the containers
else
  echo "Running docker container"
  docker run --rm -itd -v "${PWD}":/mnt --name "simcore_dev" jeffmm/simcore bash
fi
