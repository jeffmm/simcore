#! /bin/bash
#
# Run to start up docker container. Use argument 'build' to force rebuild of images.
# ./run_docker.sh build

# If 'build' argument was provided, rebuild images
if [ "$1" == "build" ]; then
  echo "Building docker image"
  docker build -t jeffmm/simcore docker
# Otherwise, just start up the containers
else
  echo "Running docker container"
  docker run --rm -itd -v ${PWD}:/mnt --name "simcore-executer" jeffmm/simcore bash
fi
