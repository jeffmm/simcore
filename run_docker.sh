#! /bin/bash
#
# Run to start up docker containers. Use argument 'build' to force rebuild of images.
# ./run.sh build

# If 'build' argument was provided, rebuild images
if [ "$1" == "build" ]; then
  echo "Building images with docker-compose"
  docker-compose -f docker/docker-compose.yml --project-name simcore up -d --build
# Otherwise, just start up the containers
else
  echo "Running containers with docker-compose"
  docker-compose -f docker/docker-compose.yml --project-name simcore up -d
fi
