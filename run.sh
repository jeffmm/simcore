#! /bin/bash
#
# Run to start up docker containers. Use argument 'build' to force rebuild of images.
# ./run.sh build

# Check that environment variables are initialized in .env
if [ -f ".env" ]; then
    echo "Found environment variables in .env"
else
  echo "Cannot find .env file"
  exit 1
fi

# If 'build' argument was provided, rebuild images
if [ "$1" == "build" ]; then
  echo "Building images with docker-compose"
  docker-compose -f docker/docker-compose.yml --project-name mdsim up -d --build
# Otherwise, just start up the containers
else
  echo "Running containers with docker-compose"
  docker-compose -f docker/docker-compose.yml --project-name mdsim up -d
fi
