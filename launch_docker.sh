#! /bin/bash
#
# Run to launch a docker container named "simcore_latest".
# Provide the flag '-b' to force rebuilding the project image:
# >> ./run_docker.sh -b
# Provide the flag '-x' to launch/build experimental image from jmm/experimental branch
# >> ./run_docker.sh -bx

show_help() {
    echo "Without additional options, $0 launches a docker container named simcore_latest to run in the background"
    echo "USAGE:"
    echo "  $0 [-hbx]"
    echo "OPTIONS:"
    echo "  -h      show this menu"
    echo "  -b      force rebuild of simcore Docker image"
    echo "  -x      build/launch experimental version of simcore Docker image"
}

# Reset in case getopts has been used previously in the shell.
OPTIND=1         
experimental=false
build=false
while getopts "h?bx" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    b)  
        build=true
        ;;
    x)  
        experimental=true
        ;;
    esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

if $build; then
    if $experimental; then
        echo "Building experimental docker image"
        docker build --no-cache -f docker/Dockerfile_experimental -t jeffmm/simcore:experimental docker
    else
        echo "Building docker image"
        docker build --no-cache -t jeffmm/simcore:latest docker
    fi
elif $experimental; then
    echo "Launching experimental simcore docker container"
    docker run --rm -itd -v "${PWD}":/mnt --name "simcore_experimental" jeffmm/simcore:experimental bash
# Otherwise, just start up the containers
else
  echo "Launching simcore docker container"
  docker run --rm -itd -v "${PWD}":/mnt --name "simcore_latest" jeffmm/simcore bash
fi
