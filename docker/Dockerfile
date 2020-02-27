FROM jeffmm/simcore_base:latest
WORKDIR /build
RUN git clone --recursive --single-branch --branch jmm/directory_tree_refactor https://github.com/jeffmm/simcore.git . &&\
    ./install.sh -otI
RUN rm -rf /build
WORKDIR /mnt
ENTRYPOINT []
