FROM vibsinglecellnf/scanpy:1.8.1

ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="git" && \
    apt-get update && \
    apt-get install -y --no-install-recommends apt-utils debconf locales && dpkg-reconfigure locales && \
    apt-get install -y --no-install-recommends $BUILDPKGS

RUN git clone https://github.com/tabdelaal/SpaGE.git \
    && cd SpaGE \
    && git checkout 7533cbf2275c3049561e8a17b9f7866e0e324743

ENV PYTHONPATH $PYTHONPATH:${PWD}/SpaGE/
