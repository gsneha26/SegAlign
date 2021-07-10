FROM nvidia/cuda:10.2-devel-ubuntu18.04
MAINTAINER Sneha D. Goenka, gsneha@stanford.edu

USER root
WORKDIR /home

RUN apt-get update && \
    apt-get -y install git cmake build-essential libboost-all-dev parallel zlib1g-dev wget && \ 
    apt-get clean && \
    apt-get purge

RUN git clone --recursive https://github.com/gsneha26/SegAlign.git 
WORKDIR SegAlign
ENV PROJECT_DIR=/home/SegAlign

RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/twoBitToFa && \
    chmod +x twoBitToFa && \
    mv twoBitToFa /usr/local/bin/ && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/faToTwoBit && \
    chmod +x faToTwoBit && \
    mv faToTwoBit /usr/local/bin/ && \
    cd ${PROJECT_DIR}/submodules/lastz/src && \
    make -j $(nproc) && \
    cp ${PROJECT_DIR}/submodules/lastz/src/lastz /usr/local/bin/ && \
    mkdir -p ${PROJECT_DIR}/build && \
    cd ${PROJECT_DIR}/build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DTBB_ROOT=${PROJECT_DIR}/submodules/TBB -DCMAKE_PREFIX_PATH=${PROJECT_DIR}/submodules/TBB/cmake ..  && \
    make -j $(nproc) && \
    cp ${PROJECT_DIR}/build/segalign* /usr/local/bin   && \
    cp ${PROJECT_DIR}/scripts/run_segalign* /usr/local/bin && \
    rm -rf ${PROJECT_DIR}/bin
