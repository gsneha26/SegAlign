FROM nvidia/cuda:10.2-devel as builder

# system dependencies are installed by ./installUbuntu.sh below, but we need sudo first
RUN apt-get -qq -y update && \
    apt-get -qq -y upgrade && \
    apt-get -qq -y install \
    sudo

# Copy build tree into place
COPY . /WGA_GPU

# build and install everything but cuda using the install script
RUN cd /WGA_GPU && rm -rf build && ./installUbuntu.sh -c

# Create a thinner final Docker image with only runtime dependencies
FROM nvidia/cuda:10.2-runtime

# Install runtime dependencies
RUN apt-get -qq -y update && \
    apt-get -qq -y upgrade && \
    apt-get -qq -y install \
    libkrb5-3 \
    libk5crypto3 \
    libboost-dev \
    libboost-program-options-dev \
    zlib1g \
    parallel

# copy all the binaries
COPY --from=builder /usr/local/bin /usr/local/bin

# copy the tbb shared library
COPY --from=builder /WGA_GPU/build/tbb_cmake_build/tbb_cmake_build_subdir_release/lib* /usr/local/lib/

# add the library path
ENV LD_LIBRARY_PATH="/usr/local/lib/:${LD_LIBRARY_PATH}"

# UCSC convention is to work in /data
RUN mkdir -p /data
WORKDIR /data

