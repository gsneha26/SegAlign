#!/bin/bash -i 

usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [Options]\n"
    printf "Options:\n"
    printf "\t-c \tDo not install CUDA\n"
    exit 1
}

while getopts "hc" o; do
    case "${o}" in
        c)
            DONT_INSTALL_CUDA=1
            ;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

CURR=$PWD

set -x

# linux essentials
sudo apt update 
sudo apt-get install \
	build-essential \
	zlib1g-dev \
	libboost-all-dev \
	cmake \
	parallel


# NVIDIA CUDA
if [ -z ${DONT_INSTALL_CUDA} ]; then
cd $CURR
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin
sudo mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget http://developer.download.nvidia.com/compute/cuda/10.2/Prod/local_installers/cuda-repo-ubuntu1804-10-2-local-10.2.89-440.33.01_1.0-1_amd64.deb
sudo dpkg -i cuda-repo-ubuntu1804-10-2-local-10.2.89-440.33.01_1.0-1_amd64.deb
sudo apt-key add /var/cuda-repo-10-2-local-10.2.89-440.33.01/7fa2af80.pub
sudo apt-get update
sudo apt-get -y install cuda
rm cuda-repo-ubuntu1804-10-2-local-10.2.89-440.33.01_1.0-1_amd64.deb
export PATH=/usr/local/cuda-10.2/bin/:$PATH
fi

# LASTZ
cd $CURR
wget http://www.bx.psu.edu/~rsharris/lastz/lastz-1.04.03.tar.gz
gunzip lastz-1.04.03.tar.gz
tar -xvf lastz-1.04.03.tar 
cd $CURR/lastz-distrib-1.04.03/src
make -j $(nproc)
cp $CURR/lastz-distrib-1.04.03/src/lastz $CURR/bin/
rm -rf $CURR/lastz-distrib-1.04.03 $CURR/lastz-1.04.03.tar

# kentUtils
cd $CURR/bin/
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/faToTwoBit
chmod +x faToTwoBit
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/faSplit
chmod +x faSplit
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/twoBitToFa
chmod +x twoBitToFa
sudo cp $CURR/bin/* /usr/local/bin/

# WGA_GPU
cd $CURR
git clone https://github.com/01org/tbb
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DTBB_ROOT=${PWD}/../tbb ..
make -j $(nproc)
cp $CURR/build/wga $CURR/bin/		
sudo cp $CURR/bin/* /usr/local/bin/
