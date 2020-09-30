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
echo $CURR

set -x

# linux essentials
sudo apt update 
for pkg in cmake build-essential libboost-all-dev parallel zlib1g-dev wget git;
do
  REQUIRED_PKG=$pkg
  PKG_OK=$(dpkg-query -W --showformat='${Status}\n' $REQUIRED_PKG|grep "install ok installed")
  echo Checking for $REQUIRED_PKG: $PKG_OK
  if [ "" = "$PKG_OK" ]; then
    sudo apt-get --yes install $REQUIRED_PKG
  fi
done

echo "${DONT_INSTALL_CUDA}"

# NVIDIA CUDA
cd $CURR
if [ -z ${DONT_INSTALL_CUDA} ]; then
  if [ -z "$(command -v nvcc)" ]; then
  currentver="$(cat /usr/local/cuda/version.txt | awk '{print $3}')"
  requiredver="10.2.89"
    if [ "$(printf '%s\n' "$requiredver" "$currentver" | uniq | wc -l )" = 2 ]; then
      if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1 )" = "$currentver" ]; then
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
    fi
  fi
fi

# LASTZ
cd $CURR
if [ -z "$(command -v lastz)" ]; then
  mkdir bin
  wget http://www.bx.psu.edu/~rsharris/lastz/lastz-1.04.03.tar.gz
  gunzip lastz-1.04.03.tar.gz
  tar -xvf lastz-1.04.03.tar 
  cd $CURR/lastz-distrib-1.04.03/src
  make -j $(nproc)
  sudo cp $CURR/lastz-distrib-1.04.03/src/lastz /usr/local/bin/
  rm -rf $CURR/lastz-distrib-1.04.03 $CURR/lastz-1.04.03.tar
fi

# kentUtils - faToTwoBit
cd $CURR
if [ -z "$(command -v faToTwoBit)" ]; then
  cd $CURR/bin/
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/faToTwoBit
  chmod +x faToTwoBit
  sudo mv faToTwoBit /usr/local/bin/
fi

# kentUtils - twoBitToFa
cd $CURR
if [ -z "$(command -v twoBitToFa)" ]; then
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/twoBitToFa
  chmod +x twoBitToFa
  sudo mv twoBitToFa /usr/local/bin/
fi

# TBB
cd $CURR
wget https://github.com/oneapi-src/oneTBB/releases/download/2019_U9/tbb2019_20191006oss_lin.tgz
tar -xvf tbb2019_20191006oss_lin.tgz
rm tbb2019_20191006oss_lin.tgz

# SegAlign 
cd $CURR
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DTBB_ROOT=${PWD}/../tbb2019_20191006oss -DCMAKE_PREFIX_PATH=${PWD}/../tbb2019_20191006oss/cmake ..
make -j $(nproc)
sudo cp $CURR/build/segalign* /usr/local/bin	
sudo cp $CURR/scripts/run_segalign* /usr/local/bin
rm -rf $CURR/bin
