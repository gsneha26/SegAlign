#!/bin/bash -i 

CURR=$PWD

#sudo apt update && sudo apt dist-upgrade
#sudo apt-get install \
#    build-essential \
#    zlib1g \
#    clang \
#    zlib1g-dev \
#    libboost-dev \
#    libboost-all-dev \
#    libboost-thread-dev \
#    cmake \
#    parallel
#
##wget https://github.com/Kitware/CMake/releases/download/v3.16.3/cmake-3.16.3-Linux-x86_64.sh
##chmod +x cmake-3.16.3-Linux-x86_64.sh
##sudo sh cmake-3.16.3-Linux-x86_64.sh
##echo "export PATH="$CURR/cmake-3.16.3-Linux-x86_64/bin/:\$PATH"" >> ~/.bashrc
##rm cmake-3.16.3-Linux-x86_64.sh
##
##source ~/.bashrc
#
#git clone --recursive https://github.com/microsoft/bond.git
#curl -sSL https://get.haskellstack.org/ | sh
#cd bond
#mkdir build
#cd build
#cmake -DBOND_ENABLE_GRPC=FALSE ..
#make -j
#sudo make install
#cd $CURR
#sudo rm -rf bond

#wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/cuda-ubuntu1604.pin
#sudo mv cuda-ubuntu1604.pin /etc/apt/preferences.d/cuda-repository-pin-600
#wget http://developer.download.nvidia.com/compute/cuda/10.2/Prod/local_installers/cuda-repo-ubuntu1604-10-2-local-10.2.89-440.33.01_1.0-1_amd64.deb
#sudo dpkg -i cuda-repo-ubuntu1604-10-2-local-10.2.89-440.33.01_1.0-1_amd64.deb
#sudo apt-key add /var/cuda-repo-10-2-local-10.2.89-440.33.01/7fa2af80.pub
#sudo apt-get update
#sudo apt-get -y install cuda
#rm cuda-repo-ubuntu1604-10-2-local-10.2.89-440.33.01_1.0-1_amd64.deb

#wget http://www.bx.psu.edu/~rsharris/lastz/lastz-1.04.03.tar.gz
#gunzip lastz-1.04.03.tar.gz
#tar -xvf lastz-1.04.03.tar 
#cd $CURR/lastz-distrib-1.04.03/src
#make -j
#cd $CURR
#cp $CURR/lastz-distrib-1.04.03/src/lastz $CURR/bin/
#rm -rf lastz-distrib-1.04.03 lastz-1.04.03.tar
#
#cd $CURR/bin/
#wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/faToTwoBit
#chmod +x faToTwoBit
#wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/faSplit
#chmod +x faSplit
#wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/twoBitToFa
#chmod +x twoBitToFa
#
#echo "export PATH="$CURR/bin/:\$PATH"" >> ~/.bashrc
#
#source ~/.bashrc

cd $CURR
git clone https://github.com/01org/tbb
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DTBB_ROOT=${PWD}/../tbb ..
make -j
cp $CURR/build/wga $CURR/bin/
