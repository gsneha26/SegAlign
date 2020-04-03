#!/bin/bash
set -x

cd $HOME
sudo yum update -y
sudo yum -y groupinstall "Development Tools"
sudo yum -y install boost-devel
sudo yum -y install parallel
sudo yum -y install kernel-devel-4.14.154-99.181.amzn1.x86_64

aws s3 cp s3://wga-gpu/WGA_GPU.tar.gz .
tar -xvf WGA_GPU.tar.gz
rm WGA_GPU.tar.gz
cd WGA_GPU

CURR=$PWD
echo $CURR

#CMake
cd $CURR
aws s3 cp s3://wga-gpu/cmake-3.8.0-Linux-x86_64.tar.gz .
tar -xvf cmake-3.8.0-Linux-x86_64.tar.gz
echo "PATH=\"$CURR/cmake-3.8.0-Linux-x86_64/bin/:\$PATH\"" >> ~/.bashrc
sudo cp -r $CURR/cmake-3.8.0-Linux-x86_64/share/cmake-3.8 /usr/local/share/
rm cmake-3.8.0-Linux-x86_64.tar.gz

#NVIDIA driver and CUDA toolkit
cd $CURR
aws s3 cp s3://wga-gpu/cuda_9.2.148_396.37_linux.run .
chmod +x cuda_9.2.148_396.37_linux.run
sudo ./cuda_9.2.148_396.37_linux.run --silent --driver --toolkit --no-opengl-libs --kernel-source-path='/usr/src/kernels/4.14.154-99.181.amzn1.x86_64'
echo "PATH=\"/usr/local/cuda/bin/:\$PATH\"" >> ~/.bashrc
rm cuda_9.2.148_396.37_linux.run
source ~/.bashrc

#LASTZ
cd $CURR
aws s3 cp s3://wga-gpu/lastz-1.04.03.tar.gz .
tar -xvf lastz-1.04.03.tar.gz
cd $CURR/lastz-distrib-1.04.03/src
make -j $(nproc)
cp $CURR/lastz-distrib-1.04.03/src/lastz $CURR/bin/
rm -rf $CURR/lastz-distrib-1.04.03 $CURR/lastz-1.04.03.tar
rm $CURR/lastz-1.04.03.tar.gz

# kentUtils
cd $CURR/bin/
aws s3 cp s3://wga-gpu/faToTwoBit .
chmod +x faToTwoBit
aws s3 cp s3://wga-gpu/faSplit .
chmod +x faSplit
aws s3 cp s3://wga-gpu/twoBitToFa .
chmod +x twoBitToFa

cd $CURR
git clone https://github.com/01org/tbb
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DTBB_ROOT=${PWD}/../tbb ..
make -j $(nproc)
cp $CURR/build/wga $CURR/bin/
sudo cp $CURR/bin/* /usr/local/bin/

cd $CURR
mkdir data
aws s3 cp --recursive s3://wga-gpu/data data/
