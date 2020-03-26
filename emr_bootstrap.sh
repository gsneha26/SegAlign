CURR=$PWD

sudo yum -y groupinstall 'Development Tools'
sudo yum -y install boost-devel \
       parallel

#CMake
wget https://github.com/Kitware/CMake/releases/download/v3.8.0/cmake-3.8.0-Linux-x86_64.tar.gz
tar -xvf cmake-3.8.0-Linux-x86_64.tar.gz
sudo cp $CURR/cmake-3.8.0-Linux-x86_64/bin/* /usr/local/bin/
sudo cp -r $CURR/cmake-3.8.0-Linux-x86_64/share/cmake-3.8 /usr/local/share/

#NVIDIA driver and CUDA toolkit
cd $CURR
wget http://developer.download.nvidia.com/compute/cuda/10.2/Prod/local_installers/cuda_10.2.89_440.33.01_linux.run
sudo ./cuda_10.2.89_440.33.01_linux.run --silent --driver --toolkit
export PATH=/usr/local/cuda-10.2/bin/:$PATH
rm cuda_10.2.89_440.33.01_linux.run

#LASTZ
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

cd $CURR
git clone https://github.com/01org/tbb
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DTBB_ROOT=${PWD}/../tbb ..
make -j $(nproc)
cp $CURR/build/wga $CURR/bin/
sudo cp $CURR/bin/* /usr/local/bin/
