cd $HOME
#sudo apt update && sudo apt dist-upgrade
#sudo apt install build-essential
#sudo apt install zlib1g-dev
#sudo apt install zlib1g
#sudo apt-get install libboost-dev
#
#mkdir bin 
#wget https://github.com/Kitware/CMake/releases/download/v3.16.3/cmake-3.16.3-Linux-x86_64.sh
#chmod +x cmake-3.16.3-Linux-x86_64.sh
#sudo sh cmake-3.16.3-Linux-x86_64.sh
#echo "export PATH="$HOME/cmake-3.16.3-Linux-x86_64/bin/:\$PATH"" >> ~/.bashrc
#rm cmake-3.16.3-Linux-x86_64.sh
#
#wget http://developer.download.nvidia.com/compute/cuda/10.2/Prod/local_installers/cuda_10.2.89_440.33.01_linux.run
#sudo sh cuda_10.2.89_440.33.01_linux.run
#echo "export PATH="/usr/local/cuda/bin/:\$PATH"" >> ~/.bashrc
#rm cuda_10.2.89_440.33.01_linux.run
#
#wget http://www.bx.psu.edu/~rsharris/lastz/lastz-1.04.03.tar.gz
#gunzip lastz-1.04.03.tar.gz
#tar -xvf lastz-1.04.03.tar 
#cd lastz-distrib-1.04.03/src
#make -j
#cp $HOME/lastz-distrib-1.04.03/src/lastz $HOME/bin/
#rm -rf lastz-distrib-1.04.03 lastz-1.04.03.tar

cd $HOME/bin/
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/faToTwoBit
chmod +x faToTwoBit
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/faSplit
chmod +x faSplit
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/twoBitToFa
chmod +x twoBitToFa

echo "export PATH="\$HOME/bin/:\$PATH"" >> ~/.bashrc

cd $HOME/WGA_GPU
git clone https://github.com/01org/tbb

source ~/.bashrc
#Download bond which Yatish sent via Gmail
