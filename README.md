# WGA_GPU 

sudo apt update && sudo apt dist-upgrade
sudo apt install build-essential
sudo apt install zlib1g-dev
sudo apt install zlib1g

git clone https://github.com/01org/tbb 

```
    $ cd build 
    $ cmake -DCMAKE_BUILD_TYPE=Release -DTBB_ROOT=${PWD}/../tbb .. 
    $ make
    $ ./wga
```
