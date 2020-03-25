# Notes for installing WGA_GPU locally on a desktop

This worked on my Ubuntu 18.04 desktop with a 4G GeForce GTX 1050 Ti.  It may be useful for others wishing to try the software on their laptop or desktop (as opposed to an AWS G3 or P3 instance). 

## Install CUDA-9 and Nvidia drivers

This was the simplest way I could find to get CUDA working ([from here](https://askubuntu.com/a/1036265)).

The steps using `cuda-repo-ubuntu1804-10-2-local-10.2.89-440.33.01_1.0-1_amd64.deb` from `install.sh` work fine on AWS instances but for some reason didn't work on my system.

```
sudo add-apt-repository ppa:graphics-drivers/ppa
sudo apt update
sudo ubuntu-drivers autoinstall
```

**Reboot**

```
sudo apt install nvidia-cuda-toolkit
nvcc --version
```

## Verfiy that CUDA is working

Building and running this small program ([from here](https://askubuntu.com/a/1215237)) is a basic installation test

```
#include <cassert>

#define N 3

__global__ void inc(int *a) {
    int i = blockIdx.x;
    if (i<N) {
        a[i]++;
    }
}

int main() {
    int ha[N], *da;
    cudaMalloc((void **)&da, N*sizeof(int));
    for (int i = 0; i<N; ++i) {
        ha[i] = i;
    }
    cudaMemcpy(da, ha, N*sizeof(int), cudaMemcpyHostToDevice);
    inc<<<N, 1>>>(da);
    cudaMemcpy(ha, da, N*sizeof(int), cudaMemcpyDeviceToHost);
    for (int i = 0; i < N; ++i) {
        assert(ha[i] == i + 1);
    }
    cudaFree(da);
    return 0;
}
```

Save as `test.cu` and compile and run:
```
nvcc test.cu -o test.out && ./test.out
```
Everything works if `./test.out` runs without printing an error. 

## Install remaining dependencies and build WGA_GPU

Using the `-c` option skips CUDA installation, as it was already done above.  Note, this script still requires `sudo` to install most packages at the system level

```
./install.sh -c
```

Whe using CUDA 9 as installed manually above, I got the following error when building wga in the `install.sh` script:
```
/WGA_GPU/build$ make
[ 11%] Building CXX object CMakeFiles/wga.dir/DRAM.cpp.o
In file included from /home/hickey/dev/WGA_GPU/DRAM.cpp:5:0:
/home/hickey/dev/WGA_GPU/tbb/include/tbb/tbb.h:21:154: note: #pragma message: TBB Warning: tbb.h contains deprecated functionality. For details, please see Deprecated Features appendix in the TBB reference manual.
 .h contains deprecated functionality. For details, please see Deprecated Features appendix in the TBB reference manual.")
                                                                                                                         ^
[ 22%] Building CUDA object CMakeFiles/wga.dir/GPU.cu.o
nvcc fatal   : redefinition of argument 'std'
CMakeFiles/wga.dir/build.make:86: recipe for target 'CMakeFiles/wga.dir/GPU.cu.o' failed
make[2]: *** [CMakeFiles/wga.dir/GPU.cu.o] Error 1
CMakeFiles/Makefile2:67: recipe for target 'CMakeFiles/wga.dir/all' failed
make[1]: *** [CMakeFiles/wga.dir/all] Error 2
Makefile:83: recipe for target 'all' failed
make: *** [all] Error 2

```
Removing the occurrences of `-std=c++11` from CMakeLists.txt resolved it:
```
sed -i 's/-std=c++11//g' CMakeLists.txt
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DTBB_ROOT=${PWD}/../tbb ..
make -j $(nproc)
```

As a developer, I find it easier to work with the local `wga` binary and scripts:

```
export PATH=$(pwd)/build:$(pwd)/bin:$PATH
```

## Run a tiny test

```
gzip -dc data/cb4.fa.gz | head -20000 > data/cb4_20k.fa
gzip -dc data/ce11.fa.gz | head -20000 > data/ce11_20k.fa
run_wga_gpu data/cb4_20k.fa data/ce11_20k.fa  > out.maf
```

Output:

```
Splitting reference chromosome
Converting chromosome wise fasta to 2bit format
Splitting query chromosome
Converting chromosome wise fasta to 2bit format
Using 8 threads
Using 1 GPU(s)

Reading query file ...

Reading target file ...

Start alignment ...

Sending reference cb4.chr1 ...

Sending query chr1 with buffer 0 ...

Starting query chr1 with buffer 0 ...
Chromosome chr1 interval 1/1 (0:999931) with buffer 0

Time elapsed (complete pipeline): 12 sec 
```


