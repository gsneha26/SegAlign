[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://opensource.org/licenses/MIT

[![License][license-badge]][license-link]
[![Build Status](https://travis-ci.com/gsneha26/SegAlign.svg?branch=master)](https://travis-ci.com/gsneha26/SegAlign)
[![Published in SC20](https://img.shields.io/badge/published%20in-SC20-blue.svg)](https://doi.ieeecomputersociety.org/10.1109/SC41405.2020.00043)

<img src="logo.png" width="300">

A Scalable GPU System for Pairwise Whole Genome Alignments based on LASTZ's seed-filter-extend paradigm.

## Table of Contents

- [Overview](#overview)
- [Dependencies](#dependencies)
- [How to run SegAlign](#run)
    - [Running a test](#test)
- [How to run SegAlign repeat masker](#run_rm)
    - [Running a test](#test_rm)
- [Running Docker Image](#docker)
    - [Running segalign](#d_segalign)
    - [Running segalign_repeat_masker](#d_segalign_rm)
- [Citing SegAlign](#cite_segalign)

## <a name="overview"></a> Overview

The system has been tested on all the AWS G3 and P3 GPU instances with AMI Ubuntu Server 18.04 LTS (HVM), SSD Volume Type (ami-0fc20dd1da406780b (64-bit x86))

* Clone SegAlign repository (https://github.com/gsneha26/SegAlign)

```
git clone https://github.com/gsneha26/SegAlign.git
export PROJECT_DIR=$PWD/SegAlign
```

## <a name="dependencies"></a> Dependencies
The following dependencies are required by SegAlign:
  * NVIDIA CUDA 10.2 toolkit
  * CMake 3.8
  * Intel TBB library
  * libboost-all-dev
  * parallel
  * zlib
  * LASTZ 1.04.15
  * faToTwoBit, twoBitToFa (from kentUtils)

The dependencies can be installed with the given script as follows, which might take a while (only installs the dependencies not present already). This script requires sudo to install most packages at the system level. Using the `-c` option skips CUDA installation [the CUDA toolkit binaries should be in `$PATH` for SegAlign]. 

```
cd $PROJECT_DIR
./scripts/installUbuntu.sh
```

## <a name="run"></a> How to run SegAlign
* Run SegAlign

```
run_segalign target query [options]
```

* For a list of options 

```
run_segalign --help
```

### <a name="test"></a> Running a test

```
cd $PROJECT_DIR
mkdir test
cd test
wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.2bit
wget https://hgdownload-test.gi.ucsc.edu/goldenPath/cb4/bigZips/cb4.2bit 
twoBitToFa ce11.2bit ce11.fa
twoBitToFa cb4.2bit cb4.fa
run_segalign ce11.fa cb4.fa --output=ce11.cb4.maf
```

## <a name="run_rm"></a> How to run SegAlign repeat masker
* Run SegAlign repeat masker

```
run_segalign_repeat_masker sequence [options]
```

* For a list of options 

```
run_segalign_repeat_masker --help
```

### <a name="test_rm"></a> Running a test

```
cd $PROJECT_DIR
mkdir test_rm
cd test_rm
wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.2bit
twoBitToFa ce11.2bit ce11.fa
run_segalign_repeat_masker ce11.fa --output=ce11.seg
```

## <a name="docker"></a> Running Docker Image
### <a name="d_segalign"></a> Running segalign 
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.2bit
wget https://hgdownload-test.gi.ucsc.edu/goldenPath/cb4/bigZips/cb4.2bit 
sudo docker run -v $(pwd):/data -it gsneha/segalign:v0.1.2 \
                           twoBitToFa \
                           /data/ce11.2bit \
                           /data/ce11.fa
sudo docker run -v $(pwd):/data -it gsneha/segalign:v0.1.2 \
                           twoBitToFa \
                           /data/cb4.2bit \
                           /data/cb4.fa
sudo docker run --ipc=host --gpus all -v $(pwd):/data -it gsneha/segalign:v0.1.2 \
                           run_segalign \
                           /data/ce11.fa \
                           /data/cb4.fa \
                           --output=/data/ce11.cb4.maf
```

### <a name="d_segalign_rm"></a> Running segalign_repeat_masker

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.2bit
sudo docker run -v $(pwd):/data -it gsneha/segalign:v0.1.2 \
                           twoBitToFa \
                           /data/ce11.2bit \
                           /data/ce11.fa
sudo docker run --ipc=host --gpus all -v $(pwd):/data -it gsneha/segalign:v0.1.2 \
                           run_segalign_repeat_masker \
                           /data/ce11.fa \
                           --output=/data/ce11.seg
```

## <a name="cite_segalign"></a> Citing SegAlign

S. Goenka, Y. Turakhia, B. Paten and M. Horowitz,  "SegAlign: A Scalable GPU-Based Whole Genome Aligner," in 2020 SC20: International Conference for High Performance Computing, Networking, Storage and Analysis (SC), Atlanta, GA, US, 2020 pp. 540-552. doi: 10.1109/SC41405.2020.00043
