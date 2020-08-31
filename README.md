[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://opensource.org/licenses/MIT

[![License][license-badge]][license-link]
[![Build Status](https://travis-ci.com/gsneha26/SegAlign.svg?branch=master)](https://travis-ci.com/gsneha26/SegAlign)

# SegAlign 

A Scalable GPU System for Pairwise Whole Genome Alignments based on LASTZ's seed-filter-extend paradigm.

## Table of Contents

- [Overview](#overview)
- [Dependencies](#dependencies)
- [How to run SegAlign](#run)
    - [Running a test](#test)
- [How to run SegAlign repeat masker](#run_rm)
    - [Running a test](#test_rm)
- [Citing SegAlign](#cite_segalign)

## <a name="overview"></a> Overview

The system has been tested on all the AWS G3 and P3 GPU instances with AMI Ubuntu Server 18.04 LTS (HVM), SSD Volume Type (ami-0fc20dd1da406780b (64-bit x86))

* Clone SegAlign repository (https://github.com/gsneha26/SegAlign)

```
    $ git clone https://github.com/gsneha26/SegAlign.git
    $ export PROJECT_DIR=$PWD/SegAlign
```

## <a name="dependencies"></a> Dependencies
The following dependencies are required by SegAlign:
  * NVIDIA CUDA 10.2 toolkit
  * CMake 3.8
  * Intel TBB library
  * libboost-all-dev
  * parallel
  * zlib
  * LASTZ 1.04.03
  * faToTwoBit, twoBitToFa (from kentUtils)

The dependencies can be installed with the given script as follows, which might take a while (only installs the dependencies not present already). This script requires sudo to install most packages at the system level. Using the `-c` option skips CUDA installation. 

```
    $ cd $PROJECT_DIR
    $ ./scripts/installUbuntu.sh
```

## <a name="run"></a> How to run SegAlign
* Run SegAlign

```
    $ run_segalign target query [options]
```

* For a list of options 

```
    $ run_segalign --help
```

### <a name="test"></a> Running a test

```
    $ cd $PROJECT_DIR
    $ mkdir test
    $ cd test
    $ wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.2bit
    $ wget https://hgdownload-test.gi.ucsc.edu/goldenPath/cb4/bigZips/cb4.2bit 
    $ twoBitToFa ce11.2bit ce11.fa
    $ twoBitToFa cb4.2bit cb4.fa
    $ run_segalign ce11.fa cb4.fa --output=ce11.cb4.maf
```

## <a name="run_rm"></a> How to run SegAlign repeat masker
* Run SegAlign repeat masker

```
    $ run_segalign_repeat_masker sequence [options]
```

* For a list of options 

```
    $ run_segalign_repeat_masker --help
```

### <a name="test_rm"></a> Running a test

```
    $ cd $PROJECT_DIR
    $ mkdir test_rm
    $ cd test_rm
    $ wget https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.2bit
    $ twoBitToFa ce11.2bit ce11.fa
    $ run_segalign_repeat_masker ce11.fa --output=ce11.seg
```

## <a name="cite_segalign"></a> Citing SegAlign

Goenka, Sneha D., Yatish Turakhia, Benedict Paten and Mark Horowitz. "SegAlign: A Scalable GPU-Based Whole Genome Aligner." In SC20: International Conference for High Performance Computing, Networking, Storage and Analysis. 
