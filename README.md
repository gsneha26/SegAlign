# SegAlign 

A Scalable GPU System for Pairwise Whole Genome Alignments

## Table of Contents

- [Overview](#overview)
- [Dependenciesi](#dependencies)
- [How to run SegAlign](#run)
    - [Running a test](#test)

## <a name="overview"></a> Overview

The system has been tested on all the AWS G3 and P3 GPU instances with AMI Ubuntu Server 18.04 LTS (HVM), SSD Volume Type (ami-0fc20dd1da406780b (64-bit x86))

* Clone SegAlign repository (https://github.com/gsneha26/SegAlign)

```
    $ git clone https://github.com/gsneha26/SegAlign.git
    $ export PROJECT_DIR=$PWD/SegAlign
```

## <a name="dependencies"></a> Dependencies
The following dependencies are required by SegAlign:
  * zlib1g-dev
  * libboost-all-dev
  * CMake 3.8
  * parallel
  * CUDA 10.2 toolkit
  * LASTZ 1.04.03
  * faToTwoBit (from kentUtils)
  * Intel TBB library

The dependencies can be installed with the given script as follows (this might take a while): 

```
    $ cd $PROJECT_DIR
    $ ./installUbuntu.sh
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
    $ gunzip $PROJECT_DIR/data/ce11.fa.gz
    $ gunzip $PROJECT_DIR/data/cb4.fa.gz
    $ run_segalign $PROJECT_DIR/data/ce11.fa $PROJECT_DIR/data/cb4.fa > ce11.cb4.maf
```
