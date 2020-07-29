# SegAlign 

A Scalable GPU System for Pairwise Whole Genome Alignments

## Table of Contents

- [Overview](#overview)
- [Dependencies](#dependencies)
- [How to run SegAlign](#run)
    - [Running a test](#test)
    - [API specific test](#api-test)

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
    $ run_segalign ce11.fa cb4.fa > ce11.cb4.maf
```

For this branch, chr1.fa file is added. Running
```
  $ cd $PROJECT_DIR
  $ mkdir test
  $ cd test
  $ segalign $PROJECT_DIR/data/chr1.fa $PROJECT_DIR/data/chr1.fa $PROJECT_DIR/data/ --nogapped --debug
```

should result in 107581 HSPs.

### <a name="api-test"></a> API specific test

In folder `api_test`,
  * example.fa consists of the input sequence
  * example_hits.txt consists of the seed hits for the plus strand (use --strand=plus with segalign) for the sequence in example.fa, each line is an entry of the type seedHit (ref_start,query_start) 
  * example_hsps.txt consists of the output HSPs resulting after UngappedExtend(), each line is an entry of the type segment (ref_start,query_start,len,score) 

For the example,
  * call `InitializeUngappedExtension()` with 
  ```
    * num_gpu = 1
    * sub_mat = {91, -114, -31, -123, -1000, -1000, -100, -9100, -114,  100, -125,  -31, -1000, -1000, -100, -9100, -31, -125,  100, -114, -1000, -1000, -100, -9100, -123,  -31, -114,  91, -1000, -1000, -100, -9100, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -9100, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -9100, -100, -100, -100, -100, -1000, -1000, -100, -9100, -9100, -9100, -9100, -9100, -9100, -9100, -9100, -9100}; 
    * input_xdrop = 910
    * input_hspthresh = 3000
    * input_noentropy = false
  ```
  * create a device char array (e.g. `seq_tmp`) with the sequence from example.fa
  * generate compressed device char array (e.g. `seq_compressed`) as the output of `CompressSeq()` with `seq_tmp` as input
  * read the hits from example_hits.txt into a device seedHit array (e.g. `input_hits`) 
  * call `UngappedExtend()` with 
    ```
    * r_seq = seq_compressed
    * q_seq = seq_compressed
    * hits = input_hits 
    ```
    hsps in hsp_out device array should match the hsps in example_hsps.txt 
  * call `ShutdownUngappedExtension()`

To run with SegAlign,
```
    $ cd $PROJECT_DIR/api_test
    $ run_segalign example.fa example.fa --strand=plus --nogapped 
```

