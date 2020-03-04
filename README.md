# WGA_GPU 

The system has been tested on all the AWS G3 and P3 GPU instances with AMI Ubuntu Server 18.04 LTS (HVM), SSD Volume Type (ami-0fc20dd1da406780b (64-bit x86))


* Clone WGA_GPU repository (https://github.com/gsneha26/WGA_GPU)

```
    $ git clone https://github.com/gsneha26/WGA_GPU.git
    $ export PROJECT_DIR=$PWD/WGA_GPU
```

* Install dependencies

```
    $ ./install.sh
```

* Run LASTZ_GPU

```
    $ run_lastz_gpu.sh target query "[options]"
```

* For a list of options 

```
    $ run_lastz_gpu.sh --help
```

* For running an example

```
    $ gunzip $PROJECT_DIR/data/ce11.fa.gz
    $ gunzip $PROJECT_DIR/data/cb4.fa.gz
    $ run_lastz_gpu.sh $PROJECT_DIR/data/ce11.fa $PROJECT_DIR/data/cb4.fa > ce11.cb4.maf
```
