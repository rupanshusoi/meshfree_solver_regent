# An Implicitly Parallel Meshfree Solver in Regent

## Overview
This project implements a meshfree solver for inviscid, compressible fluid flows using q-LSKUM in the implicitly parallel, high-level programming language Regent. This is a first step towards our aim of developing a hybrid, production-level CFD code that is able to utilize the full power of today's heterogeneous supercomputers by simultaneously executing on both CPUs and GPGPUs.

## Dependencies
Python 3.x  
[Legion](https://github.com/StanfordLegion/legion)


## Installation

### For PAW-ATM 2020
All tests were done with the `Control Replication` branch of Legion at `0baf4f44b2e744af26565e2b36dba7fea386d857` installed with OpenMP support.

## Usage
Most parameters, including number of partitions, input file path, physical constants etc, can be configured inside `src/config.rg`. To run, do ```python3 ../legion/language/regent.py src/meshfree_solver.rg``` followed by the required flags. We document some useful flags below. It is generally recommended to launch a single Regent (Legion) process per node.

| SNo. | Flag               | Effect                                       |
|------|--------------------|----------------------------------------------|
|    1 |       `-fopenmp 0` |              Suppress OpenMP code generation |
|    2 |         `-fcuda 0` |                Suppress CUDA code generation |
|    3 |        `-ll:cpu N` |                   Use `N` CPU cores per node |
|    4 |      `-ll:csize m` |                   Use `m` MB of RAM per node |
|    5 |        `-ll:gpu N` |                        Use `N` GPUs per node |
|    6 |      `-ll:fsize m` | Use `m` MB of GPU framebuffer memory per GPU |
|    7 | `-level runtime=5` |               Suppress most runtime warnings |

Further, the number of iterations and inner iterations can be set by appending `--iter X --inner-iter Y` to the flags. We generally set `X = 1000, Y = 3`.

For CPU runs (without OpenMP), the number of partitions should be equal to the number of CPU cores for optimal performance.

For Regent + OpenMP, we generally do `-ll:ocpu 1 -ll:onuma 1 -ll:othr 30` which gives us 30 OpenMP cores per socket. Additionally, we have noticed optimal performance with number of partitions equal to the number of sockets (generally two per node).

## Support
Please contact the author for any help in running the code or obtaining larger input files. Due to GitHub's file size contraints, only small input files can be provided here.

## Author
Rupanshu Soi, Department of Computer Science, Birla Institute of Technology and Science, Pilani at Hyderabad, India.

Email: <f20180294@hyderabad.bits-pilani.ac.in> or <rupanshusoi@gmail.com>
