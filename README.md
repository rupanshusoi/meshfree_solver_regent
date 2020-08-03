# An Implicitly Parallel Meshfree Solver in Regent

This project implements a meshfree solver for inviscid, compressible fluid flows using q-LSKUM in the implicitly parallel, high-level programming language Regent. This is a first step towards our aim of developing a hybrid, production-level CFD code that is able to utilize the full power of today's heterogeneous supercomputers by simultaneously executing on both CPUs and GPGPUs.

## Dependencies
Python 3.x  
[Legion](https://github.com/StanfordLegion/legion)


## Installation

### For PAW-ATM 2020
All tests were done with the `Control Replication` branch of Legion at `0baf4f44b2e744af26565e2b36dba7fea386d857` installed with OpenMP support.


## Author
Rupanshu Soi, Department of Computer Science, Birla Institute of Technology and Science, Pilani at Hyderabad, India.

Email: <f20180294@hyderabad.bits-pilani.ac.in> or <rupanshusoi@gmail.com>
