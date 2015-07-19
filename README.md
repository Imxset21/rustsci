# rustsci  [![Build Status](https://travis-ci.org/Imxset21/rustsci.svg?branch=master)](https://travis-ci.org/Imxset21/rustsci)[![Coverage Status](https://coveralls.io/repos/Imxset21/rustsci/badge.svg?branch=master&service=github)](https://coveralls.io/github/Imxset21/rustsci?branch=master)

Rust adaptation of SciPy, with a native generic array/matrix implementation.

`rustsci` uses (optional) OpenBLAS and LAPACK bindings to provide acceleration to certain operations.
Currently we focus only on the more common `f32` (single precision) and `f64` functions, with complex numbers TBD.

Besides Rust 1.1.0+, you'll currently need OpenBLAS and LAPACK libraries and headers installed, though these will be made optional in a future release.

On Ubuntu 14.04LTS, these can be installed as such:
`sudo apt-get install gfortran liblapack-dev libblas-dev libblas3gf libopenblas-dev`

To build, simply use `cargo`:
`cargo build`

To run tests:
`cargo test`
