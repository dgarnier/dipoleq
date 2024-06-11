# DipolEq

[![Documentation Status][rtd-badge]][rtd-link]
[![Actions Status][actions-badge]][actions-link]
[![codecov][codecov-badge]][codecov-link]
[![PyPI version][pypi-version]][pypi-link]
[![PyPI platforms][pypi-platforms]][pypi-link]
[![Contributor Covenant][contrib-cov-badge]][contrib-cov-link]
[![Ruff][ruff-badge]][ruff-link]


<!-- SPHINX-START -->
Dipole equilibrium design code

Based on "TokaMac v2.0" by L. Bai and M. Mauel at Columbia University.

Can be used to show high beta equilibria in a dipole confined plasma. In other words, solve the Grad-Shafranov equation in a machine with an "first closed flux surface" and no toroidal field.

The code is licensed under the MIT license. (Except for the ClibPDF code, which no longer is available but was licensed for free use in academic use. This code is not linked into the python package version of this code.)

If you use this code, please reference this paper as source.

> D. T. Garnier, J. Kesner, and M. E. Mauel, “Magnetohydrodynamic  stability in a levitated dipole,” Phys. Plasmas, vol. 6, no. 9, pp. 3431–3434, Jan. 1999, doi: 10.1063/1.873601.

## C Code - Command Line Version
### Compilation

The original code is written in C (and some C++) and gives a command line interface, that outputs a PDF, HDF5 and other output files. The code can be compiled with cmake.  There are example codes for how to process hdf5 files (into, for example an EFIT gfile) in the `postprocess` directory.

### Dependencies

#### Ubuntu

The recommended build is for hdf5, zlib and libjpeg libraries to be installed:

```
sudo apt install cmake gcc zlib1g-dev libjpeg-dev libhdf5-dev
```

#### macOS

`cmake` and `hdf5` can be installed with macports or homebrew and should be enough to compile the code.

```
brew install cmake hdf5
```

or

```
sudo port install cmake hdf5
```

### Build and Install

The code can be configured and compiled with cmake. The recommended way to compile is to create a build directory and run cmake from there:

```
cd dipoleq
mkdir build
cd build
cmake ..
cmake --build .
sudo cmake --install .
```

### Usage

The code can be run from the command line with the following command:

```
dipoleq -f <input_file>
```

See the `Testing` directories for examples of input files [`*.in`].

There is a second command line tool that can be used to generate a EFIT gfile from
the output hdf5 file:

```
dipoleq_geqdsk <output.h5>
```

## Python Module Interface

A python version of the code is now available.  This relies on the `pybind11` library to wrap the C code into a python module. This wraps the main features of the code into a python module, although the outer loop of the code is now in python.  The python module is not yet fully documented, but the main features are available, and regularly tested against the original C code. For packaging, the "scikit-build-core" package is used to build the code using a standard `pyproject.toml` file.

### Complilation

The python module can be compiled with the following command:
`pip install .`
This may install in your user directory or in a virtual environment. It will also install the `dipoleq` and `dipoleq_geqdsk` command line tools as well, though these are actually calls to the python / pybind11 version of the code and no longer include ClibPDF or HDF5 directly. This may become the preferred way to install the code in the future.

### Usage

Python usage is still in development. More on this soon.


<!-- prettier-ignore-start -->
[actions-badge]:            https://github.com/dgarnier/dipoleq/workflows/Build%20Wheels/badge.svg
[actions-link]:             https://github.com/dgarnier/dipoleq/actions
[codecov-badge]:            https://codecov.io/gh/dgarnier/dipoleq/branch/main/graph/badge.svg?token=ZLbQzIvyG8
[codecov-link]:             https://codecov.io/gh/dgarnier/dipoleq
[contrib-cov-badge]:        https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg
[contrib-cov-link]:         https://www.contributor-covenant.org/version/2/1/code_of_conduct/
[pypi-link]:                https://pypi.org/project/dipoleq/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/dipoleq
[pypi-version]:             https://badge.fury.io/py/dipoleq.svg
[rtd-badge]:                https://readthedocs.org/projects/dipoleq/badge/?version=latest
[rtd-link]:                 https://dipoleq.readthedocs.io/en/latest/?badge=latest
[ruff-badge]:               https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
[ruff-link]:                https://github.com/astral-sh/ruff
<!-- prettier-ignore-end -->
