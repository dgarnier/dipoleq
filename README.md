# dipoleq

Dipole equilibrium design code

Based on "TokaMac v2.0" by L. Bai and M. Mauel at Columbia University.

Can be used to show high beta equilibria in a dipole confined plasma. In other words, solve the Grad-Shafranov equation in a machine with an "first closed flux surface" and no toroidal field.

The code is licensed under the MIT license, except for the ClibPDF code, which no longer is available but was licensed for free use in academic use.

If you use this code, please reference this paper as source.

[1] D. T. Garnier, J. Kesner, and M. E. Mauel, “Magnetohydrodynamic stability in a levitated dipole,” Phys. Plasmas, vol. 6, no. 9, pp. 3431–3434, Jan. 1999, doi: 10.1063/1.873601.

## Command Line Usage Compilation

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

### Complilation

The python module can be compiled with the following command:
`pip install .`
This may install in your user directory or in a virtual environment. It will also install the `dipoleq` and `dipoleq_geqdsk` command line tools as well. This may become the preferred way to install the code in the future.

### Usage

Python usage is still in development. More on this soon.


[![Code style: Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json))](https://github.com/astral-sh/ruff)
