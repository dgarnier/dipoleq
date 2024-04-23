# dipoleq
Dipole equilibrium design code 

Based on "TokaMac v2.0" by L. Bai and M. Mauel at Columbia University.

Can be used to show high beta equilibria in a dipole confined plasma.  In other words, solve the Grad-Shafranov equation in a machine with an "first closed flux surface" and no toroidal field.

The code is licensed under the MIT license, except for the ClibPDF code, which no longer is available but was licensed for free use in academic use.

If you use this code, please reference this paper as source.

[1] D. T. Garnier, J. Kesner, and M. E. Mauel, “Magnetohydrodynamic stability in a levitated dipole,” Phys. Plasmas, vol. 6, no. 9, pp. 3431–3434, Jan. 1999, doi: 10.1063/1.873601.

## Compilation

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

The code can be configured and compiled with cmake.  The recommended way to compile is to create a build directory and run cmake from there:
```
cd dipoleq
mkdir build
cd build
cmake ..
cmake --build .
sudo cmake --install .
```