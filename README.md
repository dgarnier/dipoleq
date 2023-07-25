# dipoleq
Dipole equilibrium design code 

Can be used to show high beta equilibria in a dipole confined plasma.  In other words, solve the Grad-Shafranov equation in a machine with an "first closed flux surface" and no toroidal field.

The code is licensed under the MIT license, except for the ClibPDF code, which is licensed for free use in academic use.

If you use this code, please reference this paper as source.

[1] D. T. Garnier, J. Kesner, and M. E. Mauel, “Magnetohydrodynamic stability in a levitated dipole,” Phys. Plasmas, vol. 6, no. 9, pp. 3431–3434, Jan. 1999, doi: 10.1063/1.873601.

## Compilation

On ubuntu, the recommended build is for hdf4, zlib and libjpeg libraries to be installed:

`sudo apt install zlib1g-dev libjpeg-dev libhdf4-dev`

Then, the code can be compiled with:

`cmake `
