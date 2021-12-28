# fsi_piston

This is the C++ implementation of the fluid-structure interaction problem from [Emmanuel Lefrançois and Jean-Paul Boufflet (2010)](http://dx.doi.org/10.1137/090758313). (The C-Model for now).

It is based on the  [Matlab code](http://www.utc.fr/~elefra02/ifs/funcref.html) written by the authors.


<p align="center">
  <a href="https://github.com/azzeddinetiba/fsi_piston">
    <img src="References/figure_gas_chamber.jpg" alt="Logo" width="200" height="100">
  </a>

</p>

<p align="center">
  <a href="https://github.com/azzeddinetiba/fsi_piston">
    <img src="References/meshes.jpg" alt="Logo" width="320" height="40">
  </a>

</p>

## Prerequisites

- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) should be installed (i.e Eigen cmake files generated).
- [cmake](https://cmake.org/) installed.

## Building

Then, you can build:
```
cd build
cmake ..
make
```
The building process has been tested under Windows 10 with g++ 9.2.0, cmake 3.22.0 and Eigen-3.4.0.

## Usage

```
./fsi_piston.exe
```

If indeed the ``` Results/``` folder is in the project directory, results will be loaded there, i.e text files containing time history of the pressure, the velocity, timesteps and the density fields in the fluid domain. 

- If installing Eigen is ommitted, and only header files are meant to be used, You need to include them along with the header files, and indicate it to your compiler :
  ```sh
  gcc -I XXX/eigen_path/ ...
  ```
  Parameters in this case should be modified directly from the ``` include/configure.h/``` file, instead of modifying via cmake (see next part).

## Problem parameters
 You can tune the __spring mass__, __number of nodes__ in the fluid domain, and the coefficient to multiply the structural natural period to give the total simulation time (see the problem details in the paper), using respectively the three flags : 
   ```sh
  mass              nmesh              coeff
  ```

eg

   ```sh
  cmake -Dcoeff="1.2" -Dmass="1000" -Dnmesh="200" ..
  ```


