## npcf

##### A C++ program for estimation and (anisotropic) parametric inference of n-point correlation functions of 2D signals sampled on regular grids.

Author: Nicolas Venkovic

email: nvenkov1@jhu.edu

#### Dependences:

 - FFTW3

#### Compiling instructions: 

```bash
>>> g++ 
```

#### Usage:

```bash
>>> ./npcf nx ny x0 y0 verb data.in
```

#### TO DO:

 -  ERROR to fix when almost all data entries are equal. For example, try

    ./npcf 12 11 0 0 npcf.in

 -  Use symmetry to copy components of hsn using symmetry.

 -  Write out sn using symmetry

 -  Complete implementation of S4(dx1,)

 -  Add inference subroutines.