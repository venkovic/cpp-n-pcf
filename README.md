## npcf

##### A C++ program for estimation and parametric inference of n-point correlation functions of bi-dimensional random signals sampled on regular grids.

Author: Nicolas Venkovic

email: nvenkov1@jhu.edu



#### Compiler and dependencies:

 - g++, FFTW3, Eigen.
 - Enable g++ to link FFTW3 shared libraries with the option -lfftw3. 
 - Download [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) header and source files. Copy the Eigen folder as is in the current directory. 

#### Installation: 

```bash
$ make -f Makefile CONF=Release
```

#### Usage:

```bash
$ ./dist/GNU-Linux/npcf nx ny x0 y0 verb data.in
```

data.in: 2D csv data file.

nx, ny: Size of the 2D sub-slice of data to analyze.

x0, y0: upper left starting point of sub-slice in data file.

verb: Controls display settings.

 - 0: No output.
 - 1: Error messages only.
- 2: Error messages / 2D slice of data / Comparison of S2 and S3 results.

#### Example:

```bash
$ ./dist/GNU-Linux/npcf 30 30 0 0 2 im00.csv

s2(0,0) = 0.668750, s2(1,0) = 0.626250, s2(2,0) = 0.590625
s2(0,0) = 0.668750, s2(0,1) = 0.631250, s2(0,2) = 0.602500
s3(0,0,0,0) = 0.668750, s3(1,0,0,0) = 0.626250, s3(2,0,0,0) = 0.590625
s3(0,0,0,0) = 0.668750, s3(0,1,0,0) = 0.631250, s3(0,2,0,0) = 0.602500
s3(0,0,0,0) = 0.668750, s3(0,0,1,0) = 0.626250, s3(0,0,2,0) = 0.590625
s3(0,0,0,0) = 0.668750, s3(0,0,0,1) = 0.631250, s3(0,0,0,2) = 0.602500


s2(0,0) = 0.668750, s2(-1,0) = 0.626250, s2(-2,0) = 0.590625
s2(0,0) = 0.668750, s2(0,-1) = 0.631250, s2(0,-2) = 0.602500
s3(0,0,0,0) = 0.668750, s3(-1,0,0,0) = 0.626250, s3(-2,0,0,0) = 0.590625
s3(0,0,0,0) = 0.668750, s3(0,-1,0,0) = 0.631250, s3(0,-2,0,0) = 0.602500
s3(0,0,0,0) = 0.668750, s3(0,0,-1,0) = 0.626250, s3(0,0,-2,0) = 0.590625
s3(0,0,0,0) = 0.668750, s3(0,0,0,-1) = 0.631250, s3(0,0,0,-2) = 0.602500

```

#### Formats of input



#### Formats of output files

##### Anisotropic estimators:

- foo.s2 : Complete anisotropic 2-pcf estimator.

  ```
  nx,ny
  s2(-nx,0)  , s2(-nx,1)  , ...., s2(-nx,ny-1)  , s2(-nx,ny)
  s2(-nx+1,0), s2(-nx+1,1), ...., s2(-nx+1,ny-1), s2(-nx+1,ny)
  s2(-nx+2,0), s2(-nx+2,1), ...., s2(-nx+2,ny-1), s2(-nx+2,ny)
     :             :                  :               :
     :             :                  :               :
  s2(-1,0)   , s2(-1,1)   , ...., s2(-1,ny-1)   , s2(-1,ny)
  s2(0,0)    , s2(0,1)    , ...., s2(0,ny-1)    , s2(0,ny)
  s2(1,0)    , s2(1,1)    , ...., s2(1,ny-1)    , s2(1,ny)
     :             :                  :               :
     :             :                  :               :
  s2(nx-1,0) , s2(nx-1,1) , ...., s2(nx-1,ny-1) , s2(nx-1,ny)
  s2(nx,0)   , s2(nx,1)   , ...., s2(nx,ny-1)   , s2(nx,ny)
  ```


- foo.s3 : Anisotropic 3-pcf estimator of point configurations with fixed opening and rotation angles.

  ```
  nx,ny
  x1,y1,x2,y2
  s3(0,0,0,0)                  ,
  s3(dx1,dy1,0,0)              ,
  s3(2*dx1,2*dy1,0,0)          , s3(2*dx1,2*dy1,0,0),
        :
        :     
  s3((nx-1)*dx1,(nx-1)*dy1,0,0), 
  s3(nx*dx1,nx*dy1,0,0)        ,
  ```


- foo.s4 : (?)-tropic 4-pcf estimator.

  ```
  nx,ny
  x1,y1,x2,y2
  s3(0,0,0,0)                  ,
  s3(dx1,dy1,0,0)              ,
  s3(2*dx1,2*dy1,0,0)          , s3(2*dx1,2*dy1,0,0),
        :
        :     
  s3((nx-1)*dx1,(nx-1)*dy1,0,0), 
  s3(nx*dx1,nx*dy1,0,0)        ,
  ```

##### Isotropic estimators:

- foo.iso-s2 : Isotropic 2-pcf estimator.

  ```
  -
  ```


- foo.iso-s3 : Isotropic 3-pcf estimator.

  ```
  -
  ```


- foo.iso-s4 : Isotropic 4-pcf estimator.

  ```
  -
  ```

#### Pending tasks:

 -  ERROR to fix when almost all data entries are equal. For example, try

    ```bash
    $ ./dist/GNU-Linux/npcf 12 11 0 0 im00.csv
    ```

 -  Compute hs3 on minimum domain and copy values to other components instead of repeating calculations.

 -  Write subroutines to write output files.

 -  Complete implementation of S4(dx1,dy1,dx2,dy2,dx3,dy3).

 -  Verify implementation for odd nx and ny.

 -  Add inference subroutines.