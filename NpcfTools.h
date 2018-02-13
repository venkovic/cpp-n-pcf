/* 
 * File:   NpcfTools.h
 * Author: Nicolas Venkovic
 * email: nvenkov1@jhu.edu
 */

#ifndef NPCFTOOLS_H
#define NPCFTOOLS_H

#include "FftOfRealSignals.h"
#include "RealSignals.h"

#include <vector>
#include "Eigen/Dense"

class NpcfTools {
public:
    NpcfTools(string);
    virtual ~NpcfTools();
    //
    // Function 
    int read_file(int, int, int, int, int=1);
    //
    // Complete anisotropic FFT-based estimators of S2, S3 and S4
    int get_full_anistropic_s2_by_FFT();
    int get_full_anistropic_s3_by_FFT();
    int get_full_anistropic_s4_by_FFT();
    //
    // Single value anisotropic naive estimator of S2, S3, and S4
    double get_single_value_anisotropic_s2(int, int);
    double get_single_value_anisotropic_s3(int, int, int, int);  
    //
    // Wrappers of the image, S2, S3 and S4 estimates enabling user-friendly indexation
    // while actual data is stored in row-major dynamic arrays
    RealSignal2D im;
    RealSignal2D s2;
    RealSignal4D s3;
    RealSignal6D s4;
    //
    // Eigen::Array of image used for faster operations
    Eigen::ArrayXXf im_arr;    
    
private:
    string fname;
    int nx, ny, nyh;
    
    double *im_data = NULL;
    fftw_complex *him_data = NULL;
    FftOfRealSignal2D him;
    fftw_plan im_to_him;

    double *s2_data = NULL; 
    fftw_complex *hs2_data; 
    FftOfRealSignal2D hs2;
    fftw_plan hs2_to_s2;
    
    double *s3_data = NULL; 
    fftw_complex *hs3_data; 
    FftOfRealSignal4D hs3;
    fftw_plan hs3_to_s3;
    
    double *s4_data = NULL; 
    fftw_complex *hs4_data; 
    FftOfRealSignal6D hs4;
    fftw_plan hs4_to_s4;
};

#endif /* NPCFTOOLS_H */