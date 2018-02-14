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
    // Read slice of data from input file
    int read_file(int, int, int, int, int=1);
    //
    // Naive estimator of single valued anisotropic S2, S3, and S4
    double get_single_value_anisotropic_s2(int, int);
    double get_single_value_anisotropic_s3(int, int, int, int);
    //double get_single_value_anisotropic_s4(int, int, int, int, int ,int);  
    //
    // Naive estimator of anisotropic maps of S2, S3, and S4
    void get_anisotropic_map_s2(int, int);
    void get_anisotropic_map_s3(int, int, int, int, int, int);
    //double get_full_anisotropic_s4(int, int, int, int, int ,int); 
    //
    // FFT-based estimators of complete anisotropic map of S2, S3 and S4
    int get_full_anistropic_s2_by_FFT();
    int get_full_anistropic_s3_by_FFT();
    int get_full_anistropic_s4_by_FFT();
     
private:
    //
    // Data input file name
    string fname;
    //
    // Dimensions of slice of data to extract from file
    int nx, ny, nyh;
    //
    // Eigen::Array of image used for faster operations in non-FFT-based estimators
    Eigen::ArrayXXf im_arr;   
    //
    // Wrappers of im and FFT-based estimators of S2, S3 and S4 w/ user-friendly indexation enabled
    //(Data is stored in row-major dynamic arrays im_data, s2_data, ...)
    RealSignal2D im;
    RealSignal2D s2;
    RealSignal4D s3;
    RealSignal6D s4;
    //
    // Data structures used for FFT-based estimators 
    double *im_data = NULL;
    fftw_complex *him_data = NULL;
    FftOfRealSignal2D him;
    fftw_plan im_to_him;
    //
    // Data structures for FFT-based estimators of S2
    double *s2_data = NULL; 
    fftw_complex *hs2_data; 
    FftOfRealSignal2D hs2;
    fftw_plan hs2_to_s2;
    //
    // Data structures for FFT-based estimators of S3
    double *s3_data = NULL; 
    fftw_complex *hs3_data; 
    FftOfRealSignal4D hs3;
    fftw_plan hs3_to_s3;
    //
    // Data structures for FFT-based estimators of S4
    double *s4_data = NULL; 
    fftw_complex *hs4_data; 
    FftOfRealSignal6D hs4;
    fftw_plan hs4_to_s4;
};

#endif /* NPCFTOOLS_H */