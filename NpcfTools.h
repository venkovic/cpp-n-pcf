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
    void get_anisotropic_map_s2(int, int, string="");
    void get_anisotropic_map_s3(int, int, int, int, int, int, string="");
    //double get_full_anisotropic_s4(int, int, int, int, int ,int); 
    //
    // Sequential FFT-based estimators of complete anisotropic maps of S2, S3 and S4
    void get_full_anisotropic_s2_by_seq_FFT();
    void get_full_anisotropic_s3_by_seq_FFT();
    void get_full_anisotropic_s4_by_seq_FFT();
    //
    // FFT-based estimators of complete anisotropic map of S2, S3 and S4
    int get_full_anistropic_s2_by_FFT();
    int get_full_anistropic_s2_by_FFT_old();
    int get_full_anistropic_s3_by_FFT();
    int get_full_anistropic_s3_by_FFT_old();
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
    Eigen::ArrayXXf s2_arr;
    Eigen::ArrayXXf s3_arr; 
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
    fftw_plan im_to_him = NULL;
    
    double *im_i_data = NULL;    
    double *im_data_i = NULL;
    RealSignal2D im_i;    
    
    //
    // Data structures for FFT-based estimators of S2
    double *s2_data = NULL; 
    fftw_complex *hs2_sum_data = NULL; 
    FftOfRealSignal2D hs2_sum;
    
    fftw_complex *hs2_data = NULL; 
    FftOfRealSignal2D hs2;
    
    fftw_plan hs2_to_s2 = NULL;
    //
    // Data structures for FFT-based estimators of S3
    double *s3_data = NULL; 
    fftw_complex *hs3_sum_data = NULL; 
    FftOfRealSignal4D hs3_sum;
    
    fftw_complex *hs3_data = NULL; 
    FftOfRealSignal4D hs3;
    fftw_plan hs3_to_s3 = NULL;
    //
    // Data structures for FFT-based estimators of S4
    double *s4_data = NULL; 
    fftw_complex *hs4_data = NULL; 
    FftOfRealSignal6D hs4;
    fftw_plan hs4_to_s4 = NULL;
    
    
    void allocate_im_array();
    
    bool domainIsSegmented = false;
    int nDomainX, nDomainY, nDomains;
    int dnx, dny, dnyh;

};

#endif /* NPCFTOOLS_H */