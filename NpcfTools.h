/* 
 * File:   NpcfTools.h
 * Author: Nicolas Venkovic
 * email: nvenkov1@jhu.edu
 */

#ifndef NPCFTOOLS_H
#define NPCFTOOLS_H

#include "FftOfRealSignals.h"
#include "RealSignals.h"

class NpcfTools {
public:
    NpcfTools(string);
    virtual ~NpcfTools();
    //
    // Function 
    int read_file(int, int, int, int, int=1);
    //
    // Functions to compute estimates of S2, S3 and S4
    int get_s2();
    int get_s3();
    int get_s4();
    //
    // Wrappers of the image and S2, S3 and S4 estimates enabling user-friendly indexation
    RealSignal2D im;
    RealSignal2D s2;
    RealSignal4D s3;
    RealSignal6D s4;
    
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