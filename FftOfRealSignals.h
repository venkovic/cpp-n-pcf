/*
 * File:   FftOfRealSignals.h
 * Author: Nicolas Venkovic
 * email: nvenkov1@jhu.edu
 */

#ifndef FFTOFREALSIGNALS_H
#define FFTOFREALSIGNALS_H

#include <fftw3.h>

using namespace std;

class FftOfRealSignal1D {
public:
    FftOfRealSignal1D();
    double& operator() (int, size_t);
    void initialize(fftw_complex *, int);
    virtual ~FftOfRealSignal1D();
private:
    friend class NpcfTools;
    fftw_complex *hf;
    int nx;
    int nh;
};

class FftOfRealSignal2D {
public:
    FftOfRealSignal2D();
    double& operator() (int, int, size_t);
    void initialize(fftw_complex *, int, int);
    virtual ~FftOfRealSignal2D();
private:
    friend class NpcfTools;
    fftw_complex *hf;    
    int nx, ny;
    int nh;
};

class FftOfRealSignal4D {
public:
    FftOfRealSignal4D();
    double& operator() (int, int, int, int, size_t);
    void initialize(fftw_complex *, int, int, int, int);
    virtual ~FftOfRealSignal4D();
private:
    fftw_complex *hf; 
    int nx1, ny1, nx2, ny2;
    int nh;
};

class FftOfRealSignal6D {
public:
    FftOfRealSignal6D();
    double& operator() (int, int, int, int, int, int, size_t);
    void initialize(fftw_complex *, int, int, int, int, int, int);
    virtual ~FftOfRealSignal6D();
private:
    fftw_complex *hf; 
    int nx1, ny1, nx2, ny2, nx3, ny3;
    int nh;
};

#endif /* FFTOFREALSIGNALS_H */