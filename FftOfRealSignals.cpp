/* 
 * File:   FftOfRealSignals.cpp
 * Author: Nicolas Venkovic
 * email: nvenkov1@jhu.edu
 */

#include "FftOfRealSignals.h"

#include <cstdlib>
#include <iostream>

FftOfRealSignal1D::FftOfRealSignal1D() {
}

FftOfRealSignal1D::~FftOfRealSignal1D() {
}

void FftOfRealSignal1D::initialize(fftw_complex *hf, int nx) {
    this->hf=hf;
    this->nx=nx;
    nh=nx/2+1;
}

double& FftOfRealSignal1D::operator() (int w, size_t dtype) {
    if (w<0) w+=nx;
    return hf[w][dtype];
}

FftOfRealSignal2D::FftOfRealSignal2D() {
}

FftOfRealSignal2D::~FftOfRealSignal2D() {
}

void FftOfRealSignal2D::initialize(fftw_complex *hf, int nx, int ny) {
    this->hf=hf;
    this->nx=nx;
    this->ny=ny;
    nh=ny/2+1;
}

double& FftOfRealSignal2D::operator() (int wx, int wy, size_t dtype) {
    if (wx<0) wx+=nx;
    if (wy<0) wy+=ny;
    return hf[wx*nh+wy][dtype];
}

FftOfRealSignal4D::FftOfRealSignal4D() {
}

FftOfRealSignal4D::~FftOfRealSignal4D() {
}

void FftOfRealSignal4D::initialize(fftw_complex *hf, int nx1, int ny1, int nx2, int ny2) {
    this->hf=hf;
    this->nx1=nx1;
    this->ny1=ny1;
    this->nx2=nx2;
    this->ny2=ny2;
    nh=ny2/2+1;
}

double& FftOfRealSignal4D::operator() (int wx1, int wy1, int wx2, int wy2, size_t dtype) {
    if (wx1<0) wx1+=nx1;
    if (wy1<0) wy1+=ny1;    
    if (wx2<0) wx2+=nx2;
    if (wy2<0) wy2+=ny2;
    return hf[wy2+nh*(wx2+nx2*(wy1+ny1*wx1))][dtype];
}

FftOfRealSignal6D::FftOfRealSignal6D() {
}

FftOfRealSignal6D::~FftOfRealSignal6D() {
}

void FftOfRealSignal6D::initialize(fftw_complex *hf, int nx1, int ny1, int nx2, int ny2, int nx3, int ny3) {
    this->hf=hf;
    this->nx1=nx1;
    this->ny1=ny1;
    this->nx2=nx2;
    this->ny2=ny2;
    this->nx3=nx3;
    this->ny3=ny3;
    nh=ny3/2+1;
}

double& FftOfRealSignal6D::operator() (int wx1, int wy1, int wx2, int wy2, int wx3, int wy3, size_t dtype) {
    if (wx1<0) wx1+=nx1;
    if (wy1<0) wy1+=ny1;    
    if (wx2<0) wx2+=nx2;
    if (wy2<0) wy2+=ny2;
    if (wx3<0) wx3+=nx3;
    if (wy3<0) wy3+=ny3;
    return hf[wy3+nh*(wx3+nx3*(wy2+ny2*(wx2+nx2*(wy1+ny1*wx1))))][dtype];
}