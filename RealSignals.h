/* 
 * File:   RealSignals.h
 * Author: Nicolas Venkovic
 * email: nvenkov1@jhu.edu
 */

#ifndef REALSIGNALS_H
#define REALSIGNALS_H

#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

class RealSignal1D {
public:
    RealSignal1D();
    virtual ~RealSignal1D();
    void initialize(double *, int);
    double& operator() (int);
private:
    double *f;
    int nx;
};

class RealSignal2D {
public:
    RealSignal2D();
    virtual ~RealSignal2D();
    void initialize(double *, int, int);
    double& operator() (int, int);
    void read_from_file(string, int, int, int=1);
private:
    double *f;    
    int nx, ny;
};

class RealSignal4D {
public:
    RealSignal4D();
    virtual ~RealSignal4D();
    void initialize (double *, int, int, int, int);
    double& operator() (int, int, int, int);
private:
    double *f; 
    int nx1, ny1, nx2, ny2;
};

class RealSignal6D {
public:
    RealSignal6D();
    virtual ~RealSignal6D();
    void initialize (double *, int, int, int, int, int, int);
    double& operator() (int, int, int, int, int, int);
private:
    double *f; 
    int nx1, ny1, nx2, ny2, nx3, ny3;
};

#endif /* REALSIGNALS_H */