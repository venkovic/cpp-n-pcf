/* 
 * File:   RealSignals.cpp
 * Author: Nicolas Venkovic
 * email: nvenkov1@jhu.edu
 */

#include "RealSignals.h"
#include <cstdlib>
#include <iostream>

RealSignal1D::RealSignal1D() {
}

RealSignal1D::~RealSignal1D() {
}

void RealSignal1D::initialize(double *f, int nx) {
    this->f=f;
    this->nx=nx;
}

double& RealSignal1D::operator() (int x) {
    if (x<0) x+=nx;
    return f[x];
}

RealSignal2D::RealSignal2D() {
}

RealSignal2D::~RealSignal2D() {
}

void RealSignal2D::initialize(double *f, int nx, int ny) {
    this->f=f;
    this->nx=nx;
    this->ny=ny;
}

double& RealSignal2D::operator() (int x, int y) {
    if (x<0) x+=nx;
    if (y<0) y+=ny;
    return f[x*ny+y];
}

void RealSignal2D::read_from_file(string fname, int x0, int y0, int verb) {    
    /*
     READ IN DATA
     */
    int error=0;
    if (verb>1) cout << "\n2D slice of data: \n" << endl; 
    //
    // Open input file
    ifstream ip(fname);
    if (!ip.is_open()) {
        error=2;
        if (verb>0) printf("ERROR (%d) encountered while opening input file.\n",error);
    }
    else {        
        //
        // Read input file
        string line;
        int i=0;
        while (getline(ip,line)) {
            if ((i>=x0) and (i<x0+nx)) {
                int j=0;
                string token;
                int pos;
                while (j<y0+ny) {
                    pos=line.find(",");
                    if (pos<0) {
                        if (j!=y0+ny-1) {
                            error=3;
                            if (verb>0) printf("ERROR (%d) encountered while reading input file.\n",error);
                            break;
                        }
                        else {
                            pos=line.length();
                        }
                    }
                    token=line.substr(0,pos);                
                    if ((j>=y0) and (j<=y0+ny)) {
                        (*this)(i-x0,j-y0)=stof(token);
                        if (verb>1) {
                            // 
                            // Print data
                            if (j-y0+1<ny) {
                                cout << (*this)(i-x0,j-y0) << ",";
                            }
                            else {
                                cout << (*this)(i-x0,j-y0) << endl;
                            }
                        }
                    }
                    line.erase(0,pos+1);     
                    ++j;
                }
            }
            else if (i>=x0+nx) break;
            ++i;
            if (error>0) break;
        }
        if (i!=x0+nx) {
            error=4;
            if (verb>0) printf("ERROR (%d) encountered while reading input file.\n",error);
        }
    }
    cout << endl;
    //
    // Close input file
    ip.close();    
}

RealSignal4D::RealSignal4D() {
}

RealSignal4D::~RealSignal4D() {
}

void RealSignal4D::initialize(double *f, int nx1, int ny1, int nx2, int ny2) {
    this->f=f;
    this->nx1=nx1;
    this->ny1=ny1;
    this->nx2=nx2;
    this->ny2=ny2;
}

double& RealSignal4D::operator() (int x1, int y1, int x2, int y2) {
    if (x1<0) x1+=nx1;
    if (y1<0) y1+=ny1;    
    if (x2<0) x2+=nx2;
    if (y2<0) y2+=ny2;
    return f[y2+ny2*(x2+nx2*(y1+ny1*x1))];
}

RealSignal6D::RealSignal6D() {
}

RealSignal6D::~RealSignal6D() {
}

void RealSignal6D::initialize(double *f, int nx1, int ny1, int nx2, int ny2, int nx3, int ny3) {
    this->f=f;
    this->nx1=nx1;
    this->ny1=ny1;
    this->nx2=nx2;
    this->ny2=ny2;
    this->nx3=nx3;
    this->ny3=ny3;
}

double& RealSignal6D::operator() (int x1, int y1, int x2, int y2, int x3, int y3) {
    if (x1<0) x1+=nx1;
    if (y1<0) y1+=ny1;    
    if (x2<0) x2+=nx2;
    if (y2<0) y2+=ny2;
    if (x3<0) x3+=nx3;
    if (y3<0) y3+=ny3;
    return f[y3+ny3*(x3+nx3*(y2+ny2*(x2+nx2*(y1+ny1*x1))))];
}