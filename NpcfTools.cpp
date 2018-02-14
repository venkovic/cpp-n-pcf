/* 
 * File:   NpcfTools.cpp
 * Author: Nicolas Venkovic
 * email:  nvenkov1@jhu.edu.com
 */

#include "NpcfTools.h"
#include <math.h>

//#include <algorithm>

#define REAL 0
#define IMAG 1

NpcfTools::NpcfTools(string fname) {
    this->fname=fname;
}

NpcfTools::~NpcfTools() { 
    // Free memory
    if (im_data != NULL) fftw_free(im_data);
    if (him_data != NULL) fftw_free(him_data);
    if (s2_data != NULL) fftw_free(s2_data);
    if (s3_data != NULL) fftw_free(s3_data);
    if (s4_data != NULL) fftw_free(s4_data);
}

int NpcfTools::read_file(int nx, int ny, int x0, int y0, int verb) {
    int error=0;
 
    this->nx=nx;
    this->ny=ny;
    nyh=ny/2+1;
    
    // Allocate memory for 2D real signal and read data
    im_data=(double *) fftw_malloc(sizeof(double)*nx*ny);  
    im.initialize(im_data,nx,ny);
    im.read_from_file(fname,x0,y0,verb);
    
    // Allocate memory for Eigen::Array and copy data
    im_arr.resize(nx,ny);
    for (int i=0;i<nx;i++) {
        for (int j=0;j<ny;j++) {
            im_arr(i,j)=im(i,j);
        }
    }
    
    // Allocate memory for DFT of 2D real signal
    him_data=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nx*nyh);
    him.initialize(him_data,nx,ny);
    
    // Define FFTW plan
    im_to_him=fftw_plan_dft_r2c_2d(nx,ny,im_data,him_data,FFTW_ESTIMATE);
    
    // Execute FFTW plan
    fftw_execute(im_to_him); 
    
    // Destroy FFTW plans
    fftw_destroy_plan(im_to_him);
    
    return error;
}

double NpcfTools::get_single_value_anisotropic_s2(int i, int j) {
    if (i>=0 && j>=0) {
        return (im_arr.block(i,j,nx-i,ny-j)*im_arr.block(0,0,nx-i,ny-j)).sum()/(nx-i)/(ny-j);
    }
    else if (i>=0 && j<0) {  
        return (im_arr.block(i,0,nx-i,ny+j)*im_arr.block(0,-j,nx-i,ny+j)).sum()/(nx-i)/(ny+j);
    }
    else if (i<0 && j>=0) {
        return (im_arr.block(0,j,nx+i,ny-j)*im_arr.block(-i,0,nx+i,ny-j)).sum()/(nx+i)/(ny-j);
    }
    else {
        return (im_arr.block(0,0,nx+i,ny+j)*im_arr.block(-i,-j,nx+i,ny+j)).sum()/(nx+i)/(ny+j);
    }
}

double NpcfTools::get_single_value_anisotropic_s3(int i, int j, int k, int l) {
    double s3=0;
    int x0l=max(0,max(i,k));
    int x1l=x0l-i;
    int x2l=x0l-k;
    int dx=min(nx,min(nx+i,nx+k))-x0l;
    int y0l=max(0,max(j,l));
    int y1l=y0l-j;
    int y2l=y0l-l;
    int dy=min(ny,min(ny+j,ny+l))-y0l; 
    return (im_arr.block(x0l,y0l,dx,dy)*im_arr.block(x1l,y1l,dx,dy)*im_arr.block(x2l,y2l,dx,dy)).sum()/dx/dy;
}

void NpcfTools::get_anisotropic_map_s2(int nx, int ny) {
    
    s2_arr.resize(2*nx+1,ny+1);
    
    for (int i=-nx;i<=nx;i++) {
        for (int j=0;j<=ny;j++) {
            s2_arr(nx+i,j)=get_single_value_anisotropic_s2(i,j);
        }
    }
}

void NpcfTools::get_anisotropic_map_s3(int nx, int ny, int dx1, int dy1, int dx2, int dy2) {

    
    s3_arr.resize(nx+1,ny+1);

    for (int i=0;i<=nx;i++) {
        for (int j=0;j<=ny;j++) {
            s3_arr(i,j)=get_single_value_anisotropic_s3(i*dx1,i*dy1,j*dx2,j*dy2);
        }
    }
}


int NpcfTools::get_full_anistropic_s2_by_FFT() {
    int error=0;

    // Allocate memory for estimate S2
    s2_data=(double *) fftw_malloc(sizeof(double)*nx*ny);
    s2.initialize(s2_data,nx,ny);
  
    // Allocate memory for DFT of estimate S2, i.e. hS2
    hs2_data=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nx*nyh);
    hs2.initialize(hs2_data,nx,ny);

    // Define FFTW plan for transformation hS2 -> S2
    hs2_to_s2=fftw_plan_dft_c2r_2d(nx,ny,hs2_data,s2_data,FFTW_ESTIMATE);

    // Compute hS2 from him, component by component   
    for (int i=-nx/2+1;i<=nx/2;i++) {
        for (int j=0;j<=ny/2+1;j++) {
            hs2(i,j,REAL)=(him(i,j,REAL)*him(i,j,REAL)+him(i,j,IMAG)*him(i,j,IMAG))/pow(nx,2)/pow(ny,2);
            hs2(i,j,IMAG)=0;
        }
    }    
    
    // Execute FFTW plan
    fftw_execute(hs2_to_s2);  
    
    // Destroy FFTW plan
    fftw_destroy_plan(hs2_to_s2);    
    
    // Free memory
    fftw_free(hs2_data);   
    
    return error;
}

int NpcfTools::get_full_anistropic_s3_by_FFT() {
    int error=0;
    
    // Allocate memory for estimate S3
    s3_data=(double *) fftw_malloc(sizeof(double)*nx*ny*nx*ny);
    s3.initialize(s3_data,nx,ny,nx,ny);
    
    // Allocate memory for DFT of estimate S3, i.e. hS3
    hs3_data=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nx*ny*nx*nyh);
    hs3.initialize(hs3_data,nx,ny,nx,ny);
    
    // Define FFTW plan for transformation hS3 -> S3
    int *n= (int *)fftw_malloc(4*sizeof(int)); 
    n[0]=nx;
    n[1]=ny;
    n[2]=nx;
    n[3]=ny;
    hs3_to_s3=fftw_plan_dft_c2r(4,n,hs3_data,s3_data,FFTW_ESTIMATE);
    
    // Compute hS3 from him, component by component
    double a,b,c,d,e,f;
    for (int i=-nx/2+1;i<=nx/2;i++) {
      for (int j=-ny/2+1;j<=ny/2;j++) {
          if (j>=0) {
              a=him(i,j,REAL);
              b=him(i,j,IMAG);
          }
          else {
              if (i!=nx/2) {
                  a=him(-i,-j,REAL);
                  b=-him(-i,-j,IMAG);
              }
              else {
                  a=him(nx-i,-j,REAL);
                  b=-him(nx-i,-j,IMAG);
              }
          }
          for (int k=-nx/2+1;k<=nx/2;k++) {
              for (int l=0;l<=ny/2;l++) {
                  c=him(k,l,REAL);
                  d=him(k,l,IMAG);
                  int ik=i+k;
                  int jl=j+l;
                  if (ik<-nx/2+1) {
                      ik+=nx;
                  }
                  else if (ik>nx/2) {
                      ik-=nx;
                  }                
                  if (jl>ny/2) {
                      jl-=ny;
                  }                
                  if (jl>=0) {
                      e=him(ik,jl,REAL);
                      f=him(ik,jl,IMAG);
                  }
                  else {
                      if (ik!=nx/2) {
                          e=him(-ik,-jl,REAL);
                          f=-him(-ik,-jl,IMAG); 
                      }
                      else {
                          e=him(nx-ik,-jl,REAL);
                          f=-him(nx-ik,-jl,IMAG); 
                      }
                  }
                  hs3(i,j,k,l,REAL)=(a*c*e+a*d*f+b*c*f-b*d*e)/pow(nx,3)/pow(ny,3);
                  hs3(i,j,k,l,IMAG)=(-a*c*f+a*d*e+b*c*e+b*d*f)/pow(nx,3)/pow(ny,3);
              }
          }
      }
    }            
    
    // Execute FFTW plan
    fftw_execute(hs3_to_s3);    
    
    // Destroy FFTW plan
    fftw_destroy_plan(hs3_to_s3);
    
    // Free memory
    fftw_free(hs3_data);  
    
    return error;  
}

int NpcfTools::get_full_anistropic_s4_by_FFT() {
    int error=0;
    return error;
}