/* 
 * File:   NpcfTools.cpp
 * Author: Nicolas Venkovic
 * email:  nvenkov1@jhu.edu.com
 */

#include "NpcfTools.h"
#include <math.h>
#include <fstream>

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
    
    
    fftw_cleanup();
    
}

int NpcfTools::read_file(int nx, int ny, int x0, int y0, int verb) {
    int error=0;
 
    this->nx=nx;
    this->ny=ny;
    nyh=ny/2+1;
    
    dnx=nx;
    dny=ny;
    
    // Allocate memory for 2D real signal and read data
    im_data=(double *) fftw_malloc(sizeof(double)*nx*ny);  
    im.initialize(im_data,nx,ny);
    im.read_from_file(fname,x0,y0,verb);
    
    // Allocate memory for Eigen::Array and copy data
    if (im_arr.rows()!=nx || im_arr.cols()!=ny) allocate_im_array();
    
    //cout << "test:" << im(nx-1,ny-1) << endl;
    //cout << "test:" << im(-1,-1) << endl;
    //cout << "test:" << im(nx,ny) << endl;  
        
    // Allocate memory for DFT of 2D real signal
    //him_data=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nx*nyh);
    //him.initialize(him_data,nx,ny);
    
    // Define FFTW plan
    //im_to_him=fftw_plan_dft_r2c_2d(nx,ny,im_data,him_data,FFTW_ESTIMATE);
    
    // Execute FFTW plan
    //fftw_execute(im_to_him); 
    
    // Destroy FFTW plans
    //fftw_destroy_plan(im_to_him);
    
    return error;
}

void NpcfTools::allocate_im_array(){
    im_arr.resize(nx,ny);
    for (int i=0;i<nx;i++) {
        for (int j=0;j<ny;j++) {
            //im_arr(i,j)=im(i,j);
            im_arr(i,j)=im_data[j+i*ny];
        }
    }
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

void NpcfTools::get_anisotropic_map_s2(int nx, int ny, string fname) {
    
    s2_arr.resize(2*nx+1,ny+1);
    
    ofstream fout;
    if (fname!="") fout.open(fname);
    
    for (int i=-nx;i<=nx;i++) {
        for (int j=0;j<=ny;j++) {
            s2_arr(nx+i,j)=get_single_value_anisotropic_s2(i,j);
            if (fout) {
                fout << s2_arr(nx+i,j);
                if (j<ny) fout << ",";
            }
        }
        if (fout) fout << endl;
    }
    fout.close();
}

void NpcfTools::get_anisotropic_map_s3(int nx, int ny, int dx1, int dy1, int dx2, int dy2, string fname) {

    
    s3_arr.resize(nx+1,ny+1);

    ofstream fout;
    if (fname!="") fout.open(fname);    
    
    for (int i=0;i<=nx;i++) {
        for (int j=0;j<=ny;j++) {
            s3_arr(i,j)=get_single_value_anisotropic_s3(i*dx1,i*dy1,j*dx2,j*dy2);
            if (fout) {
                fout << s3_arr(i,j);
                if (j<ny) fout << ",";
            }
        }
        if (fout) fout << endl;
    }
}



int NpcfTools::get_full_anisotropic_s2_by_seq_FFT(int dnx, int dny, string fname) {        
    
    this->dnx=dnx;
    this->dny=dny;
    dnyh=dny/2+1;
    //
    nDomainX=nx/dnx;
    nDomainY=ny/dny;
    nDomains=nDomainX*nDomainY;   
    //
    // Allocate memory for anisotropic FFT-based estimator of S2
    s2_data=(double *) fftw_malloc(sizeof(double)*dnx*dny);
    hs2_data=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*dnx*dnyh);
    him_data=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*dnx*dnyh);    

    if (nDomains>1) {
        domainIsSegmented=true; 
        im_i_data=(double *) fftw_malloc(sizeof(double)*dnx*dny);
        hs2_sum_data=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*dnx*dnyh); 
    }
 
    if (domainIsSegmented)  { 
        //
        // Recursively attempt allocating memory for smaller domains
        if (s2_data==NULL || hs2_data==NULL || im_i_data==NULL || him_data==NULL || hs2_sum_data==NULL) {
            return 1;
        }
        s2.initialize(s2_data,dnx,dny);
        hs2.initialize(hs2_data,dnx,dny);
        //im_i.initialize(im_i_data,dnx,dny);
        him.initialize(him_data,dnx,dny);
        hs2_sum.initialize(hs2_sum_data,dnx,dny);
        //
        // Create FFTW plan
        im_to_him=fftw_plan_dft_r2c_2d(dnx,dny,im_i_data,him_data,FFTW_ESTIMATE);
    }        
   else {
        if (s2_data==NULL || hs2_data==NULL || him_data==NULL) {
            return 1;
        }        
        //
        // Memory successfully allocated
        s2.initialize(s2_data,dnx,dny);
        hs2.initialize(hs2_data,dnx,dny);
        him.initialize(him_data,dnx,dny);
        //
        // Create FFTW plan
        im_to_him=fftw_plan_dft_r2c_2d(dnx,dny,im_data,him_data,FFTW_ESTIMATE);
    }
       
    cout << nDomainX << " " << dnx << endl;
    cout << nDomainY << " " << dny << endl;
    cout << nDomains << endl;        
    
    //
    // Compute and add the contribution to hS2 of each sub-domain
    if (domainIsSegmented) {
        for (int iDomainX=0;iDomainX<nDomainX;iDomainX++) {
            for (int iDomainY=0;iDomainY<nDomainY;iDomainY++) { 
                
                cout << nDomainY*iDomainX+iDomainY+1 << " / " << nDomains << endl;
                
                
                cout << "stage 1" << endl;
                // Extract im_i, get him_i
                for (int i=0;i<dnx;i++) {
                    for (int j=0;j<dny;j++) {
                        //im_i(i,j)=im(iDomainX*nDomainX+i,iDomainY*nDomainY+j);
                        //im_i_data[j+i*dny]=im_arr(iDomainX*nDomainX+i,iDomainY*nDomainY+j);
                        im_i_data[j+i*dny]=im_data[(iDomainX*dnx+i)*ny+iDomainY*dny+j];
                    }
                }
                
                cout << "stage 2" << endl;
                
                fftw_execute(im_to_him);
                //
                // Compute
                
                cout << "stage 3" << endl;
                
                get_full_anistropic_s2_by_FFT();

                cout << "stage 4" << endl;
                //
                // Add contribution
                for (int k=0;k<dnx*dnyh;k++) {
                    hs2_sum_data[k][REAL]+=hs2_data[k][REAL];
                    hs2_sum_data[k][IMAG]+=hs2_data[k][IMAG];
                }
                cout << "stage 5" << endl;
            }            
        }
        //
        // Destroy FFTW plan and free memory
        fftw_destroy_plan(im_to_him);    
        fftw_free(im_i_data);    
    }
    else {
        //
        // Execute FFTW plan
        fftw_execute(im_to_him); 
        //
        // Compute
        get_full_anistropic_s2_by_FFT();
    }
  
    // Define FFTW plan for transformation hS2 -> S2
    if (domainIsSegmented) {     
        // Define FFTW plan for transformation hS3 -> S3
        //hs3_to_s3=fftw_plan_dft_c2r(4,dn,hs3_sum_data,s3_data,FFTW_ESTIMATE);    
        hs2_to_s2=fftw_plan_dft_c2r_2d(dnx,dny,hs2_sum_data,s2_data,FFTW_ESTIMATE);
    }
    else {
        // Define FFTW plan for transformation hS3 -> S3
        //hs3_to_s3=fftw_plan_dft_c2r(4,dn,hs3_data,s3_data,FFTW_ESTIMATE);  
        hs2_to_s2=fftw_plan_dft_c2r_2d(nx,ny,hs2_data,s2_data,FFTW_ESTIMATE);
    }    
    //
    // Execute FFTW plan
    fftw_execute(hs2_to_s2);  
    //    
    // Write output file
    
    double frac0=0;
    for (int i=0;i<nx;i++) {
        for (int j=0;j<ny;j++) {
            frac0+=im_arr(i,j);
        }
    }
    
    ofstream fout;
    if (fname!="") fout.open(fname);    
    
    for (int i=-dnx/2+1;i<=dnx/2;i++) {
        for (int j=0;j<=dny/2;j++) {
            fout << s2(i,j)/nDomains;
            if (j<dny/2) fout << ",";
        }
        if (i<dnx/2) fout << "\n";
    }
    fout.close();
    cout << endl << s2(0,0)/nDomains << endl;    
    
    //
    // Destroy FFTW plan and free memory
    fftw_destroy_plan(hs2_to_s2);    
    fftw_free(hs2_data);       
    if (hs2_sum_data != NULL) fftw_free(hs2_sum_data);        
    
    return 0;
}

int NpcfTools::get_full_anisotropic_s3_by_seq_FFT(int dnx, int dny, string fname) {
    this->dnx=dnx;
    this->dny=dny;
    dnyh=dny/2+1;
    //
    nDomainX=nx/dnx;
    nDomainY=ny/dny;
    nDomains=nDomainX*nDomainY;    
    //
    //
    // Allocate memory for anisotropic FFT-based estimator of S2
    s3_data=(double *) fftw_malloc(sizeof(double)*dnx*dny*dnx*dny);
    hs3_data=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*dnx*dny*dnx*dnyh);
    him_data=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*dnx*dny*dnx*dnyh);
    
    if (nDomains>1) {
        domainIsSegmented=true; 
        im_i_data=(double *) fftw_malloc(sizeof(double)*dnx*dny);
        hs3_sum_data=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*dnx*dny*dnx*dnyh); 
    }
 
    
    if (domainIsSegmented)  { 
        //
        // Recursively attempt allocating memory for smaller domains
        if (s3_data==NULL || hs3_data==NULL || im_i_data==NULL || him_data==NULL || hs3_sum_data==NULL) {
            return 1;
        }
        s3.initialize(s3_data,dnx,dny,dnx,dny);
        hs3.initialize(hs3_data,dnx,dny,dnx,dny);
        //im_i.initialize(im_i_data,dnx,dny);
        him.initialize(him_data,dnx,dny);
        hs3_sum.initialize(hs3_sum_data,dnx,dny,dnx,dny);
        //
        // Create FFTW plan
        im_to_him=fftw_plan_dft_r2c_2d(dnx,dny,im_i_data,him_data,FFTW_ESTIMATE);
    }    
    else {
        if (s3_data==NULL || hs3_data==NULL || him_data==NULL) {
            return 1;
        }
        //
        // Memory successfully allocated
        s3.initialize(s3_data,dnx,dny,dnx,dny);
        hs3.initialize(hs3_data,dnx,dny,dnx,dny);
        him.initialize(him_data,dnx,dny);
        //
        // Create FFTW plan
        im_to_him=fftw_plan_dft_r2c_2d(dnx,dny,im_data,him_data,FFTW_ESTIMATE);
    }
       
    cout << nDomainX << " " << dnx << endl;
    cout << nDomainY << " " << dny << endl;
    cout << nDomains << endl;    
    
    //
    // Compute and add the contribution to hS2 of each sub-domain
    if (domainIsSegmented) {
        for (int iDomainX=0;iDomainX<nDomainX;iDomainX++) {
            for (int iDomainY=0;iDomainY<nDomainY;iDomainY++) { 
                
                cout << nDomainY*iDomainX+iDomainY+1 << " / " << nDomains << endl;
                
                
                cout << "stage 1" << endl;
                //
                // Extract im_i, get him_i
                for (int i=0;i<dnx;i++) {
                    for (int j=0;j<dny;j++) {
                        //im_i(i,j)=im(iDomainX*nDomainX+i,iDomainY*nDomainY+j);
                        //im_i_data[j+i*dny]=im_arr(iDomainX*nDomainX+i,iDomainY*nDomainY+j);
                        im_i_data[j+i*dny]=im_data[(iDomainX*dnx+i)*ny+iDomainY*dny+j];
                    }
                }
                
                cout << "stage 2" << endl;
                
                fftw_execute(im_to_him);
                //
                // Compute
                
                cout << "stage 3" << endl;
                
                get_full_anistropic_s3_by_FFT();

                cout << "stage 4" << endl;
                //
                // Add contribution
                for (int k=0;k<dnx*dny*dnx*dnyh;k++) {
                    hs3_sum_data[k][REAL]+=hs3_data[k][REAL];
                    hs3_sum_data[k][IMAG]+=hs3_data[k][IMAG];
                }
                cout << "stage 5" << endl;
            }            
        }
        
        
        //
        // Destroy FFTW plan and free memory
        fftw_destroy_plan(im_to_him);    
        fftw_free(im_i_data);    

    }
    else {
        //
        // Execute FFTW plan
        fftw_execute(im_to_him); 
        //
        // Compute
        get_full_anistropic_s3_by_FFT();
    }

    
    
    int *dn= (int *)fftw_malloc(4*sizeof(int)); 
    dn[0]=dnx;
    dn[1]=dny;
    dn[2]=dnx;
    dn[3]=dny;   
    //
    // Define FFTW plan for transformation hS2 -> S2
    if (domainIsSegmented) {     
        // Define FFTW plan for transformation hS3 -> S3
        hs3_to_s3=fftw_plan_dft_c2r(4,dn,hs3_sum_data,s3_data,FFTW_ESTIMATE);        
    }
    else {
        // Define FFTW plan for transformation hS3 -> S3
        hs3_to_s3=fftw_plan_dft_c2r(4,dn,hs3_data,s3_data,FFTW_ESTIMATE);  
    }
    //
    // Execute FFTW plan
    fftw_execute(hs3_to_s3);  
    //    
    // Write output file
    
    
    ofstream fout;
    if (fname!="") fout.open(fname);    
    
    for (int i=0;i<=dnx/2;i++) {
        for (int j=0;j<=dny/2;j++) {
            fout << s3(i,0,0,j)/nDomains;
            if (j<dny/2) fout << ",";
        }
        if (i<dnx/2) fout << "\n";
    }
    fout.close();
    cout << endl << s3(0,0,0,0)/nDomains << endl;    
    
    //
    // Destroy FFTW plan and free memory
    fftw_destroy_plan(hs3_to_s3);    
    fftw_free(hs3_data);       
    if (hs3_sum_data != NULL) fftw_free(hs3_sum_data);    
}

void NpcfTools::get_full_anisotropic_s4_by_seq_FFT() {
}



int NpcfTools::get_full_anistropic_s2_by_FFT() {
    int error=0;

    for (int i=-dnx/2+1;i<=dnx/2;i++) {
        //for (int j=0;j<=dny/2+1;j++) {
        for (int j=0;j<=dny/2;j++) {
            hs2(i,j,REAL)=(him(i,j,REAL)*him(i,j,REAL)+him(i,j,IMAG)*him(i,j,IMAG))/pow(dnx,2)/pow(dny,2);
            hs2(i,j,IMAG)=0;
        }
    }    
    
    return error;
}


int NpcfTools::get_full_anistropic_s2_by_FFT_old() {
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
    
    
    // Compute hS3 from him, component by component
    double a,b,c,d,e,f;
    for (int i=-dnx/2+1;i<=dnx/2;i++) {
      for (int j=-dny/2+1;j<=dny/2;j++) {
          if (j>=0) {
              a=him(i,j,REAL);
              b=him(i,j,IMAG);
          }
          else {
              if (i!=dnx/2) {
                  a=him(-i,-j,REAL);
                  b=-him(-i,-j,IMAG);
              }
              else {
                  a=him(dnx-i,-j,REAL);
                  b=-him(dnx-i,-j,IMAG);
              }
          }
          for (int k=-dnx/2+1;k<=dnx/2;k++) {
              for (int l=0;l<=dny/2;l++) {
                  c=him(k,l,REAL);
                  d=him(k,l,IMAG);
                  int ik=i+k;
                  int jl=j+l;
                  if (ik<-dnx/2+1) {
                      ik+=dnx;
                  }
                  else if (ik>dnx/2) {
                      ik-=dnx;
                  }                
                  if (jl>dny/2) {
                      jl-=dny;
                  }                
                  if (jl>=0) {
                      e=him(ik,jl,REAL);
                      f=him(ik,jl,IMAG);
                  }
                  else {
                      if (ik!=dnx/2) {
                          e=him(-ik,-jl,REAL);
                          f=-him(-ik,-jl,IMAG); 
                      }
                      else {
                          e=him(dnx-ik,-jl,REAL);
                          f=-him(dnx-ik,-jl,IMAG); 
                      }
                  }
                  hs3(i,j,k,l,REAL)=(a*c*e+a*d*f+b*c*f-b*d*e)/pow(dnx,3)/pow(dny,3);
                  hs3(i,j,k,l,IMAG)=(-a*c*f+a*d*e+b*c*e+b*d*f)/pow(dnx,3)/pow(dny,3);
              }
          }
      }
    }            
    
    return error;  
}

int NpcfTools::get_full_anistropic_s3_by_FFT_old() {
    int error=0;
    
    // Allocate memory for estimate S3
    s3_data=(double *) fftw_malloc(sizeof(double)*nx*ny*nx*ny);
    if (s3_data == NULL) {
        printf("Could not allocate enough memory for s3_data\n");
    }
    s3.initialize(s3_data,nx,ny,nx,ny);
    
    // Allocate memory for DFT of estimate S3, i.e. hS3
    hs3_data=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nx*ny*nx*nyh);
    if (hs3_data == NULL) printf("Could not allocate enough memory for hs3_data\n");
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