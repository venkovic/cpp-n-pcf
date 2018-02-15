/* 
 * File:   main.cpp
 * Author: Nicolas Venkovic
 * email:  nvenkov1@jhu.edu.com
 */

#include "NpcfTools.h"

using namespace std;

int main(int argc, char** argv) {
    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    int x0 = atoi(argv[3]);
    int y0 = atoi(argv[4]);
    int verb = atoi(argv[5]);
    string fname = argv[6];


    int nyh=ny/2+1;
    double *s3_data = NULL; 
    fftw_complex *hs3_data = NULL; 
    FftOfRealSignal4D hs3;
    fftw_plan hs3_to_s3;

    // Test allocation
    s3_data=(double *) fftw_malloc(sizeof(double)*nx*ny*nx*ny);

    if (s3_data==NULL) cout << "Could not allocate memory for s3_data\n";
    
    fftw_plan test_plan = NULL;
   
    if (test_plan == NULL) cout << "plan not created" << endl;
    //s3.initialize(s3_data,nx,ny,nx,ny);

    // Allocate memory for DFT of estimate S3, i.e. hS3
    hs3_data=(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*nx*ny*nx*nyh);

    if (hs3_data==NULL) cout << "Could not allocate memory for hs3_data\n";

    
    if (s3_data != NULL) fftw_free(s3_data);
    if (hs3_data != NULL) fftw_free(hs3_data);
    
    //hs3.initialize(hs3_data,nx,ny,nx,ny);        
        



    
    
    NpcfTools npcf(fname);
    //int error = npcf.read_file(nx,ny,x0,y0,verb);
    int error =1;
    
    
    if (!error) {
        int j=0;
        int k=0;     
        
        npcf.get_anisotropic_map_s2(60,60,"test.s2");
        npcf.get_anisotropic_map_s3(60,60,0,1,1,0,"test.s3");

        /*
        for (int i=-60;i<=60;i++) {
            for (int l=0;l<=60;l++) {
                //printf("%d %d %d %d\n",i,j,k,l);
                s2=npcf.get_single_value_anisotropic_s2(i,l);
                //s3=npcf.get_single_value_anisotropic_s3(i,j,k,l);
            }               
        }*/
    }
    return 0;
}