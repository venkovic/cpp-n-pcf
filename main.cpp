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
    
    NpcfTools npcf(fname);
    int error = npcf.read_file(nx,ny,x0,y0,verb);
        
    if (!error) {
        int j=0;
        int k=0;     
        
        
        
        
        int error_3pcf=npcf.get_full_anisotropic_s3_by_seq_FFT(100,100,"test_FFT_70by70.s3");
        //int error_2pcf=npcf.get_full_anisotropic_s2_by_seq_FFT(280,ny,"test_FFT_280by902_1001.s2");
        //int error_2pcf=npcf.get_full_anisotropic_s2_by_seq_FFT(nx,ny,"test_FFT_560by902_1001.s2");
        
        //npcf.get_anisotropic_map_s2(60,60,"test.s2");
        //npcf.get_anisotropic_map_s3(60,60,0,1,1,0,"test.s3");

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