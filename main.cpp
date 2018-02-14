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

        npcf.get_anisotropic_map_s2(60,60);
        //npcf.get_anisotropic_map_s3(60,60,0,1,1,0);

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