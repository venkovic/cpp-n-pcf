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
    string fname = argv[5];
    
    NpcfTools npcf(fname);
    int verb=2;
    int error = npcf.read_file(nx,ny,x0,y0,verb);
    if (!error) {
        npcf.get_s2();     
        printf("\ns2(0,0) = %f, s2(1,0) = %f, s2(2,0) = %f\n",npcf.s2(0,0),npcf.s2(1,0),npcf.s2(2,0));
        printf("\ns2(0,0) = %f, s2(0,1) = %f, s2(0,2) = %f\n",npcf.s2(0,0),npcf.s2(0,1),npcf.s2(0,2));
        
        npcf.get_s3();     
        printf("\ns3(0,0,0,0) = %f, s3(1,0,0,0) = %f, s3(2,0,0,0) = %f\n",npcf.s3(0,0,0,0),npcf.s3(1,0,0,0),npcf.s3(2,0,0,0));
        printf("\ns3(0,0,0,0) = %f, s3(0,1,0,0) = %f, s3(0,2,0,0) = %f\n",npcf.s3(0,0,0,0),npcf.s3(0,1,0,0),npcf.s3(0,2,0,0));
        printf("\ns3(0,0,0,0) = %f, s3(0,0,1,0) = %f, s3(0,0,2,0) = %f\n",npcf.s3(0,0,0,0),npcf.s3(0,0,1,0),npcf.s3(0,0,2,0));
        printf("\ns3(0,0,0,0) = %f, s3(0,0,0,1) = %f, s3(0,0,0,2) = %f\n\n",npcf.s3(0,0,0,0),npcf.s3(0,0,0,1),npcf.s3(0,0,0,2));
        
        printf("\ns2(0,0) = %f, s2(-1,0) = %f, s2(-2,0) = %f\n",npcf.s2(0,0),npcf.s2(-1,0),npcf.s2(-2,0));
        printf("\ns2(0,0) = %f, s2(0,-1) = %f, s2(0,-2) = %f\n",npcf.s2(0,0),npcf.s2(0,-1),npcf.s2(0,-2));
        
        printf("\ns3(0,0,0,0) = %f, s3(-1,0,0,0) = %f, s3(-2,0,0,0) = %f\n",npcf.s3(0,0,0,0),npcf.s3(-1,0,0,0),npcf.s3(-2,0,0,0));
        printf("\ns3(0,0,0,0) = %f, s3(0,-1,0,0) = %f, s3(0,-2,0,0) = %f\n",npcf.s3(0,0,0,0),npcf.s3(0,-1,0,0),npcf.s3(0,-2,0,0));
        printf("\ns3(0,0,0,0) = %f, s3(0,0,-1,0) = %f, s3(0,0,-2,0) = %f\n",npcf.s3(0,0,0,0),npcf.s3(0,0,-1,0),npcf.s3(0,0,-2,0));
        printf("\ns3(0,0,0,0) = %f, s3(0,0,0,-1) = %f, s3(0,0,0,-2) = %f\n\n",npcf.s3(0,0,0,0),npcf.s3(0,0,0,-1),npcf.s3(0,0,0,-2));
    }
    return 0;
}