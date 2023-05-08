#include "problem.hpp"
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include "problem.hpp"
#include <cmath>
using namespace std;
double Pt0(double x,double y){
    return sin(M_PI*x)*sin(M_PI*y);
}
double Vxt0(double x,double y){
    return 0;
}
double Vyt0(double x,double y){
    return 0;
}
double F(double x, double y, double t){
	return 0;
}
double Csq(double x,double y){
    return 1;
}
double Sigma_x(double x){
    return 0;
}
double P(double x,double y,double t){
    return cos(M_PI*M_SQRT2*t)*sin(M_PI*x)*sin(M_PI*y);
}
void init(){
    cerr << "Test 1" <<endl; 
    int dir_status = mkdir("data", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(dir_status != 0 && errno != EEXIST){
        cerr << "Failed to create data directory" <<endl;
        exit(1);
    }
}
typedef std::vector<std::vector<double>> grid;
double max_diff = 0;
int M = 10000; // Time 
int N = 100; // X Y
double Time = 1.0;
void output_layer(int m,double t, const grid &P1,const grid &Vx,const grid &Vy){
    double dt = 1.0/M;
    double dx = 1.0/N;
    double dy = dx;
    for(int j=0;j<N;j++){
        for(int i=0;i<N;i++){
            double diff = fabs(P1[i][j] - P(i*dx,j*dy,t) );
            max_diff = fmax(max_diff,diff);
        } 
    }
    //std::cout <<std::endl<<std::endl;
    int T=100;
    int print_T = 5*T;
    if(m%T)
        return;
    char fname[128];
    snprintf(fname,128, "data/%04d.pgm",m/T);
    if(m%print_T==0)
        fprintf(stderr,"[%d/%d] Computing image %s\033[0K\r",m/T,M/T,fname);
    std::ofstream data_stream(fname);
    if(!data_stream.is_open()){
        cerr <<"\nCan't open data file" <<endl;
        exit(1);
    }
    data_stream<<"P5"<<endl;
    data_stream<<N<<' '<<N<<endl;
    data_stream<<255<<endl;
    for(int j=0;j<N;j++){
        for(int i=0;i<N;i++){
            char byte = fmax(fmin(((P1[i][j]+1)/2)*255,255),0);
            data_stream<< byte;
        } 
    }
    //std::cout <<std::endl<<std::endl;
    data_stream.close();
    //if(m%100==0)
    //    fprintf(stderr,"[%d/%d] max_diff = %lf\033[0K\r",m,M,max_diff);
    
}
void finished(){

        printf("max_diff = %lf\n",max_diff);
}
