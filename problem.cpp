#include "problem.hpp"
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;
double Pt0(double x,double y){
    x-=0.5;
    y-=0.5;
	
    if(x*x + y*y < 1.0/100)
        return 0;
    return 0;
}
double Vxt0(double x,double y){
    return 0;
}
double Vyt0(double x,double y){
    return 0;
}
double F(double x, double y, double t){
	x-=0.5;
	y-=0.5;
    double R = 1.0/20;
	if(x*x + y*y<R*R)
		return cos(M_PI*2*10*t+sqrt(x*x+y*y)/R*M_PI)*10;
	return 0;
}
double Csq(double x,double y){
	x-=0.5;
	y-=0.5;
	if(fabs(y)<0.25)
		return 1;
//	if(fabs(y)>0.25)
//		return 2;
    return 0.15;
}
double Sigma_x(double x){
    return (x<0.2||x>0.8)*exp(8*fabs(x-0.5));
}
void init(){
    int dir_status = mkdir("data", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(dir_status != 0 && errno != EEXIST){
        cerr << "Failed to create data directory" <<endl;
        exit(1);
    }
}
typedef std::vector<std::vector<double>> grid;
int M = 90000; // Time 
int N = 300; // X Y
void output_layer(int m,double t, const grid &P1){
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
        cerr <<"Can't open data file" <<endl;
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
    
}
void finished(){

}
