#include "solver.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
using namespace std;
/* u_t = v_x
 * v_t = c^2*u_x + F(x,t)
 * u(x,0) = m1
 * v(x,0) = \int m2(\xi) d \xi
 * u(0,t) = m3(t)
 * u(a,t) = m4(t)
 */
#include "problem.hpp"
int main(){
    int M = 90000; // Time 
    int N = 300; // X Y
    
    double dt = 1.0/M;
    double dx = 1.0/N;
    double dy = dx;
    Solver solver(M,N);


    int dir_status = mkdir("data", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(dir_status != 0 && errno != EEXIST){
        cerr << "Failed to create data directory" <<endl;
        return 1;
    }


    int T=100;
    int print_T = 5*T;
    for(int m=0;m<M;m++){
        double t = m * dt;
        grid *P1_ptr = solver.next_frame(t);
        grid &P1 = *P1_ptr;
        if(m%T)
            continue;
        char fname[128];
        snprintf(fname,128, "data/%04d.pgm",m/T);
        if(m%print_T==0)
            fprintf(stderr,"[%d/%d] Computing image %s\033[0K\r",m/T,M/T,fname);
        std::ofstream data_stream(fname);
        if(!data_stream.is_open()){
            cerr <<"Can't open data file" <<endl;
            return 1;
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
    return 0;

}

