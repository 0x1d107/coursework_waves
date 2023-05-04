#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
#include "problem.hpp"
using namespace std;
/* u_t = v_x
 * v_t = c^2*u_x + F(x,t)
 * u(x,0) = m1
 * v(x,0) = \int m2(\xi) d \xi
 * u(0,t) = m3(t)
 * u(a,t) = m4(t)
 */
int main(int argc,char *argv[]){
    if(argc>1){
        M = atol(argv[1]);
        N = sqrt(M);
    }
    if(argc>2){
        N = atol(argv[2]);
    }
    grid P_a(N,std::vector<double>(N,0));
    grid P_b(N,std::vector<double>(N,0));
    grid Vx_a(N,std::vector<double>(N,0));
    grid Vx_b(N,std::vector<double>(N,0));
    grid Vy_a(N,std::vector<double>(N,0));
    grid Vy_b(N,std::vector<double>(N,0));
    grid Psi_a(N,std::vector<double>(N,0));
    grid Psi_b(N,std::vector<double>(N,0));


    grid &P0 = P_a;
    grid &P1 = P_b;
    grid &Vx0 = Vx_a;
    grid &Vx1 = Vx_b;
    grid &Vy0 = Vy_a;
    grid &Vy1 = Vy_b;
    grid &Psi0 = Psi_a;
    grid &Psi1 = Psi_b;
    
    double dt = 1.0/M;
    double dx = 1.0/N;
    double dy = dx;
    init();
    for(int j=0;j<N;j++){
        double y = j*dy;
        for(int i = 0; i < N; i++){
            double x = i*dx;
            P0[i][j] = Pt0(x,y);
            Vx0[i][j] = Vxt0(x,y);
            Vy0[i][j] = Vyt0(x,y);
        }
    }
    for(int m=0;m<M;m++){
        double t = m * dt;
        #pragma omp parallel for
        for(int j=0;j<N;j++){
            double y = j*dy;
            #pragma omp parallel for
            for(int i = 0; i < N; i++){
                double x = i*dx;
                if(i>=1&&j>=1&&i<N-1&&j<N-1){
                    P1[i][j] = P0[i][j] - Csq(x,y) * dt*(
                            (Vx0[i][j]-Vx0[i-1][j])/dx
                            + (Vy0[i][j] - Vy0[i][j-1])/dy 
                            ) + dt*F(x,y,t) - dt * Sigma_x(x)*P0[i][j] + dt*Psi0[i][j];
                }
                if(i<N-1&&j<N-1){
                    Vx1[i][j] = Vx0[i][j] - dt/dx *( P1[i+1][j] - P1[i][j]) - dt*Sigma_x(x)*Vx0[i][j];
                    Vy1[i][j] = Vy0[i][j] - dt/dy *( P1[i][j+1] - P1[i][j]);
                }
                if(i>=1&&j>=1&&i<N-1&&j<N-1){
                    Psi1[i][j] = Psi0[i][j] - dt/dy * (Vy0[i][j] - Vy0[i][j-1])*Csq(x,y)*Sigma_x(x);
                }

            }

        }
        std::swap(P0,P1);
        std::swap(Vx0,Vx1);
        std::swap(Vy0,Vy1);
        std::swap(Psi0,Psi1);
        output_layer(m, t, P1);
    }
    finished();

    return 0;

}

