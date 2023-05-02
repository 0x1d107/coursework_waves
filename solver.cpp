#include "solver.hpp"
#include <omp.h>
Solver::Solver(int M,int N){

    dt = 1.0/M;
    dx = 1.0/N;
    dy = dx;
    P0   = new grid(N,std::vector<double>(N,0));
    P1   = new grid(N,std::vector<double>(N,0));
    Vx0  = new grid(N,std::vector<double>(N,0));
    Vx1  = new grid(N,std::vector<double>(N,0));
    Vy0  = new grid(N,std::vector<double>(N,0));
    Vy1  = new grid(N,std::vector<double>(N,0));
    Psi0 = new grid(N,std::vector<double>(N,0));
    Psi1 = new grid(N,std::vector<double>(N,0));
    for(int j=0;j<N;j++){
        double y = j*dy;
        for(int i = 0; i < N; i++){
            double x = i*dx;
            (*P0)[i][j] = Pt0(x,y);
            (*Vx0)[i][j] = Vxt0(x,y);
            (*Vy0)[i][j] = Vyt0(x,y);
        }
    }

}
 grid *Solver::next_frame(double &t){
        t = m*dt;
        #pragma omp parallel for
        for(int j=1;j<N-1;j++){
            double y = j*dy;
            #pragma omp parallel for
            for(int i = 1; i < N-1; i++){
                double x = i*dx;
                (*P1)[i][j] = (*P0)[i][j] - Csq(x,y) * dt*(
						((*Vx0)[i][j]-(*Vx0)[i-1][j])/dx
						+ ((*Vy0)[i][j] - (*Vy0)[i][j-1])/dy) + dt*F(x,y,t) - dt * Sigma_x(x)*(*P0)[i][j] + dt*(*Psi0)[i][j];
                (*Vx1)[i][j] = (*Vx0)[i][j] - dt/dx *( (*P1)[i+1][j] - (*P1)[i][j]) - dt*Sigma_x(x)*(*Vx0)[i][j];
                (*Vy1)[i][j] = (*Vy0)[i][j] - dt/dy *( (*P1)[i][j+1] - (*P1)[i][j]);
                (*Psi1)[i][j] = (*Psi0)[i][j] - dt/dy * ((*Vy0)[i][j] - (*Vy0)[i][j-1])*Csq(x,y)*Sigma_x(x);

            }

        }
        grid *g;
        g = P0;
        P0 = P1;
        P1 = g;
        g = Vx0;
        Vx0 = Vx1;
        Vx1 = g;
        g = Vy0;
        Vy0 = Vy1;
        Vy1 = g;
        g = Psi0;
        Psi0 = Psi1;
        Psi1 = g;
        m++;
        return P1;
}
Solver::~Solver(){
    delete P0;
    delete P1;
    delete Vx0;
    delete Vx1;
    delete Vy0;
    delete Vy1;
    delete Psi0;
    delete Psi1;
}
