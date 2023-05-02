#pragma once
#include <vector>
#include "problem.hpp"
class Solver{
public:
    Solver(int M,int N);
    grid * next_frame(double &t);
    virtual ~Solver();
private:
    int N;
    int M;
    double dx;
    double dy;
    double dt;
    int m=0;

    grid *P0;  
    grid *P1; 
    grid *Vx0; 
    grid *Vx1; 
    grid *Vy0; 
    grid *Vy1; 
    grid *Psi0;
    grid *Psi1;

};
