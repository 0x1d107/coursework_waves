#pragma once
#include <vector>
typedef std::vector<std::vector<double>> grid;
double Pt0(double x,double y);
double Vxt0(double x,double y);
double Vyt0(double x,double y);
double F(double x, double y, double t);
double Csq(double x,double y);
double Sigma_x(double x);
void init();
void output_layer(int m,double t, const grid &P1, const grid &Vx,const grid &Vy);
void finished();
extern int M;
extern int N;
extern double Time;
