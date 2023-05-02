#pragma once
#include <vector>
double Pt0(double x,double y);
double Vxt0(double x,double y);
double Vyt0(double x,double y);
double F(double x, double y, double t);
double Csq(double x,double y);
double Sigma_x(double x);
typedef std::vector<std::vector<double>> grid;
