#include "problem.hpp"
#include <cmath>
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
