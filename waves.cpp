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
typedef std::vector<std::vector<double>> grid;
int main(){
    int M = 90000; // Time 
    int N = 300; // X Y
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


    int dir_status = mkdir("data", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(dir_status != 0 && errno != EEXIST){
        cerr << "Failed to create data directory" <<endl;
        return 1;
    }


    for(int j=0;j<N;j++){
        double y = j*dy;
        for(int i = 0; i < N; i++){
            double x = i*dx;
            P0[i][j] = Pt0(x,y);
            Vx0[i][j] = Vxt0(x,y);
            Vy0[i][j] = Vyt0(x,y);
        }
    }
    int T=100;
    int print_T = 5*T;
    for(int m=0;m<M;m++){
        double t = m * dt;
        #pragma omp parallel for
        for(int j=1;j<N-1;j++){
            double y = j*dy;
            #pragma omp parallel for
            for(int i = 1; i < N-1; i++){
                double x = i*dx;
                P1[i][j] = P0[i][j] - Csq(x,y) * dt*(
						(Vx0[i][j]-Vx0[i-1][j])/dx
						+ (Vy0[i][j] - Vy0[i][j-1])/dy 
						) + dt*F(x,y,t) - dt * Sigma_x(x)*P0[i][j] + dt*Psi0[i][j];
                Vx1[i][j] = Vx0[i][j] - dt/dx *( P1[i+1][j] - P1[i][j]) - dt*Sigma_x(x)*Vx0[i][j];
                Vy1[i][j] = Vy0[i][j] - dt/dy *( P1[i][j+1] - P1[i][j]);
                Psi1[i][j] = Psi0[i][j] - dt/dy * (Vy0[i][j] - Vy0[i][j-1])*Csq(x,y)*Sigma_x(x);

            }

        }
        std::swap(P0,P1);
        std::swap(Vx0,Vx1);
        std::swap(Vy0,Vy1);
        std::swap(Psi0,Psi1);
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

