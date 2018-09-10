#include <stdio.h>
#include <stdlib.h>
using namespace std; 

const int Nx = 64;
const int Ny = 64;
const int Nz = 64;
const int Nx1 = Nx+1;
const int Ny1 = Ny+1;
const int Nz1 = Nz+1;
const int L = Nx + 1;
const int Q = 19;
const int dt = 1;

const double rho0 = 1.0;
const double ux0 = 0.0;// Initial conditions velocity
const double uy0 = 0.0;
const double uz0 = 0.0;
const double Jx0 = 0.0;// Initial conditions density current vector 
const double Jy0 = 0.0;
const double Jz0 = 0.0;
const double Bx0 = 0.0;// Initial conditions magnetic field
const double By0 = 0.0;
const double Bz0 = 0.0;

int cx[Q] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0};
int cy[Q] = {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1};
int cz[Q] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1};

double f[Nx1][Ny1][Nz1][Q]; 
double f_post[Nx1][Ny1][Nz1][Q]; 
double rho[Nx1][Ny1][Nz1];
double ux[Nx1][Ny1][Nz1];
double uy[Nx1][Ny1][Nz1];
double uz[Nx1][Ny1][Nz1];

double Jx[Nx1][Ny1][Nz1];
double Jy[Nx1][Ny1][Nz1];
double Jz[Nx1][Ny1][Nz1];

double Bx[Nx1][Ny1][Nz1];
double By[Nx1][Ny1][Nz1];
double Bz[Nx1][Ny1][Nz1];

double tau;


double w[Q] = {1.0/3 ,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,
			   1.0/18,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,
			   1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36};

void Init_Eq(void);
double feq(double RHO, double U, double V, double W, double J,double B,int q);
void Coll_BGK(void);
void Streaming(void); 
void Macro(void); 
//void BBOS(void);
double u0[Nx1][Ny1][Nz1],v0[Nx1][Ny1][Nz1],w0[Nx1][Ny1][Nz1];
void Data_Output(void);
double Fi(double RHO, double U, double V, double W,int q);
double af(double RHO, double U, double V, double W, double J,double B,int q);




int main(){
	printf("Done !\n");
	return 0;
}


double feq(double RHO, double U, double V, double W,int q)
{
	double cu, U2;
	cu=cx[q]*U+cy[q]*V+cz[q]*W; 
	U2=U*U+V*V+W*W; 
	return w[q]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2);
}


double af(double RHO, double U, double V, double W, double Jx,double Jy, double Jz, double Bx, double By, double Bz, int q)
{
	double cu;
	cu = cx[q]*U + cy[q]*V + cz[q]*W;
	return 3*w[q]*((cu*cx[q]-U)*(Jy*Bz-Jz*By) + (cu*cy[q]-V)*(Jz*Bx-Jx*Bz) + 
		(cu*cz[q]-W)*(Jy*Bx-Jx*By));
}


void Coll_BGK()
{
	int j, i, k, q;
	double FEQ;
	for (i=0;i<=Nx1;i++) for(j=0;j<=Ny1;j++) for(k=0;k<=Nz1;k++) for(q=0;q<Q;q++)
	{
		FEQ=feq(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q); 
		f_post[i][j][k][q] = (1 - 0.5*(dt/tau))*f[i][j][k][q] + 0.5*(dt/tau)*FEQ
		+(0.5*dt)*af(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],Jx[i][j][k],Jy[i][j][k],Jz[i][j][k],Bx[i][j][k],By[i][j][k],Bz[i][j][k],q);
	}
}
