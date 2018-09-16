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
const int Q = 19
const int M = 7;
const int Qm = 7;
const int dt = 1;
const int a = 3;

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

int cmx[Qm] = {0, 1, 0,-1, 0, 0, 0};
int cmy[Qm] = {0, 0, 0, 0, 0, 1,-1};
int cmz[Qm] = {0, 0, 1, 0,-1, 0, 0};

double f[Nx1][Ny1][Nz1][Q]; 
double f_post[Nx1][Ny1][Nz1][Q];
double g[a][Nx1][Ny1][Nz1][Q]; 
double g_post[a][Nx1][Ny1][Nz1][Q]; 
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


double W[M] = {1/4.0,1/8.0,1/8.0,1/8.0,1/8.0,1/8.0,1/8.0};


void Init(void);
double feq(double RHO, double U, double V, double W,int q);
double geqx(double Bx, double By, double Bz, double U, double V, double W,int m);
double geqy(double Bx, double By, double Bz, double U, double V, double W,int m);
double geqz(double Bx, double By, double Bz, double U, double V, double W,int m);
void Coll_BGK(void);
void Streaming(void); 
void StreamingM(void); 
void Macro(void); 
void BBOS(void);
void Data_Output(void);
double af(double RHO, double U, double V, double W, double J,double B,int q);



int main(){
	printf("Done !\n");
	return 0;
}


void Init_Eq()
{
	int j, i, k, q;
	for(i=0;i<=Nx;i++) for(j=0;j<=Ny;j++) for(k = 0; k<=Nz; k++)
	{
		rho[i][j][k] = rho0;
		ux[i][j][k] = ux0;
		uy[i][j][k] = uy0;
		uz[i][j][k] = uz0;
		Jx[i][j][k] = ux0;
		Jy[i][j][k] = uy0;
		Jz[i][j][k] = uz0;
		Bx[i][j][k] = ux0;
		By[i][j][k] = uy0;
		Bz[i][j][k] = uz0;
		
		for(q=0;q<Q;q++){
			f[i][j][k][q]=feq(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q);
		}

		for(m = 0; m<M;m++){
			g[0][i][j][k][m] =geqx(double Bx[i][j][k], double By[i][j][k], double Bz[i][j][k], ux[i][j][k],uy[i][j][k],uz[i][j][k],m);
			g[1][i][j][k][m] =geqy(double Bx[i][j][k], double By[i][j][k], double Bz[i][j][k], ux[i][j][k],uy[i][j][k],uz[i][j][k],m);
			g[2][i][j][k][m] =geqz(double Bx[i][j][k], double By[i][j][k], double Bz[i][j][k], ux[i][j][k],uy[i][j][k],uz[i][j][k],m);
		}

	}
}

// Funciones de equilibrio 
double feq(double RHO, double U, double V, double W,int q)
{
	double cu, U2;
	cu=cx[q]*U+cy[q]*V+cz[q]*W; 
	U2=U*U+V*V+W*W; 
	return w[q]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2);
}


double geqx(double Bx, double By, double Bz, double U, double V, double W,int m)
{
	return W[m]*(Bx + 4.0*cmx[m]*(U*Bx-Bx*U + V*Bx-By*U + W*Bx - Bz*U));
}

double geqy(double Bx, double By, double Bz, double U, double V, double W,int m)
{
	return W[m]*(By + 4.0*cmy[m]*(U*By-Bx*V + V*By-By*V + W*By - Bz*V));
}

double geqz(double Bx, double By, double Bz, double U, double V, double W,int m)
{
	return W[m]*(Bz + 4.0*cmz[m]*(U*Bz-Bx*W + V*Bz-Bz*V + W*By - Bz*W));
}


// Termino de forzamiento electrodinamico
double af(double RHO, double U, double V, double W, double Jx,double Jy, double Jz, double Bx, double By, double Bz, int q)
{
	double cu;
	cu = cx[q]*U + cy[q]*V + cz[q]*W;
	return 3*w[q]*RHO*((cu*cx[q]-U)*(Jy*Bz-Jz*By) + (cu*cy[q]-V)*(Jz*Bx-Jx*Bz) + 
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

void Streaming()
{
	int j, i, k, jd, id, kd, q;
	for (i=0;i<Nx1;i++) for(j=0;j<Ny1;j++) for(k=0;k<Nz1;k++) for(q=0;q<Q;q++){
	jd=j-cy[q]; id=i-cx[q]; kd=k-cz[q]; 
	if(jd>=0 && jd<=Ny && id>=0 && id<=Nx && kd>=0 && kd<=Nz){
		f[i][j][k][q]=f_post[id][jd][kd][q]; // streaming
		}
	}
}

void StreamingM()
{
	int j, i, k, jd, id, kd, m;
	for (i=0;i<Nx1;i++) for(j=0;j<Ny1;j++) for(k=0;k<Nz1;k++) for(m=0;m<M;m++){
	jd=j-cy[m]; id=i-cx[m]; kd=k-cz[m]; 
	if(jd>=0 && jd<=Ny && id>=0 && id<=Nx && kd>=0 && kd<=Nz){
		g[0][i][j][k][m]=g_post[1][id][jd][kd][m]; // streaming
		g[1][i][j][k][m]=g_post[1][id][jd][kd][m]; // streaming
		g[2][i][j][k][m]=g_post[2][id][jd][kd][m]; // streaming
		}
	}
}
