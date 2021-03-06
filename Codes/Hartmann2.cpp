#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Nx 64
#define Ny 64
#define Nz 64
#define Nx1 (Nx+1)
#define Ny1 (Ny+1)
#define Nz1 (Nz+1)
#define dx 1/Nx
#define dy 1/Ny
#define dz 1/Nz
#define L (Ny+1)
#define Q 19 		
#define rho0 1.0  
#define ux0 0.0
#define uy0 0.0   
#define uz0 0.0   
#define go 4.06901041667e-7
#define dt 1
#define M 7
#define Qm 7
#define a 3

const double Jx0 = 0.0;// Initial conditions density current vector 
const double Jy0 = 0.0;
const double Jz0 = 0.0;
double Bx0;// Initial conditions magnetic field
const double By0 = 0.0;
const double Bz0 = 0.0;



int cx[Q] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0};
int cy[Q] = {0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1};
int cz[Q] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1};
//			 0  1  2  3. 4. 5. 6. 7. 8. 9. 10 11 12 13 14 15 16 17 18

int cmx[Qm] = {0, 1, 0,-1, 0, 0, 0};
int cmy[Qm] = {0, 0, 0, 0, 0,-1, 1};
int cmz[Qm] = {0, 0, 1, 0,-1, 0, 0};
//             0. 1. 2. 3. 4. 5. 6.   

double Fx, Fy, Fz;
double f[Nx1][Ny1][Nz1][Q]; 
double f_post[Nx1][Ny1][Nz1][Q]; 
double g[a][Nx1][Ny1][Nz1][Q]; 
double g_post[a][Nx1][Ny1][Nz1][Q]; 
double rho[Nx1][Ny1][Nz1], ux[Nx1][Ny1][Nz1], uy[Nx1][Ny1][Nz1],uz[Nx1][Ny1][Nz1];

double Jx[Nx1][Ny1][Nz1];
double Jy[Nx1][Ny1][Nz1];
double Jz[Nx1][Ny1][Nz1];

double Bx[Nx1][Ny1][Nz1];
double By[Nx1][Ny1][Nz1];
double Bz[Nx1][Ny1][Nz1];
double Btest[Nx1][Ny1][Nz1];

bool EsFrontera[Nx1][Ny1][Nz1];
bool Inlety[Nx1][Ny1][Nz1];
bool Inletz[Nx1][Ny1][Nz1];
bool Outlety[Nx1][Ny1][Nz1];
bool Outletz[Nx1][Ny1][Nz1];

double tau; 
double taum;
double w[Q] = {1.0/3 ,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,
			   1.0/18,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,
			   1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36};
			//   0      1       2        3      4       5      
  			//   6      7       8        9     10      11      
			//   12     13     14       15     16      17  	18 

double wm[M] = {1/4.0,1/8.0,1/8.0,1/8.0,1/8.0,1/8.0,1/8.0};
//               0.     1.    2.    3.    4.    5.    6
void Init_Eq(void);
double feq(double RHO, double U, double V, double W,int q);
double geqx(double Bx, double By, double Bz, double U, double V, double W,int m);
double geqy(double Bx, double By, double Bz, double U, double V, double W,int m);
double geqz(double Bx, double By, double Bz, double U, double V, double W,int m);

void Coll_BGK(void);
void Coll_BGKP(void);
void Coll_BGK_M(void);

void Streaming(void); 
void StreamingM(void); 

void Den_Vel_Mag(void); 

void BBOS(void);
void BCMagnetic(void);

double u0[Nx1][Ny1][Nz1],v0[Nx1][Ny1][Nz1],w0[Nx1][Ny1][Nz1];
void Data_Output(void);
double af(double RHO, double U, double V, double W, double J,double B,int q);
double Si(double RHO, double U, double V, double W,int q);
void divB();



int main(int argc, char* argv[])
{
	int n,M2,N2,O2;
	double uold,error, tol= 1e-9;
	double Bx0;
	Bx0 = atoi(argv[1]);
	M2=Nx/2; N2=Ny/2; O2 = Nz/2;
	n=0;
	tau=0.6;
	taum = 0.6;
	uold = 0.0;
	error = 1e-1;
	printf("tolerancia=%e\n",tol);
	Init_Eq();
	while(error > tol)
	{
		n++;
		uold = uy[M2][N2][O2];
		Coll_BGK(); 
		Coll_BGK_M();
		Streaming(); 
		StreamingM(); 
		BBOS(); 
		BCMagnetic();
		Den_Vel_Mag();
		divB();
		error = abs(uy[M2][N2][O2]-uold);
		//printf("rho=%e ux_c=%e uy_c=%e uz_c=%e k=%d\n",
		//	rho[M2][N2][O2],ux[M2][N2][O2],uy[M2][N2][O2],
		//	uz[M2][N2][O2], n); 
		if (n%500==0)
			{
				printf("div B =%e uy=%.10e error=%.10e k=%d\n",Btest[M2][N2][O2],uy[M2][N2][O2],error, n); 
				Data_Output();
			}	
	}
	printf("%d\n", n);
	printf("Done !\n"); 
}



void Init_Eq()
{
	int j, i, k, q, m;
	for(i=0;i<=Nx;i++) for(j=0;j<=Ny;j++) for(k = 0; k<=Nz; k++)
	{
		rho[i][j][k] = rho0;
		ux[i][j][k] = ux0;
		uy[i][j][k] = uy0;
		uz[i][j][k] = uz0;

		Jx[i][j][k] = Jx0;
		Jy[i][j][k] = Jy0;
		Jz[i][j][k] = Jz0;
		Bx[i][j][k] = Bx0;
		By[i][j][k] = By0;
		Bz[i][j][k] = Bz0;

		EsFrontera[i][j][k] = false;
		EsFrontera[0][j][k] = true;
		EsFrontera[Nx][j][k] = true;
		Inlety[i][Ny][k] = true;
		Outlety[i][0][k] = true;
		Inletz[i][j][Nz] = true;
		Outletz[i][j][0] = true;
		
		for(q=0;q<Q;q++)
		f[i][j][k][q]=feq(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q);
	
		for(m = 0; m<M;m++){
			g[0][i][j][k][m] =geqx(Bx[i][j][k], 
								   By[i][j][k], 
								   Bz[i][j][k], 
								   ux[i][j][k],
								   uy[i][j][k],
								   uz[i][j][k],m);
			g[1][i][j][k][m] =geqy(Bx[i][j][k], 
								   By[i][j][k], 
								   Bz[i][j][k], 
								   ux[i][j][k],
								   uy[i][j][k],
								   uz[i][j][k],m);
			g[2][i][j][k][m] =geqz(Bx[i][j][k], 
								   By[i][j][k], 
								   Bz[i][j][k], 
								   ux[i][j][k],
								   uy[i][j][k],
								   uz[i][j][k],m);
		}

	}
}
 
double feq(double RHO, double U, double V, double W,int q)
{
	double cu, U2;
	cu=cx[q]*U + cy[q]*V + cz[q]*W; 
	U2=U*U + V*V + W*W; 
	return w[q]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2);
}


double geqx(double Bx, double By, double Bz, double U, double V, double W,int m)
{
	return wm[m]*(Bx + 4.0*cmx[m]*(U*Bx-Bx*U + V*Bx-By*U + W*Bx-Bz*U));
}

double geqy(double Bx, double By, double Bz, double U, double V, double W,int m)
{
	return wm[m]*(By + 4.0*cmy[m]*(U*By-Bx*V + V*By-By*V + W*By-Bz*V));
}

double geqz(double Bx, double By, double Bz, double U, double V, double W,int m)
{
	return wm[m]*(Bz + 4.0*cmz[m]*(U*Bz-Bx*W + V*Bz-By*W + W*Bz - Bz*W));
}



// Electrodynamics force
double af(double RHO, double U, double V, double W, 
		  double Jx,double Jy, double Jz, double 
		  Bx, double By, double Bz, int q)
{
	double Me;
	Me = 9.0*(cx[q]*U + cy[q]*V + cz[q]*W + 1);
	return -9.0*w[q]*RHO*((M*cx[q]-U)*(Jx*Bz-Jz*By)+(M*cy[q]-V)*(Jz*Bx-Jx*Bz)+(M*cz[q]-W)*(Jy*Bx-Jx*By));
}



double Si(double RHO, double U, double V, double W,int q)
{
	double Fx = 0.0, Fy = -RHO*go, Fz = 0.0;
	double t1 = cx[q]*Fx + cy[q]*Fy + cz[q]*Fz;
	double t2 = cx[q]*cx[q]*U*Fx + cx[q]*cy[q]*(V*Fx + U*Fy) + 
				cx[q]*cz[q]*(W*Fx + U*Fz) + cy[q]*cy[q]*V*Fy +
				cy[q]*cz[q]*(V*Fz + W*Fy) + cz[q]*cz[q]*W*Fz;
	double t3 = U*Fx + V*Fy + W*Fz;
	return (1-((dt)/(2*tau)))*w[q]*(3.0*t1+ 9.0*t2- 3.0*t3);	
}

void Coll_BGKP()
{
	int j, i, k, q;
	double FEQ;
	for (i=0;i<=Nx1;i++) for(j=0;j<=Ny1;j++) for(k=0;k<=Nz1;k++) for(q=0;q<Q;q++)
	{
		FEQ=feq(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q); 
		f_post[i][j][k][q] = (1 - (dt/tau))*f[i][j][k][q]+ (dt/tau)*FEQ
		+dt*Si(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q);
	}
}

void Coll_BGK()
{
	int j, i, k, q;
	double FEQ,acel,P;
	double FBAR;
	for (i=0;i<=Nx1;i++) for(j=0;j<=Ny1;j++) for(k=0;k<=Nz1;k++) for(q=0;q<Q;q++)
	{
		FEQ=feq(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q);
		acel =  af(rho[i][j][k], ux[i][j][k], uy[i][j][k], uz[i][j][k], 
		  Jx[i][j][k],Jy[i][j][k], Jz[i][j][k], Bx[i][j][k], By[i][j][k], Bz[i][j][k], q);
		P = Si(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],q);
		FBAR = f[i][j][k][q] + (0.5*dt/tau)*(f[i][j][k][q]-FEQ)-(0.5*dt/tau)*acel;

		f_post[i][j][k][q] = FBAR - (dt/(tau+0.5*dt))*(FBAR-FEQ)+(tau/(tau+0.5*dt))*acel
		+dt*P;
	}
}


void Coll_BGK_M()
{
	int j, i, k, m;
	double GEQX,GEQY,GEQZ,GBARX,GBARY,GBARZ;
	for (i=0;i<=Nx1;i++) for(j=0;j<=Ny1;j++) for(k=0;k<=Nz1;k++) for(m=0;m<M;m++)
	{
		GEQX=geqx(Bx[i][j][k],By[i][j][k],Bz[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],m); 
		GEQY=geqy(Bx[i][j][k],By[i][j][k],Bz[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],m); 
		GEQZ=geqz(Bx[i][j][k],By[i][j][k],Bz[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],m);
		GBARX=g[0][i][j][k][m] + 0.5*(dt/taum)*(g[0][i][j][k][m]-GEQX);
		GBARY=g[1][i][j][k][m] + 0.5*(dt/taum)*(g[1][i][j][k][m]-GEQY);
		GBARZ=g[2][i][j][k][m] + 0.5*(dt/taum)*(g[2][i][j][k][m]-GEQZ);
		g_post[0][i][j][k][m] = GBARX - (dt/(taum+0.5*dt))*(GBARX-GEQX);
		g_post[1][i][j][k][m] = GBARY - (dt/(taum+0.5*dt))*(GBARY-GEQY);
		g_post[2][i][j][k][m] = GBARZ - (dt/(taum+0.5*dt))*(GBARZ-GEQZ);
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
		g[0][i][j][k][m]=g_post[0][id][jd][kd][m]; // streaming
		g[1][i][j][k][m]=g_post[1][id][jd][kd][m]; // streaming
		g[2][i][j][k][m]=g_post[2][id][jd][kd][m]; // streaming
		}
	}
}

void BBOS(){
	int i,j,k;
	for (i=0;i<=Nx;i++) for(j=0;j<=Ny;j++) for(k=0;k<=Nz;k++){ 
		
		if(Inlety[i][j][k] == true){
			//plano arriba 
			f[i][j][k][4]=f_post[i-cx[4]*dt][0][k-cz[4]*dt][4];
			f[i][j][k][8]=f_post[i-cx[8]*dt][0][k-cz[8]*dt][8];
			f[i][j][k][12]=f_post[i-cx[12]*dt][0][k-cz[12]*dt][12];
			f[i][j][k][13]=f_post[i-cx[13]*dt][0][k-cz[13]*dt][13];
			f[i][j][k][18]=f_post[i-cx[18]*dt][0][k-cz[18]*dt][18];
		}

		if (Outlety[i][j][k] == true){
			//plano abajo
			f[i][j][k][3]=f_post[i-cx[3]*dt][Ny][k-cz[3]*dt][3];
			f[i][j][k][7]=f_post[i-cx[7]*dt][Ny][k-cz[7]*dt][7];
			f[i][j][k][11]=f_post[i-cx[11]*dt][Ny][k-cz[11]*dt][11];
			f[i][j][k][14]=f_post[i-cx[14]*dt][Ny][k-cz[14]*dt][14];
			f[i][j][k][17]=f_post[i-cx[17]*dt][Ny][k-cz[17]*dt][17];
		}
		if(Inletz[i][j][k] == true){
			//plano adelante
			f[i][j][k][6]=f_post[i-cx[6]*dt][j-cy[6]*dt][Nz][6];
			f[i][j][k][10]=f_post[i-cx[10]*dt][j-cy[10]*dt][Nz][10];
			f[i][j][k][12]=f_post[i-cx[12]*dt][j-cy[12]*dt][Nz][12];
			f[i][j][k][15]=f_post[i-cx[15]*dt][j-cy[15]*dt][Nz][15];
			f[i][j][k][17]=f_post[i-cx[17]*dt][j-cy[17]*dt][Nz][17];
		}

		if (Outletz[i][j][k] == true){
			//plano atras
			f[i][j][k][3]=f_post[i-cx[3]*dt][j-cy[3]*dt][0][3];
			f[i][j][k][7]=f_post[i-cx[7]*dt][j-cy[7]*dt][Nz][7];
			f[i][j][k][11]=f_post[i-cx[11]*dt][j-cy[11]*dt][Nz][11];
			f[i][j][k][14]=f_post[i-cx[14]*dt][j-cy[14]*dt][Nz][14];
			f[i][j][k][17]=f_post[i-cx[17]*dt][j-cy[17]*dt][Nz][17];
		}
		
		if(EsFrontera[i][j][k]==true){
			//for (int q = 1; q < Q; q++)
				f[i][j][k][1]=f_post[i][j][k][2];
				f[i][j][k][2]=f_post[i][j][k][1];
				f[i][j][k][3]=f_post[i][j][k][4];
				f[i][j][k][4]=f_post[i][j][k][3];
				f[i][j][k][5]=f_post[i][j][k][6];
				f[i][j][k][6]=f_post[i][j][k][5];
				f[i][j][k][7]=f_post[i][j][k][8];
				f[i][j][k][8]=f_post[i][j][k][7];
				f[i][j][k][9]=f_post[i][j][k][10];
				f[i][j][k][10]=f_post[i][j][k][9];
				f[i][j][k][11]=f_post[i][j][k][12];
				f[i][j][k][12]=f_post[i][j][k][11];
				f[i][j][k][13]=f_post[i][j][k][14];
				f[i][j][k][14]=f_post[i][j][k][13];
				f[i][j][k][15]=f_post[i][j][k][16];
				f[i][j][k][16]=f_post[i][j][k][15];
				f[i][j][k][17]=f_post[i][j][k][18];
				f[i][j][k][18]=f_post[i][j][k][17];		
		}
	}
}

void BCMagnetic(){
	int i,j,k,m;
	for (i=0;i<=Nx;i++) for(j=0;j<=Ny;j++) for(k=0;k<=Nz;k++){

		if(Inlety[i][j][k] == true){
			//plano arriba 
			g[0][i][j][k][5] = g_post[0][i-cmx[5]*dt][0][k-cmz[5]*dt][5];
			g[1][i][j][k][5] = g_post[1][i-cmx[5]*dt][0][k-cmz[5]*dt][5];
			g[2][i][j][k][5] = g_post[2][i-cmx[5]*dt][0][k-cmz[5]*dt][5];
		}

		if (Outlety[i][j][k] == true){
			//plano abajo

			g[0][i][j][k][6] = g_post[0][i-cmx[6]*dt][Ny][k-cmz[6]*dt][6];
			g[1][i][j][k][6] = g_post[0][i-cmx[6]*dt][Ny][k-cmz[6]*dt][6];
			g[2][i][j][k][6] = g_post[2][i-cmx[6]*dt][Ny][k-cmz[6]*dt][6];

		}
		
		if (Inletz[i][j][k] == true){
			//plano adelante 
			g[0][i][j][k][4] = g_post[0][i-cmx[4]*dt][j-cmy[4]*dt][Nz][4];
			g[1][i][j][k][4] = g_post[1][i-cmx[4]*dt][j-cmy[4]*dt][Nz][4];
			g[2][i][j][k][4] = g_post[2][i-cmx[4]*dt][j-cmy[4]*dt][Nz][4];
		}
		
		if (Outletz[i][j][k] == true){
			//plano atras
			g[0][i][j][k][2] = g_post[0][i-cmx[2]*dt][j-cmy[2]*dt][0][2];
			g[1][i][j][k][2] = g_post[1][i-cmx[2]*dt][j-cmy[2]*dt][0][2];
			g[2][i][j][k][2] = g_post[2][i-cmx[2]*dt][j-cmy[2]*dt][0][2];
		}
		
		if(EsFrontera[i][j][k]==true){
			for(m=0;m<M;m++)
			{
				
					g[0][i][j][k][m] =geqx(0.0,0.0,0.0,0.0,0.0,0.0,m);
					g[1][i][j][k][m] =geqy(0.0,0.0,0.0,0.0,0.0,0.0,m);
					g[2][i][j][k][m] =geqz(0.0,0.0,0.0,0.0,0.0,0.0,m);
			}
		}							
	}
}



void Den_Vel_Mag()
{
	int j, i, k;
	double Fx = 0.0;
	double Fy = -go;
	double Fz = 0.0;
	for(i = 0; i <= Nx; i++) for(j = 0; j <= Ny; j++) for(k = 0; k <= Nz; k++)
		{
	
			rho[i][j][k] = f[i][j][k][0]+f[i][j][k][1]+f[i][j][k][2]+f[i][j][k][3]+
			f[i][j][k][4]+f[i][j][k][5]+f[i][j][k][6]+f[i][j][k][7]+f[i][j][k][8]+
			f[i][j][k][9]+f[i][j][k][10]+f[i][j][k][11]+f[i][j][k][12]+f[i][j][k][13]+
			f[i][j][k][14]+f[i][j][k][15]+f[i][j][k][16]+f[i][j][k][17]+f[i][j][k][18];

			ux[i][j][k] = (f[i][j][k][1]+f[i][j][k][7]+f[i][j][k][9]+f[i][j][k][13]+f[i][j][k][15]
			    -f[i][j][k][2]-f[i][j][k][8]-f[i][j][k][10]-f[i][j][k][14]-f[i][j][k][16])/rho[i][j][k]
				+ 0.5*dt*(Jy[i][j][k]*Bz[i][j][k]-Jz[i][j][k]*By[i][j][k])+0.5*dt*Fx;

			uy[i][j][k] = (f[i][j][k][3] + f[i][j][k][7] + f[i][j][k][11] + 
				f[i][j][k][14] + f[i][j][k][17] - f[i][j][k][4] - f[i][j][k][8]
				-f[i][j][k][12] - f[i][j][k][13] - f[i][j][k][18])/rho[i][j][k]
				+ 0.5*dt*(Jz[i][j][k]*Bx[i][j][k]-Jx[i][j][k]*Bz[i][j][k])+0.5*dt*Fy;

			uz[i][j][k] = (f[i][j][k][5] + f[i][j][k][9] + f[i][j][k][11]
				+f[i][j][k][16] + f[i][j][k][18] - f[i][j][k][6] - f[i][j][k][10]
				-f[i][j][k][12] - f[i][j][k][15] -f[i][j][k][17])/rho[i][j][k]
				+ 0.5*dt*(Jx[i][j][k]*By[i][j][k]-Jy[i][j][k]*Bx[i][j][k])+0.5*dt*Fz;

			Bx[i][j][k] = g[0][i][j][k][0]+g[0][i][j][k][1]+g[0][i][j][k][2]+g[0][i][j][k][3]
			+g[0][i][j][k][4]+g[0][i][j][k][5]+g[0][i][j][k][6];

			By[i][j][k] = g[1][i][j][k][0]+g[1][i][j][k][1]+g[1][i][j][k][2]+g[1][i][j][k][3]
			+g[1][i][j][k][4]+g[1][i][j][k][5]+g[1][i][j][k][6];
				
			Bz[i][j][k] = g[2][i][j][k][0]+g[2][i][j][k][1]+g[2][i][j][k][2]+g[2][i][j][k][3]
				+g[2][i][j][k][4]+g[2][i][j][k][5]+g[2][i][j][k][6];

				

			if ( i == 0)
			{
				Jx[i][j][k] = 0.5*(-3.0*Bz[i][j][k]+4.0*Bz[i][j+1][k]-Bz[i][j+2][k])/dy - 0.5*(-3.0*By[i][j][k]+4.0*By[i][j][k+1]-By[i][j][k+2])/dz ;			
			}
			else if (i == Nx)
			{
				Jx[i][j][k] = 0.5*(3.0*Bz[i][j][k]-4.0*Bz[i][j-1][k]+Bz[i][j-2][k])/dy - 0.5*(3.0*By[i][j][k]-4.0*By[i][j][k-1]+By[i][j][k-2])/dz ;
			}
			else
			{
				Jx[i][j][k] = (Bz[i][j+1][k]-Bz[i][j-1][k])/dy -(By[i][j][k+1]-By[i][j][k-1])/dz;
			}	



			if ( j == 0 )
			{
				Jy[i][j][k] = 0.5*(-3.0*Bx[i][j][k]+4.0*Bx[i][j][k+1]-Bx[i][j][k+2])/dz - 0.5*(-3.0*Bz[i][j][k]+4.0*Bz[i+1][j][k]-Bz[i+2][j][k])/dx ;			
			}
			else if ( j == Ny )
			{
				Jy[i][j][k] = 0.5*(3.0*Bx[i][j][k]-4.0*Bx[i][j][k-1]+Bx[i][j][k-2])/dz - 0.5*(3.0*Bz[i][j][k]-4.0*Bz[i+1][j][k]+Bz[i-2][j][k])/dx ;			
			}
			else
			{
				Jy[i][j][k] = (Bx[i][j][k+1]-Bx[i][j][k-1])/dz -(Bz[i+1][j][k]-Bz[i-1][j][k])/dx;
			}




			if ( k == 0 )
			{
				Jz[i][j][k] = 0.5*(-3.0*By[i][j][k]+4.0*By[i+1][j][k]-By[i+2][j][k])/dz - 0.5*(-3.0*Bx[i][j][k]+4.0*Bx[i][j+1][k]-Bx[i][j+2][k])/dx ;			
			}

			else if ( k == Nz )
			{
				Jz[i][j][k] = 0.5*(3.0*By[i][j][k]-4.0*By[i-1][j][k]+By[i-2][j][k])/dz - 0.5*(3.0*Bx[i][j][k]-4.0*Bx[i][j-1][k]+Bx[i][j-2][k])/dx ;			
			}
			else
			{
				Jz[i][j][k] = (By[i+1][j][k]-By[i-1][j][k])/dx -(Bx[i][j+1][k]-Bx[i][j-1][k])/dy;	
			}
			
		
		}
}

void divB()
{
	int i,j,k;
	for (i=0;i<=Nx;i++) for(j=0;j<=Ny;j++) for(k=0;k<=Nz;k++)
	{
		if (i == 0)
		{
			0.5*(-3.0*Bx[i][j][k]+4.0*Bx[i+1][j][k]-Bx[i+2][j][k])/dx;
		}
		else if (i == Nx)
		{
			0.5*(3.0*Bx[i][j][k]-4.0*Bx[i-1][j][k]+Bx[i-2][j][k])/dx;	
		}
		else if (j==0)
		{
			0.5*(-3.0*Bx[i][j][k]+4.0*Bx[i][j+1][k]+Bx[i][j+2][k])/dy;	
		}
		else if (j==Ny)
		{
			0.5*(3.0*Bx[i][j][k]-4.0*Bx[i][j-1][k]+Bx[i][j-2][k])/dy;
		}
		else if (k==0)
		{
			0.5*(-3.0*Bx[i][j][k]+4.0*Bx[i][j][k+1]-Bx[i][j][k+2])/dz;
		}
		else if (k==Nz)
		{
			0.5*(3.0*Bx[i][j][k]-4.0*Bx[i][j][k-1]+Bx[i][j][k-2])/dz;
		}
		else
		{
		Btest[i][j][k] = (Bx[i+1][j][k]-Bx[i-1][j][k])/dx + (By[i][j-1][k]-By[i][j-1][k])/dy +(Bz[i][j][k+1]-Bz[i][j][k-1])/dz;
		}
	}
}

void Data_Output() 
{
	int i,j,z,k;
	z = Nz/2;
	FILE *fp;
	fp=fopen("x.dat","w+");
	for(i=0;i<=Nx;i++) fprintf(fp,"%e \n", float(i)/L);
	fclose(fp);
	fp=fopen("y.dat","w+");
	for(j=0;j<=Ny;j++) fprintf(fp,"%e \n", float(j)/L);
	fclose(fp);
	fp=fopen("z.dat","w+");
	for(k=0;k<=Nz;k++) fprintf(fp,"%e \n", float(k)/L);
	fclose(fp);


	fp=fopen("vx.dat","w");
	for(i=0;i<=Nx;i++){
	for (j=0; j<=Ny; j++) fprintf(fp,"%e ",ux[i][j][z]);
	fprintf(fp,"\n");
	}
	fclose(fp);
	
	fp=fopen("vy.dat","w");
	for(i=0;i<=Nx;i++){
	for (j=0; j<=Ny; j++) fprintf(fp,"%e ",uy[i][j][z]);
	fprintf(fp,"\n");
	}
	fclose(fp);

	fp=fopen("vz.dat","w");
	for(i=0;i<=Nx;i++){
	for (j=0; j<=Ny; j++) fprintf(fp,"%e ",uz[i][j][z]);
	fprintf(fp,"\n");
	}
	fclose(fp);

	fp=fopen("rho.dat","w");
	for(i=0;i<=Nx;i++){
	for (j=0; j<=Ny; j++) fprintf(fp,"%e ",rho[j][i][z]);
	fprintf(fp,"\n");
	}
	fclose(fp);

	fp=fopen("Bx.dat","w");
	for(i=0;i<=Nx;i++){
	for (j=0; j<=Ny; j++) fprintf(fp,"%e ",Bx[i][j][z]);
	fprintf(fp,"\n");
	}
	fclose(fp);

	fp=fopen("By.dat","w");
	for(i=0;i<=Nx;i++){
	for (j=0; j<=Ny; j++) fprintf(fp,"%e ",By[i][j][z]);
	fprintf(fp,"\n");
	}
	fclose(fp);


	fp=fopen("Bz.dat","w");
	for(i=0;i<=Nx;i++){
	for (j=0; j<=Ny; j++) fprintf(fp,"%e ",Bz[i][j][z]);
	fprintf(fp,"\n");
	}
	fclose(fp);
}