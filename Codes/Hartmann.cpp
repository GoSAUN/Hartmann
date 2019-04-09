#include <stdio.h>
#include <stdlib.h>
using namespace std; 

const int Nx = 32;
const int Ny = 32;
const int Nz = 32;
const int Nx1 = Nx+1;
const int Ny1 = Ny+1;
const int Nz1 = Nz+1;
const int L = Nx + 1;
const int Q = 19;
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
const double Bx0 = 1.0;// Initial conditions magnetic field
const double By0 = 1.0;
const double Bz0 = 1.0;

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

bool EsFrontera[Nx1][Ny1][Nz1];
bool Inlet[Nx1][Ny1][Nz1];
bool Outlet[Nx1][Ny1][Nz1];

double tau;


double  w[Q] = {1.0/3 ,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,
			   1.0/18,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,
			   1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36};


double wm[M] = {1/4.0,1/8.0,1/8.0,1/8.0,1/8.0,1/8.0,1/8.0};


void Initialize(void);
double feq(double RHO, double U, double V, double W,double Bx, double By, double Bz,int q);
double geqx(double Bx, double By, double Bz, double U, double V, double W,int m);
double geqy(double Bx, double By, double Bz, double U, double V, double W,int m);
double geqz(double Bx, double By, double Bz, double U, double V, double W,int m);
void Coll_BGK(void);
void Coll_BGK_M(void);
void Streaming(void); 
void StreamingM(void); 
void Macro_Fluid(void); 
void BBOS(void);
void Data_Output(void);
double af(double RHO, double U, double V, double W, double J,double B,int q);



int main(){
	int n;
	Initialize();
	n=0;
	while(n<=100)
	{
		n++;
		Coll_BGK(); 
		Coll_BGK_M();
		Streaming(); 
		StreamingM(); 
		//BBOS();
		Macro_Fluid();
		printf("rho=%e ux_c=%e uy_c=%e uz_c=%e n=%d\n",
			rho[Nx/2][Ny/2][Nz/2],
			ux[Nx/2][Ny/2][Nz/2],
			uy[Nx/2][Ny/2][Nz/2],
			uz[Nx/2][Ny/2][Nz/2], n);
	}
	Data_Output();
	printf("Done !\n");
}


void Initialize()
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
		
		for(q=0;q<Q;q++){
			f[i][j][k][q]=feq(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],Bx[i][j][k],By[i][j][k],Bz[i][j][k],q);
		}

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

// Equilibrium functions 
double feq(double RHO, double U, double V, double W,int q)
{
	double cu, U2;
	cu=cx[q]*U+cy[q]*V+cz[q]*W; 
	U2=U*U+V*V+W*W; 
	return w[q]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2);
}


double feq(double RHO, double U, double V, double W,double Bx, double By, double Bz,int q)
{
	double cu, U2,c2,B2,cB;
	cu=cx[q]*U+cy[q]*V+cz[q]*W; 
	U2=U*U+V*V+W*W; 
	c2 = cx[q]*cx[q]+cy[q]*cy[q]+cz[q]*cz[q];
	B2 = Bx*Bx + By*By + Bz*Bz; 
	cB = cx[q]*Bx + cy[q]*By +cz[q]*Bz; 
	return w[q]*RHO*(1.0+3.0*cu+4.5*cu*cu-1.5*U2)-4.5*w[q]*(0.5*B2*c2-cB*cB);
}



double geqx(double Bx, double By, double Bz, double U, double V, double W,int m)
{
	return wm[m]*(Bx + 4.0*cmx[m]*(U*Bx-Bx*U + V*Bx-By*U + W*Bx - Bz*U));
}

double geqy(double Bx, double By, double Bz, double U, double V, double W,int m)
{
	return wm[m]*(By + 4.0*cmy[m]*(U*By-Bx*V + V*By-By*V + W*By - Bz*V));
}

double geqz(double Bx, double By, double Bz, double U, double V, double W,int m)
{
	return wm[m]*(Bz + 4.0*cmz[m]*(U*Bz-Bx*W + V*Bz-Bz*V + W*By - Bz*W));
}


// Electrodynamics force
double af(double RHO, double U, double V, double W, 
		  double Jx,double Jy, double Jz, double 
		  Bx, double By, double Bz, int q)
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
		FEQ=feq(rho[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],Bx[i][j][k],By[i][j][k],Bz[i][j][k],q); 
		f_post[i][j][k][q] = (1 + 0.5*(dt/tau))*f[i][j][k][q] - 0.5*(dt/tau)*FEQ
		+(0.5*dt)*af(rho[i][j][k],
					 ux[i][j][k],uy[i][j][k],uz[i][j][k],
					 Jx[i][j][k],Jy[i][j][k],Jz[i][j][k],
					 Bx[i][j][k],By[i][j][k],Bz[i][j][k],q);
	}
}



void Coll_BGK_M()
{
	int j, i, k, m;
	double GEQX,GEQY,GEQZ;
	for (i=0;i<=Nx1;i++) for(j=0;j<=Ny1;j++) for(k=0;k<=Nz1;k++) for(m=0;m<M;m++)
	{
		GEQX=geqx(Bx[i][j][k],By[i][j][k],Bz[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],m); 
		GEQY=geqy(Bx[i][j][k],By[i][j][k],Bz[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],m); 
		GEQZ=geqz(Bx[i][j][k],By[i][j][k],Bz[i][j][k],ux[i][j][k],uy[i][j][k],uz[i][j][k],m); 
		g_post[0][i][j][k][m] = (1 + 0.5*(dt/tau))*g[0][i][j][k][m] - 0.5*(dt/tau)*GEQX;
		g_post[1][i][j][k][m] = (1 + 0.5*(dt/tau))*g[1][i][j][k][m] - 0.5*(dt/tau)*GEQY;
		g_post[2][i][j][k][m] = (1 + 0.5*(dt/tau))*g[2][i][j][k][m] - 0.5*(dt/tau)*GEQZ;
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



void BBOS(){
	int i,j,k,q;
	for (i=0;i<=Nx;i++) for(j=0;j<=Ny;j++) for(k=0;k<=Nz;k++){ 
		
		if(Inlet[i][j][k] == true){
			//plano arriba 
			f[i][j][k][4]=f_post[i-cx[4]*dt][0][k-cz[4]*dt][4];
			f[i][j][k][8]=f_post[i-cx[8]*dt][0][k-cz[8]*dt][8];
			f[i][j][k][12]=f_post[i-cx[12]*dt][0][k-cz[12]*dt][12];
			f[i][j][k][13]=f_post[i-cx[13]*dt][0][k-cz[13]*dt][13];
			f[i][j][k][18]=f_post[i-cx[18]*dt][0][k-cz[18]*dt][18];
		}

		if (Outlet[i][j][k] == true){
			//plano abajo
			f[i][j][k][3]=f_post[i-cx[3]*dt][Ny][k-cz[3]*dt][3];
			f[i][j][k][7]=f_post[i-cx[7]*dt][Ny][k-cz[7]*dt][7];
			f[i][j][k][11]=f_post[i-cx[11]*dt][Ny][k-cz[11]*dt][11];
			f[i][j][k][14]=f_post[i-cx[14]*dt][Ny][k-cz[14]*dt][14];
			f[i][j][k][17]=f_post[i-cx[17]*dt][Ny][k-cz[17]*dt][17];
			}
		
		if(EsFrontera[i][j][k]==true){
			//for (int q = 1; q < Q; q++)
				//f[i][j][k][q]=f_post[i][j][k][int(q+pow(-1,q))];
				
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






void Macro_Fluid()
{
	int i,j,k;
	for(i=0;i <= Nx;i++)
		for(i=0;i <= Nx;i++)
			for(i=0;i <= Nx;i++)
			{
				rho[i][j][k] = f[i][j][k][0]+f[i][j][k][1]+f[i][j][k][2]+f[i][j][k][3]+
				f[i][j][k][4]+f[i][j][k][5]+f[i][j][k][6]+f[i][j][k][7]+f[i][j][k][8]+
				f[i][j][k][9]+f[i][j][k][10]+f[i][j][k][11]+f[i][j][k][12]+f[i][j][k][13]+
				f[i][j][k][14]+f[i][j][k][15]+f[i][j][k][16]+f[i][j][k][17]+f[i][j][k][18];

				ux[i][j][k] = (f[i][j][k][1]+f[i][j][k][7]+f[i][j][k][9]+f[i][j][k][13]+f[i][j][k][15]
			    -f[i][j][k][2]-f[i][j][k][8]-f[i][j][k][10]-f[i][j][k][14]-f[i][j][k][16]
			    +(Jx[i][j][k]*Bz[i][j][k] - Jz[i][j][k]*Bx[i][j][k])*dt*0.5 )/rho[i][j][k];

				uy[i][j][k] = (f[i][j][k][3] + f[i][j][k][7] + f[i][j][k][11] + 
				f[i][j][k][14] + f[i][j][k][17] - f[i][j][k][4] - f[i][j][k][8]
			    -f[i][j][k][12] - f[i][j][k][13] - f[i][j][k][18]
			    +(Jz[i][j][k]*Bx[i][j][k] - Jx[i][j][k]*Bz[i][j][k])*dt*0.5
			    )/rho[i][j][k];

				uz[i][j][k] = (f[i][j][k][5] + f[i][j][k][9] + f[i][j][k][11]
				+f[i][j][k][16] + f[i][j][k][18] - f[i][j][k][6] - f[i][j][k][10]
				-f[i][j][k][12] - f[i][j][k][15] -f[i][j][k][17]
				+(Jx[i][j][k]*By[i][j][k] - Jy[i][j][k]*Bx[i][j][k])*dt*0.5
				)/rho[i][j][k];

				Bx[i][j][k] = g[0][i][j][k][0]+g[0][i][j][k][1]+g[0][i][j][k][2]+g[0][i][j][k][3]
				+g[0][i][j][k][4]+g[0][i][j][k][5]+g[0][i][j][k][6];

				By[i][j][k] = g[1][i][j][k][0]+g[1][i][j][k][1]+g[1][i][j][k][2]+g[1][i][j][k][3]
				+g[1][i][j][k][4]+g[1][i][j][k][5]+g[1][i][j][k][6];

				Bz[i][j][k] = g[2][i][j][k][0]+g[2][i][j][k][1]+g[2][i][j][k][2]+g[2][i][j][k][3]
				+g[2][i][j][k][4]+g[2][i][j][k][5]+g[2][i][j][k][6];
				// curl of magnetic field


			}

}


void Data_Output() 
{
	int i,j,k,z;
	z = 1;
	FILE *fp;
	fp=fopen("x.dat","w+");
	for(i=0;i<=Nx;i++) fprintf(fp,"%e \n", float(i)/L);
	fclose(fp);
	fp=fopen("y.dat","w+");
	for(j=0;j<=Ny;j++) fprintf(fp,"%e \n", float(j)/L);
	fclose(fp);
	fp=fopen("vx.dat","w");
	for(i=0;i<=Nx;i++) {
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

	fp=fopen("rho.dat","w");
	for(i=0;i<=Nx;i++){
	for (j=0; j<=Ny; j++) fprintf(fp,"%e ",rho[j][i][z]);
	fprintf(fp,"\n");
	}
	fclose(fp);
}






