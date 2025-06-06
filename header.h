// Functions for jammed active particle distribution in 2D
// Note that PBC may not work properly for fluid systems.

using namespace std;
#include <math.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include "time.h"

typedef struct{
double x,y;
int t;
}vec;

//Parameters, variables and arrays
int steps;
int frame;
int nT;
double side;
double side_r;      //inverse of side
double m=1;
double damp=0.1;					           //Is this too much damping for inertial FIRE to work efficiently?
double deltaT=0.2;					   //deltaT=0.2 works good for MD
double totF_cutoff=1e-11;
double Gamma=exp(-damp*deltaT/m);
double c1=m/damp*(1-Gamma);
double c2=m/(damp*damp)*(damp*deltaT/m-1+Gamma);
double rA=1.3999999999999999/2.;
double rB=0.5;						
double rv,rv2,rc,rc2,Vgap,Vgap2;
vec position[20000];
vec velocity[20000];
vec forceInteraction[20000];
vec total_force[20000];
double totF;
double activeDirector[20000][2];
double activeSpeed;
double ke;
double ke_COM;                      
double pe;
double pe_Eff;
double fxCOM,fyCOM;
double Vx[20000];    // Positions of monomers in past
double Vy[20000];    // Positions of monomers in past
int Nverl[20000];    // The number of monomers in neighborhood of i^th monomer
int verl[20000][20]; // Second index identities monomers in neighborhood of i^th monomer in first index. There can't be more than about 12 monomers in neighborhood. 20 is a sufficiently big enough list. The first Nverl[i] entries in second index hold the identities of monomers of i^th monomer in first index. All succeeding entries in second index will remain zeroes.
int no_neighbors[20000];
double virial_pressure=0.0;
char * config;

//FIRE parameters
double alpha;
double alphaStart=0.1;
double fInc=1.1;
double deltaTmax=10*deltaT;
double fAlpha=0.99;
int Nmin=5;
int powerPositiveSteps;
double fDec=0.5;


//Function declarations
void mdrun();
void FIRE();
void updatePositionBrownian(int i); 
void brownianRun();
void forceInt(int i);
void totalForce(int i);
void averageTotalForce();
void updatePosition(int i);
void updateVelocity(int i);
void energy();
void print(int k);
void setVerlet();
void checkVerlet(int k);
void remove_rattlers();
void count_neighbors();
double calculate_deltaZ();
void write_config();
void printGapDist(int i);
void printHessian();
double doubleDerivPotDiag(double dx,double r,double d0);
double doubleDerivPotNonDiag(double dx,double dy,double r,double d0);
void dynamic_printHessian();

#include "./initialize.h" //initializations such as reading Hyderabad configs, LAMMPS dump config and generating random packed configurations

//Functions
void brownianRun(){

        for(int i=1; i<=nT; i++){
        totalForce(i);                                     //Calculate interaction+active force
        updatePositionBrownian(i);
        }
}

void updatePositionBrownian(int i){
//Position Update
position[i].x=position[i].x+ total_force[i].x*deltaT;
position[i].y=position[i].y+ total_force[i].y*deltaT;

checkVerlet(i);
}

void mdrun(){

	for(int i=1; i<=nT; i++){
	totalForce(i);					   //Calculate interaction+active force
	updatePosition(i);
	updateVelocity(i);
	}
}

void updatePosition(int i){
//Position Update
position[i].x=position[i].x+ c1*velocity[i].x+c2*total_force[i].x;
position[i].y=position[i].y+ c1*velocity[i].y+c2*total_force[i].y;

checkVerlet(i);
}

void updateVelocity(int i){
//Velocity Update
velocity[i].x=Gamma*velocity[i].x+1./damp*(1.-Gamma)*total_force[i].x;
velocity[i].y=Gamma*velocity[i].y+1./damp*(1.-Gamma)*total_force[i].y;
}

void FIRE(){
	for(int i=1; i<=nT; i++){
	updatePosition(i);
	totalForce(i);
	updateVelocity(i);
	//Check for convergence?
	}

//Calculate power
double power=0.;
	for(int i=1; i<=nT; i++){
	power += total_force[i].x*velocity[i].x;
	power += total_force[i].y*velocity[i].y;
	}

	if (power>0.) powerPositiveSteps += 1;
	else 	      powerPositiveSteps  = 0;


//Calculate unit vector of total force
double force_mag=0.;
	for(int i=1; i<=nT; i++){
	force_mag += total_force[i].x*total_force[i].x;
	force_mag += total_force[i].y*total_force[i].y;
	}
force_mag=sqrt(force_mag);

vec unit_total_force[20000];
	for(int i=1; i<=nT; i++){
	unit_total_force[i].x = total_force[i].x/force_mag;
	unit_total_force[i].y = total_force[i].y/force_mag;
	}

//Calculate magnitude of velocity
double modV=0.;
	for(int i=1; i<=nT; i++){
	modV += velocity[i].x*velocity[i].x;
	modV += velocity[i].y*velocity[i].y;
	}
modV=sqrt(modV);

//Rescale velocity
	for(int i=1; i<=nT; i++){
	velocity[i].x = (1-alpha)*velocity[i].x+alpha*unit_total_force[i].x*modV;
	velocity[i].y = (1-alpha)*velocity[i].y+alpha*unit_total_force[i].y*modV;
	}

	if (power>0. && powerPositiveSteps>Nmin){
	deltaT = min(deltaT*fInc,deltaTmax);
	alpha *= fAlpha;
	}

	if (power<=0.){
	deltaT *= fDec;
		for(int i=1; i<=nT; i++){
		velocity[i].x=0.;
		velocity[i].y=0.;
		}
	alpha=alphaStart;
	}
}

void forceInt(int i){
//Calculates interaction forces and pressure
//Note that pressure is calculated considring the system as quasi 2D, so d=3
double fxi,fyi;
double r,dx,dy;
int ti,tj;
double d0,ks=1.0;

	if (i==0) virial_pressure = 0.;

forceInteraction[i].x=0.0;
forceInteraction[i].y=0.0;

ti=position[i].t;

	for (int j=1;j<=Nverl[i];j++){     
	
	tj=position[verl[i][j]].t;
                        if      (ti+tj==2) d0=2.*rA;
                        else if (ti+tj==3) d0=rA+rB;
                        else if (ti+tj==4) d0=2.*rB;
	dx=position[i].x-position[ verl[i][j] ].x;
        dy=position[i].y-position[ verl[i][j] ].y;

	//PBC
	//dx -= side*nearbyint(dx*side_r);
	//dy -= side*nearbyint(dy*side_r);
		if (dx>side/2) dx-=side;
		if (dx<-side/2) dx+=side;
		if (dy>side/2) dy-=side;
		if (dy<-side/2) dy+=side;


        r=(dx*dx+dy*dy);
	r=sqrt(r);
		if (r<d0){
		//fxi= -ks*(r-d0)*(dx)/r;             //Check for 1/2
                //fyi= -ks*(r-d0)*(dy)/r;
		fxi= ks/d0*(1.-r/d0)*dx/r;             
                fyi= ks/d0*(1.-r/d0)*dy/r;

		forceInteraction[i].x += fxi;
		forceInteraction[i].y += fyi;

		virial_pressure += dx*fxi + dy*fyi;
		}
	}

	if (i==nT) virial_pressure /= 2.*2*side*side;  //Accounting for overcounting too
}

void totalForce(int i){
forceInt(i);  					   //Calculate interaction force
total_force[i].x=forceInteraction[i].x+activeSpeed*activeDirector[i][0];
total_force[i].y=forceInteraction[i].y+activeSpeed*activeDirector[i][1];
}

/*
void averageTotalForce(){
totF=0.;
double act_forcex=0.,act_forcey=0.;

        for (int i=1;i<=nT;i++){
	act_forcex += activeSpeed*activeDirector[i][0];
	act_forcey += activeSpeed*activeDirector[i][1];
	}

double f;
double fx,fy;
        for (int i=1;i<=nT;i++){
	fx=total_force[i].x-act_forcex;
	fy=total_force[i].y-act_forcey;
	//f=total_force[i].x*total_force[i].x+total_force[i].y*total_force[i].y;
	f=fx*fx+fy*fy;
	f=sqrt(f);
	totF += f;
	}
totF /= nT;
}
*/

void averageTotalForce(){

double f;
fxCOM=0.;fyCOM=0.;

        for (int i=1;i<=nT;i++){
	fxCOM+=total_force[i].x;
	fyCOM+=total_force[i].y;
	}

fxCOM /= nT;
fyCOM /= nT;

double fx,fy;
totF=0.;

        for (int i=1;i<=nT;i++){
	fx=total_force[i].x-fxCOM;
	fy=total_force[i].y-fyCOM;
	f=fx*fx+fy*fy;
	f=sqrt(f);
	totF += f;
	}

totF /= nT;
}

void energy(){
ke=0.;
ke_COM=0.;                      //KE in the COM frame
double vx_COM=0.,vy_COM=0.;
double dx,dy,r;
pe=0;
pe_Eff=0;

        for (int i=1;i<=nT;i++){
	vx_COM += velocity[i].x;
	vy_COM += velocity[i].y;
	}
vx_COM/= nT;
vy_COM/= nT;

        for (int i=1;i<=nT;i++){
	ke += velocity[i].x*velocity[i].x+velocity[i].y*velocity[i].y;
	ke_COM += (velocity[i].x-vx_COM)*(velocity[i].x-vx_COM)+(velocity[i].y-vy_COM)*(velocity[i].y-vy_COM);
	}
ke/=(2.*nT);
ke_COM/=(2.*nT);

int ti,tj;
double d0,ks=1.0;
double sum_act_vec=0.;
pe_Eff=0.;

        for (int i=1;i<=nT;i++){
	//pe_Eff += fxCOM*position[i].x;
	//pe_Eff += fyCOM*position[i].y;
	ti=position[i].t;
	sum_act_vec += activeDirector[i][0]*position[i].x+activeDirector[i][1]*position[i].y;
	//Use PBC position or unwrapped position? 

		for (int j=i+1;j<=nT;j++){
		tj=position[j].t;
                        if      (ti+tj==2) d0=2.*rA; 
                        else if (ti+tj==3) d0=rA+rB; 
                        else if (ti+tj==4) d0=2.*rB; 
		dx=position[i].x-position[j].x;
		dy=position[i].y-position[j].y;

		//PBC
		//dx -= side*round(dx*side_r);
		//dy -= side*round(dy*side_r);
			if (dx>side/2) dx-=side;
		        if (dx<-side/2) dx+=side;
			if (dy>side/2) dy-=side;
			if (dy<-side/2) dy+=side;


		r=(dx*dx+dy*dy);
		r=sqrt(r);
			if (r<d0){
			//pe += (r-d0)*(r-d0);                           
			pe += (1.-r/d0)*(1.-r/d0);                           
			}
		}
	}

pe *= 0.5*ks;
//pe_Eff *= nT;
pe_Eff +=  -activeSpeed*sum_act_vec;
//pe_Eff += -totF*sum_act_vec;
//pe/=(2*nT);   //why is there a two here? you already multiplied by 0.5 above and there is no overcounting in i,j for loop above
pe_Eff += pe;
pe_Eff /=nT;
pe /=nT;
}

void print(int k){
double x,y;
FILE *out;
out=fopen("out.dump","a");
fprintf(out,"ITEM: TIMESTEP\n");
fprintf(out,"%d\n",k);
fprintf(out,"ITEM: NUMBER OF ATOMS\n");
fprintf(out,"%d\n",nT);
fprintf(out,"ITEM: BOX BOUNDS pp pp pp\n");
fprintf(out,"-%.16f %.16f \n",side/2.,side/2.);
fprintf(out,"-%.16f %.16f \n",side/2.,side/2.);
fprintf(out,"-%.16f %.16f \n",0.7,0.7);
fprintf(out,"ITEM: ATOMS id type x y\n");

        for (int i=1;i<=nT;i++){
	//Printing PBC
	//x=position[i].x-floor(position[i].x/side)*side -side/2.;
	//y=position[i].y-floor(position[i].y/side)*side -side/2.;
	x=position[i].x;
	y=position[i].y;
	if (x>side/2.) x-= side;
	if (x<-side/2.) x+= side;
	if (y>side/2.) y-= side;
	if (y<-side/2.) y+= side;
	fprintf(out,"%d %d %.16f %.16f\n",i,position[i].t,x,y);
        }
fclose(out);
}

void setVerlet(){
int i,k;
double dx,dy;
double r2;

//To compare whether the monomers have left the Verlet radius, you need to know prior positions of monomers. Vx stores prior positions.
	for (i=1;i<=nT;i++){
	Nverl[i]=0;   //initialization
	Vx[i]=position[i].x;
	Vy[i]=position[i].y;
	}


	for (i=1;i<=nT;i++){
		for (k=i+1;k<=nT;k++){
		dx=position[i].x-position[k].x;
		dy=position[i].y-position[k].y;

		//PBC
		//dx -= side*nearbyint(dx*side_r);
		//dy -= side*nearbyint(dy*side_r);
			if (dx>side/2) dx-=side;
		        if (dx<-side/2) dx+=side;
			if (dy>side/2) dy-=side;
			if (dy<-side/2) dy+=side;


		r2=dx*dx+dy*dy;
			if (r2<rv2){
			Nverl[i]++;
			verl[i][Nverl[i]]=k;
			Nverl[k]++;
			verl[k][Nverl[k]]=i;
			}
		}
	}
}

void checkVerlet(int k){
//Checks Verlet list and updates it if necessary.
double dx,dy,r2;

dx=position[k].x-Vx[k];  
dy=position[k].y-Vy[k];
//PBC
//dx -= side*nearbyint(dx*side_r);
//dy -= side*nearbyint(dy*side_r);
	if (dx>side/2) dx-=side;
        if (dx<-side/2) dx+=side;
	if (dy>side/2) dy-=side;
	if (dy<-side/2) dy+=side;

r2=dx*dx+dy*dy;

	if(r2>Vgap2) setVerlet();
}

void count_neighbors(){
long double dx,dy,d0,separation;
int ti,tj,n;
//initialization
        for (int i=1; i<=nT; i++){
        no_neighbors[i]=0;
        }

//Find number of neighbors
        for (int i=1; i<=nT; i++){
        ti=position[i].t;
                for (int j=i+1; j<=nT; j++){
                tj=position[j].t;
                        if (ti+tj==2) d0=2.*rA;
                        else if (ti+tj==3) d0=rA+rB;
                        else if (ti+tj==4) d0=2.*rB;

                dx=fabs(position[i].x-position[j].x);
                dy=fabs(position[i].y-position[j].y);

		if (dx>side/2) dx-=side;
		if (dy>side/2) dy-=side;
                separation = dx*dx + dy*dy;
                separation = sqrt(separation);
    
                        if (separation<1.00*d0) {
			//if (d0-separation>1e-11) {
                        no_neighbors[i]++; 
                        no_neighbors[j]++;
                        }
                }
        }

}

void remove_rattlers(){

int n=0;
vec position_noRattlers[nT+1];
vec total_force_noRattlers[nT+1];
vec forceInteraction_noRattlers[nT+1];
vec velocity_noRattlers[nT+1];
double activeDirector_noRattlers[20000][2];
int k=0;

	while (n!=nT){
		//if (k>0) {nT=n; printf("%d\n",nT);}
		if (k>0) nT=n; 
	count_neighbors();
	k++;
	n=0;
		for (int i=1; i<=nT; i++){
			if (no_neighbors[i]>2){
			n++;
			position_noRattlers[n].t = position[i].t;
			position_noRattlers[n].x = position[i].x;
			position_noRattlers[n].y = position[i].y;

			activeDirector_noRattlers[n][0] = activeDirector[i][0];
			activeDirector_noRattlers[n][1] = activeDirector[i][1];

			velocity_noRattlers[i].x = velocity[i].x;
			velocity_noRattlers[i].y = velocity[i].y;

			total_force_noRattlers[i].x = total_force[i].x;
			total_force_noRattlers[i].y = total_force[i].y;

			forceInteraction_noRattlers[i].x = forceInteraction[i].x;
			forceInteraction_noRattlers[i].y = forceInteraction[i].y;
			}
		}

	//Remove rattlers from position list and create postions list using only no rattlers
		for (int i=1; i<=nT; i++){
		position[i].t=0;
		position[i].x=0.;
		position[i].y=0.;

		activeDirector[i][0]=0.;
		activeDirector[i][1]=0.;

		velocity[i].x=0.;
		velocity[i].y=0.;

		total_force[i].x=0.;
		total_force[i].y=0.;

		forceInteraction[i].x=0.;
		forceInteraction[i].y=0.;
		}


		for (int i=1; i<=n; i++){
		position[i].t=position_noRattlers[i].t;
		position[i].x=position_noRattlers[i].x;
		position[i].y=position_noRattlers[i].y;
		
		activeDirector[i][0]=activeDirector_noRattlers[i][0];
		activeDirector[i][1]=activeDirector_noRattlers[i][1];

		velocity[i].x=velocity_noRattlers[i].x ;
		velocity[i].y=velocity_noRattlers[i].y ;

		total_force[i].x=total_force_noRattlers[i].x;
		total_force[i].y=total_force_noRattlers[i].y;

		forceInteraction[i].x=forceInteraction_noRattlers[i].x;
		forceInteraction[i].y=forceInteraction_noRattlers[i].y;
		}

	}

setVerlet();
}

double calculate_deltaZ(){
//Calculate neighbors
		
		//initialization
		for (int i=1; i<=nT; i++){
		no_neighbors[i]=0;
		}

int ti,tj;
double dx,dy,d0,separation;

		for (int i=1; i<=nT; i++){
		ti=position[i].t;
	
			for (int j=1;j<=Nverl[i];j++){     
			tj=position[verl[i][j]].t;
				if (ti+tj==2) d0=2.*rA;
				else if (ti+tj==3) d0=rA+rB;
				else if (ti+tj==4) d0=2.*rB;

			dx=fabs(position[i].x-position[verl[i][j]].x);
			dy=fabs(position[i].y-position[verl[i][j]].y);

			if (dx>side/2) dx-=side;
			if (dy>side/2) dy-=side;
			separation = dx*dx + dy*dy;
			separation = sqrt(separation);
		
				if (separation<1.00*d0) {
				//if (d0-separation>1e-11) {
				no_neighbors[i]++; 
				}
			}
		}


	
float z_average_no_rattlers=0;
int n=0;
		for (int i=1; i<=nT; i++){
			if (no_neighbors[i]>2){
			n++;
			z_average_no_rattlers += no_neighbors[i];
			}
		}
	z_average_no_rattlers /= n;

return z_average_no_rattlers - 4.; 
}

void write_config(){
FILE *out;
out=fopen("config.in","w");
fclose(out);

int k=1;
out=fopen("config.in","a");
fprintf(out,"%d %.16f\n",nT,side);
	for (int i=1;i<=nT;i++){
 	fprintf(out,"%d %d %.16f %.16f\n",i,position[i].t,position[i].x,position[i].y);
	} 
fprintf(out,"\n");

	for (int i=1;i<=nT;i++){
	fprintf(out,"%d %.16f %.16f\n",i,activeDirector[i][0],activeDirector[i][1]);
	}
fclose(out);
}

void printGapDist(int i){
double r,dx,dy;
int ti,tj;
double d0;

ti=position[i].t;

	for (int j=1;j<=Nverl[i];j++){     
	
	tj=position[verl[i][j]].t;
                        if      (ti+tj==2) d0=2.*rA;
                        else if (ti+tj==3) d0=rA+rB;
                        else if (ti+tj==4) d0=2.*rB;
	dx=position[i].x-position[ verl[i][j] ].x;
        dy=position[i].y-position[ verl[i][j] ].y;

		if (dx>side/2.) dx-=side;
		if (dx<-side/2.) dx+=side;
		if (dy>side/2.) dy-=side;
		if (dy<-side/2.) dy+=side;

        r=(dx*dx+dy*dy);
	r=sqrt(r);

		if (i<verl[i][j]){
		printf("%.16f\n",r-d0);
		}
	}
}


void dynamic_printHessian(){
double r,dx,dy;
int ti,tj;
double d0;
int k;
int l,m;
double diagx,diagy,nonDiag;

//initialize
double* hessian;
hessian = (double*)malloc((int)(2*nT*2*nT)*sizeof(double));
	for (int i=0; i<(int)(2*nT*2*nT); i++){
	hessian[i]=0.;
	}

	for (int i=1; i<=nT; i++){
	ti=position[i].t;

		for (int j=1;j<=Nverl[i];j++){     
		k=verl[i][j];
		tj=position[k].t;
				if      (ti+tj==2) d0=2.*rA;
				else if (ti+tj==3) d0=rA+rB;
				else if (ti+tj==4) d0=2.*rB;
		dx=position[i].x-position[k].x;
		dy=position[i].y-position[k].y;

			if (dx>side/2) dx-=side;
			if (dx<-side/2) dx+=side;
			if (dy>side/2) dy-=side;
			if (dy<-side/2) dy+=side;

		r=(dx*dx+dy*dy);
		r=sqrt(r);

			if (r<d0 && i<k){
			l=(int)(2*(i-1));  //I am converting from particle index starting from 1 to starting from 0
			m=(int)(2*(k-1));  //l,m indices identify the derivative of potential energy with xcoord of particles i,k. ycoord is l+1

			//Diagonal 
			diagx=doubleDerivPotDiag(dx,r,d0);
			diagy=doubleDerivPotDiag(dy,r,d0);
			//hessian[l][l] += diagx;
			hessian[(int)(l*2*nT)+l] += diagx;
			
			//hessian[l+1][l+1] += diagy;
			hessian[(int)((l+1)*2*nT)+l+1] += diagy;

			//hessian[m][m] += diagx;
			hessian[(int)(m*2*nT)+m] += diagx;
			//hessian[m+1][m+1] += diagy;
			hessian[(int)((m+1)*2*nT)+m+1] += diagy;

			//Non-diagonal
			//xixk,yiyk
			//hessian[l][m] -= diagx;;
			hessian[(int)(l*2*nT)+m] += -diagx;;
			//hessian[m][l] -= diagx;;
			hessian[(int)(m*2*nT)+l] += -diagx;;

			//hessian[l+1][m+1] -= diagy;
			hessian[(int)((l+1)*2*nT)+m+1] += -diagy;
			//hessian[m+1][l+1] -= diagy;
			hessian[(int)((m+1)*2*nT)+l+1] += -diagy;
			
			//xiyi,xkyk
			nonDiag=doubleDerivPotNonDiag(dx,dy,r,d0);
			//hessian[l][l+1] += nonDiag;
			hessian[(int)(l*2*nT)+l+1] += nonDiag;
			//hessian[l+1][l] += nonDiag;
			hessian[(int)((l+1)*2*nT)+l] += nonDiag;

			//hessian[m][m+1] += nonDiag;
			hessian[(int)(m*2*nT)+m+1] += nonDiag;
			//hessian[m+1][m] += nonDiag;
			hessian[(int)((m+1)*2*nT)+m] += nonDiag;

			//xiyk,yixk
			//hessian[l][m+1] -= nonDiag;
			hessian[(int)(l*2*nT)+m+1] += -nonDiag;
			//hessian[m+1][l] -= nonDiag;
			hessian[(int)((m+1)*2*nT)+l] += -nonDiag;

			//hessian[l+1][m] -= nonDiag;
			hessian[(int)((l+1)*2*nT)+m] += -nonDiag;
			//hessian[m][l+1] -= nonDiag;
			hessian[(int)(m*2*nT)+l+1] += -nonDiag;
			}
		}
	}


	for (int i=0; i<int(2*nT*2*nT); i++){
		if (i%(int)(2*nT)==0 && i!=0) printf("\n");
		if (hessian[i]==0.0) printf("%.2f ",hessian[i]);
		else printf("%.16f ",hessian[i]);

	}

free(hessian);
}

void printHessian(){
double r,dx,dy;
int ti,tj;
double d0;
int k;
int l,m;
double diagx,diagy,nonDiag;

//initialize
int si=(int)(2*nT+1);
double hessian[si][si];

	for (int i=1; i<=(int)(2*nT); i++){
		for (k=1; k<=(int)(2*nT); k++){
		hessian[i][k]=0.;
		}
	}

	for (int i=1; i<=nT; i++){
	ti=position[i].t;

		for (int j=1;j<=Nverl[i];j++){     
		k=verl[i][j];
		tj=position[k].t;
				if      (ti+tj==2) d0=2.*rA;
				else if (ti+tj==3) d0=rA+rB;
				else if (ti+tj==4) d0=2.*rB;
		dx=position[i].x-position[k].x;
		dy=position[i].y-position[k].y;

			if (dx>side/2) dx-=side;
			if (dx<-side/2) dx+=side;
			if (dy>side/2) dy-=side;
			if (dy<-side/2) dy+=side;

		r=(dx*dx+dy*dy);
		r=sqrt(r);

			if (r<d0 && i<k){	 //For every contact between particles
			l=(int)(2*i-1);  //I am converting from particle index starting from 1 to starting from 0
			m=(int)(2*k-1);  //l,m indices identify the derivative of potential energy with xcoord of particles i,k. ycoord is l+1

			//Diagonal 
			diagx=doubleDerivPotDiag(dx,r,d0);
			diagy=doubleDerivPotDiag(dy,r,d0);
			hessian[l][l] += diagx;
			
			hessian[l+1][l+1] += diagy;

			hessian[m][m] += diagx;
			hessian[m+1][m+1] += diagy;

			//Non-diagonal
			//xixk,yiyk
			hessian[l][m] += -diagx;;
			hessian[m][l] += -diagx;;

			hessian[l+1][m+1] += -diagy;
			hessian[m+1][l+1] += -diagy;
			
			//xiyi,xkyk
			nonDiag=doubleDerivPotNonDiag(dx,dy,r,d0);
			hessian[l][l+1] += nonDiag;
			hessian[l+1][l] += nonDiag;

			hessian[m][m+1] += nonDiag;
			hessian[m+1][m] += nonDiag;

			//xiyk,yixk
			hessian[l][m+1] += -nonDiag;
			hessian[m+1][l] += -nonDiag;

			hessian[l+1][m] += -nonDiag;
			hessian[m][l+1] += -nonDiag;
			}

		}
	}


	for (int i=1; i<=int(2*nT); i++){
	for (int k=1; k<=int(2*nT); k++){
	printf("%.16f ",hessian[i][k]);
	}
	printf("\n");
	}

}

double doubleDerivPotDiag(double dx,double r,double d0){
double dderiv=0.0;

dderiv += dx*dx/pow(r*d0,2);
dderiv += dx*dx*(1.0-r/d0)/(pow(r,3)*d0);
dderiv += -(1.0-r/d0)/(r*d0);

return dderiv;
}

double doubleDerivPotNonDiag(double dx,double dy,double r,double d0){
double dderiv=0.0;

dderiv += dx*dy/(pow(r*d0,2));
dderiv += dx*dy*(1.0-r/d0)/(pow(r,3)*d0);

return dderiv;
}
