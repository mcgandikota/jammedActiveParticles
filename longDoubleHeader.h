//Energy minimization of a jammed active particle distribution in 2D
//
using namespace std;
#include <math.h>
#include <tgmath.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include "time.h"

typedef struct{
long double x,y;
int t;
}vec;

//Parameters, variables and arrays
int steps;
int frame;
int nT;
long double side;
long double side_r;      //inverse of side
long double m=1;
long double damp=0.01;					           //Is this too much damping for inertial FIRE to work efficiently?
long double deltaT=0.2;					   //deltaT=0.2 works good for MD
long double totF_cutoff=14e-20;
//long double deltaT=0.0001;					   //deltaT=0.2 works good for MD
long double Gamma=exp(-damp*deltaT/m);
long double c1=m/damp*(1-Gamma);
long double c2=m/(damp*damp)*(damp*deltaT/m-1+Gamma);
long double rA=1.3999999999999999/2.;
long double rB=0.5;						
long double rv,rv2,rc,rc2,Vgap,Vgap2;
vec position[20000];
vec velocity[20000];
vec forceInteraction[20000];
vec total_force[20000];
long double totF;
long double activeDirector[20000][2];
long double activeSpeed;
long double ke;
long double ke_COM;                      
long double pe;
long double pe_Eff;
long double Vx[20000];    // Positions of monomers in past
long double Vy[20000];    // Positions of monomers in past
int Nverl[20000];    // The number of monomers in neighborhood of i^th monomer
int verl[20000][20]; // Second index identities monomers in neighborhood of i^th monomer in first index. There can't be more than about 12 monomers in neighborhood. 20 is a sufficiently big enough list. The first Nverl[i] entries in second index hold the identities of monomers of i^th monomer in first index. All succeeding entries in second index will remain zeroes.
long double virial_pressure=0.0;
char * config;

//FIRE parameters
long double alpha;
long double alphaStart=0.1;
long double fInc=1.1;
long double deltaTmax=10*deltaT;
long double fAlpha=0.99;
int Nmin=5;
int powerPositiveSteps;
long double fDec=0.5;

//Function declarations
void initialize();
void initialize_random();
void initialize_LAMMPS_config();
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


//Functions
void initialize(){
FILE *inp;
char lines[1000];
char st;
long double v;
//get nT and side
inp=fopen(config,"r");  
fgets(lines,1000,inp);
sscanf(lines, "%d %lf ",&nT,&side);
fgets(lines,1000,inp);

//Initialize Positions
int q,t;
long double x,y,z;
        for (int i=1; i<=nT; i++){
        fgets(lines,1000,inp);
	sscanf(lines,"%d %d %lf %lf",&q,&t,&x,&y);
	position[i].t=t;
	//position[i].x=x+side/2.;  //The box is now [0,side]
	//position[i].y=y+side/2.;
	position[i].x=x;  
	position[i].y=y;
	}
        fgets(lines,1000,inp);
        for (int i=1; i<=nT; i++){
        fgets(lines,1000,inp);
	sscanf(lines, "%d %lf %lf",&q,&x,&y);
     	activeDirector[i][0]=x;
	activeDirector[i][1]=y;
        }
fclose(inp);

//Initialize Velocities and Forces
        for (int i=1; i<=nT; i++){
	velocity[i].x=0.;
	velocity[i].y=0.;

	forceInteraction[i].x=0.;
	forceInteraction[i].y=0.;
	total_force[i].x=0.;
	total_force[i].x=0.;
	}

//Initialize Verlet List
//!!!!!!!!!!!!!!Can rv be smaller 1.2* ...?
rv=1.8*2.*rA; 			       			   //Outer diameter in Verlet list
	      						   //Choosing the bigger particle's radius
rv2=pow(rv,2);
rc= 2.*rA;						   //Inner diameter in Verlet list
rc2=rc*rc;
Vgap=0.5*(rv-rc);  					   //Radius of annulus in Verlet method
Vgap2=pow(Vgap,2);

setVerlet();


FILE *out;
out=fopen("out.dump","w");
fclose(out);

//Initialize FIRE
powerPositiveSteps=0;
alpha=alphaStart;
}

void initialize_random(){
side = sqrtl(nT*M_PI*(rA*rA+rB*rB)/(2*density[0]));
double dx,dy,dr;

//Initialize Positions
srand48(time(NULL)); 					   //takes current time as seed for random number generator

	for (int i=1;i<=nT;i++){
		if (i<=(int)(nT/2)) position[i].t=1;
		else     	    position[i].t=2;
	here:
	position[i].x=(0.5-drand48())*0.99*side; 			   //0.99 so that it wont end up on the edge of PBC
	position[i].y=(0.5-drand48())*0.99*side;

		for (int j=1;j<i;j++){
		dx=position[i].x-position[j].x;
		dy=position[i].y-position[j].y; 
		dr=sqrtl(dx*dx+dy*dy);
			if (dr<0.4) goto here;
		} 
	activeDirector[i][0]=0.;			  //We will not use these directors for these simulations
	activeDirector[i][1]=0.;

	//printf("%d ",i);
	}
printf("\n");

//Initialize Velocities and Forces
        for (int i=1; i<=nT; i++){
	velocity[i].x=0.;
	velocity[i].y=0.;

	forceInteraction[i].x=0.;
	forceInteraction[i].y=0.;
	total_force[i].x=0.;
	total_force[i].x=0.;
	}

//Initialize Verlet List
//!!!!!!!!!!!!!!Can rv be smaller 1.2* ...?
rv=1.8*2.*rA; 			       			   //Outer diameter in Verlet list
	      						   //Choosing the bigger particle's radius
rv2=pow(rv,2);
rc= 2.*rA;						   //Inner diameter in Verlet list
rc2=rc*rc;
Vgap=0.5*(rv-rc);  					   //Radius of annulus in Verlet method
Vgap2=pow(Vgap,2);

setVerlet();


FILE *out;
out=fopen("out.dump","w");
fclose(out);

//Initialize FIRE
powerPositiveSteps=0;
alpha=alphaStart;
}

void initialize_LAMMPS_config(){
FILE *inp;
char lines[1000];
char st;
long double v;
inp=fopen(config,"r");  //get nT and side
        for (int i=1; i<=3; i++){
        fgets(lines,1000,inp);
        }
sscanf(lines, "%d %s",&nT,&st);
        for (int i=1; i<=11; i++){
        fgets(lines,1000,inp);
        }
sscanf(lines,"%lf %lf %c %c",&side,&v,&st,&st);
side=-2.0*side;
side_r=1./side;
fclose(inp);

//Initialize Positions
int q,t;
long double x,y,z;
inp=fopen(config,"r");  
        for (int i=1; i<=23; i++){
        fgets(lines,1000,inp);
        }
        for (int i=1; i<=nT; i++){
        fgets(lines,1000,inp);
	sscanf(lines,"%d %d %lf %lf %lf %lf %lf %lf",&q,&t,&x,&y,&z,&v,&v,&v);
	position[i].t=t;
	//position[i].x=x+side/2.;  //The box is now [0,side]
	//position[i].y=y+side/2.;
	position[i].x=x;  
	position[i].y=y;
	}
        for (int i=1; i<=3; i++){
        fgets(lines,1000,inp);
	}
        for (int i=1; i<=nT; i++){
        fgets(lines,1000,inp);
	sscanf(lines, "%lf %lf %lf %lf %lf %lf %lf %lf",&z,&z,&z,&z,&x,&z,&z,&y);
     	activeDirector[i][0]=x*x-y*y;
	activeDirector[i][1]=2*x*y;
        }
fclose(inp);

//Initialize Velocities and Forces
        for (int i=1; i<=nT; i++){
	velocity[i].x=0.;
	velocity[i].y=0.;

	forceInteraction[i].x=0.;
	forceInteraction[i].y=0.;
	total_force[i].x=0.;
	total_force[i].x=0.;
	}

//Initialize Verlet List
//!!!!!!!!!!!!!!Can rv be smaller 1.2* ...?
rv=1.8*2.*rA; 			       			   //Outer diameter in Verlet list
	      						   //Choosing the bigger particle's radius
rv2=pow(rv,2);
rc= 2.*rA;						   //Inner diameter in Verlet list
rc2=rc*rc;
Vgap=0.5*(rv-rc);  					   //Radius of annulus in Verlet method
Vgap2=pow(Vgap,2);

setVerlet();


FILE *out;
out=fopen("out.dump","w");
fclose(out);

//Initialize FIRE
powerPositiveSteps=0;
alpha=alphaStart;
}


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
long double power=0.;
	for(int i=1; i<=nT; i++){
	power += total_force[i].x*velocity[i].x;
	power += total_force[i].y*velocity[i].y;
	}

	if (power>0.) powerPositiveSteps += 1;
	else 	      powerPositiveSteps  = 0;


//Calculate unit vector of total force
long double force_mag=0.;
	for(int i=1; i<=nT; i++){
	force_mag += total_force[i].x*total_force[i].x;
	force_mag += total_force[i].y*total_force[i].y;
	}
force_mag=sqrtl(force_mag);

vec unit_total_force[20000];
	for(int i=1; i<=nT; i++){
	unit_total_force[i].x = total_force[i].x/force_mag;
	unit_total_force[i].y = total_force[i].y/force_mag;
	}

//Calculate magnitude of velocity
long double modV=0.;
	for(int i=1; i<=nT; i++){
	modV += velocity[i].x*velocity[i].x;
	modV += velocity[i].y*velocity[i].y;
	}
modV=sqrtl(modV);

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
long double fxi,fyi;
long double r,dx,dy;
int ti,tj;
long double d0,ks=1.0;

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
	r=sqrtl(r);
		if (r<d0){
		fxi= -ks*(r-d0)*(dx)/r;             //Check for 1/2
                fyi= -ks*(r-d0)*(dy)/r;

		forceInteraction[i].x += fxi;
		forceInteraction[i].y += fyi;

		virial_pressure += dx*fxi + dy*fyi;
		}
	}

	if (i==nT) virial_pressure /= 3.*2*side*side;  //Accounting for overcounting too
}

void totalForce(int i){
forceInt(i);  					   //Calculate interaction force
total_force[i].x=forceInteraction[i].x+activeSpeed*activeDirector[i][0];
total_force[i].y=forceInteraction[i].y+activeSpeed*activeDirector[i][1];
}

void averageTotalForce(){
totF=0.;
long double f;
        for (int i=1;i<=nT;i++){
	f=total_force[i].x*total_force[i].x+total_force[i].y*total_force[i].y;
	f=sqrtl(f);
	totF += f;
	}
totF /= nT;
}

void energy(){
ke=0.;
ke_COM=0.;                      //KE in the COM frame
long double vx_COM=0.,vy_COM=0.;
long double dx,dy,r;
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
long double d0,ks=1.0;
        for (int i=1;i<=nT;i++){
	ti=position[i].t;
	pe_Eff -= activeDirector[i][0]*position[i].x+activeDirector[i][1]*position[i].y;
	//Use PBC position or unwrapped position? 

		for (int j=i+1;j<=nT;j++){
		tj=position[j].t;
                        if      (ti+tj==2) d0=2.*rA; 
                        else if (ti+tj==3) d0=rA+rB; 
                        else if (ti+tj==4) d0=2.*rB; 
		dx=position[i].x-position[j].x;
		dy=position[i].y-position[j].y;

		//PBC
		//dx -= side*nearbyint(dx*side_r);
		//dy -= side*nearbyint(dy*side_r);
			if (dx>side/2) dx-=side;
		        if (dx<-side/2) dx+=side;
			if (dy>side/2) dy-=side;
			if (dy<-side/2) dy+=side;


		r=(dx*dx+dy*dy);
			if (r<d0){
			pe += ks*(r-d0)*(r-d0);                           
			}
		}
	}
pe *= 0.5;
pe_Eff *= activeSpeed;
pe_Eff += pe;
pe_Eff /=(2*nT);
pe/=(2*nT);
}

void print(int k){
long double x,y;
FILE *out;
out=fopen("out.dump","a");
fprintf(out,"ITEM: TIMESTEP\n");
fprintf(out,"%d\n",k);
fprintf(out,"ITEM: NUMBER OF ATOMS\n");
fprintf(out,"%d\n",nT);
fprintf(out,"ITEM: BOX BOUNDS pp pp pp\n");
fprintf(out,"-%.16Lf %.16Lf \n",side/2.,side/2.);
fprintf(out,"-%.16Lf %.16Lf \n",side/2.,side/2.);
fprintf(out,"-%.16Lf %.16Lf \n",0.7,0.7);
fprintf(out,"ITEM: ATOMS id type x y\n");

        for (int i=1;i<=nT;i++){
	//Printing PBC
	x=position[i].x-floor(position[i].x/side)*side -side/2.;
	y=position[i].y-floor(position[i].y/side)*side -side/2.;
	fprintf(out,"%d %d %.10Le %.10Le \n",i,position[i].t,x,y);
        }
fclose(out);
}

void setVerlet(){
int i,k;
long double dx,dy;
long double r2;

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
long double dx,dy,r2;

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
