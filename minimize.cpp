//Energy minimization of a jammed active particle distribution in 2D
//
using namespace std;
#include <math.h>
#include <cstdio>
#include <iostream>

typedef struct{
double x,y;
int t;
}vec;

void initialize();
void mdrun();
void forceInt(int i);
void print();
void setVerlet(void);
void checkVerlet(int k);

int steps=100000;
int frame=1000;
int nT;
double side;
double m=1;
double damp=1;					           //Is this too much damping for inertial FIRE to work efficiently?
double deltaT=0.2;
double Gamma=exp(-damp*deltaT/m);
double c1=m/damp*(1-Gamma);
double c2=m/(damp*damp)*(damp*deltaT/m-1+Gamma);
double rA=1.3999999999999999/2.;
double rB=0.5;
double rv,rv2,rc,rc2,Vgap,Vgap2;
vec position[20000];
vec velocity[20000];
vec forceInteraction[20000];
double activeDirector[20000][2];
double activeSpeed;
double Vx[20000];    // Positions of monomers in past
double Vy[20000];    // Positions of monomers in past
int Nverl[20000];    // The number of monomers in neighborhood of i^th monomer
int verl[20000][20]; // Second index identities monomers in neighborhood of i^th monomer in first index. There can't be more than about 12 monomers in neighborhood. 20 is a sufficiently big enough list. The first Nverl[i] entries in second index hold the identities of monomers of i^th monomer in first index. All succeeding entries in second index will remain zeroes.
char * config;


int main(int argc, char *argv[]){
	if (argc!=3) {printf("./a.out config.in activeSpeed\n");exit(0);}
config=argv[1];
activeSpeed=atof(argv[2]);

initialize(); 					           //store positions, active velocity directors
print();
        for (int i=1; i<=steps; i++){
	mdrun();
		if (i%frame==0){
		printf("%d\n",i);
		print();		   //Print movie
		}
	}


//while(ftol<FTOL)  {mdrun();}  //ftol is the change in contact forces between consecutive (10000) mdruns.
return 0;
}

void initialize(){
FILE *inp;
char lines[1000];
char st;
double v;
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
fclose(inp);

int q,t;
double x,y,z;
inp=fopen(config,"r");  //get nT and side
        for (int i=1; i<=23; i++){
        fgets(lines,1000,inp);
        }
        for (int i=1; i<=nT; i++){
        fgets(lines,1000,inp);
	sscanf(lines,"%d %d %lf %lf %lf %lf %lf %lf",&q,&t,&x,&y,&z,&v,&v,&v);
	position[i].t=t;
	position[i].x=x+side/2.;  //The box is now [0,side]
	position[i].y=y+side/2.;
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
}

void mdrun(){

	for(int i=1; i<=nT; i++){
	forceInt(i);
	//Position Update
	position[i].x=position[i].x+ c1*velocity[i].x+c2*(forceInteraction[i].x+activeSpeed*activeDirector[i][0]);
	position[i].y=position[i].y+ c1*velocity[i].y+c2*(forceInteraction[i].y+activeSpeed*activeDirector[i][1]);
	//Velocity Update
	velocity[i].x=Gamma*velocity[i].x+1./damp*(1.-Gamma)*(forceInteraction[i].x+activeSpeed*activeDirector[i][0]);
	velocity[i].y=Gamma*velocity[i].y+1./damp*(1.-Gamma)*(forceInteraction[i].y+activeSpeed*activeDirector[i][1]);

	checkVerlet(i);
	}
}

void forceInt(int i){
double fxi,fyi;
double r,dx,dy;
int ti,tj;
double d0,ks;

forceInteraction[i].x=0.0;
forceInteraction[i].y=0.0;

ti=position[i].t;

	for (int j=1;j<=Nverl[i];j++){
	tj=position[j].t;
                        if      (ti+tj==2) {d0=2.*rA; ks=0.255;}
                        else if (ti+tj==3) {d0=rA+rB; ks=0.347;}
                        else if (ti+tj==4) {d0=2.*rB; ks=0.5;}

        dx=fmod(position[i].x,side)-fmod(position[ verl[i][j] ].x,side);
        dy=fmod(position[i].y,side)-fmod(position[ verl[i][j] ].y,side);
        
        	if 	(dx < -side/2.) dx += side;
                else if (dx > side/2.)  dx -= side;

                if      (dy < -side/2.) dy += side;
                else if (dy > side/2.) dy -= side;


        r=(dx*dx+dy*dy);
		if (r<d0){
		fxi= -2.*ks*(r-d0)*(dx)/r;
                fyi= -2.*ks*(r-d0)*(dy)/r;

		forceInteraction[i].x += fxi;
		forceInteraction[i].y += fyi;
		}
    	}
}

void print(){
double x,y;
FILE *out;
out=fopen("out.dump","a");
fprintf(out,"%d\n\n",nT);
        for (int i=1;i<=nT;i++){
	//Printing PBC
	x=position[i].x-floor(position[i].x/side)*side;
	y=position[i].y-floor(position[i].y/side)*side;
	fprintf(out,"%d %d %lf %lf 0.000000 %lf %lf\n",i,position[i].t,x,y,position[i].x,position[i].y);
        }
fclose(out);
}

void setVerlet(void){
int i,k;
double dx,dy;
double r2;

//To compare whether the monomers have left the Verlet radius, you need to know prior positions of monomers. Vx stores prior positions.
	for (i=1;i<=nT;i++){
	Nverl[i]=0;   //initialization
	Vx[i]=fmod(position[i].x,side);
	Vy[i]=fmod(position[i].y,side);
	}


	for (i=1;i<=nT;i++){
		for (k=i+1;k<=nT;k++){
		dx=fmod(position[i].x,side)-fmod(position[k].x,side);
		dy=fmod(position[i].y,side)-fmod(position[k].y,side);

			if 	(dx < -side/2.) dx += side;
			else if (dx > side/2.)  dx -= side;

			if      (dy < -side/2.) dy += side;
			else if (dy > side/2.) dy -= side;

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

dx=fmod(position[k].x,side)-Vx[k];  //Vx,Vy is already in fmod form
dy=fmod(position[k].y,side)-Vy[k];
	if 	(dx < -side/2.) dx += side;
	else if (dx > side/2.)  dx -= side;

	if      (dy < -side/2.) dy += side;
	else if (dy > side/2.) dy -= side;
r2=dx*dx+dy*dy;

	if(r2>Vgap2) setVerlet();
}
