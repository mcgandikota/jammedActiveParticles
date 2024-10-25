//Energy minimization of a jammed active particle distribution in 2D

long double density_c[10000],density[10000];
long double U[10000];
long double delta;
int nT_noRattlers;
int contacts_noRattlers;
#include "./longDoubleHeader.h"
//#include "./header.h"

void calculate_density(int k);
void dilate(int k);
void minimizePE(int k);
void calculate_critical_density(int k);
void remove_rattlers();

int main(int argc, char *argv[]){
	if (argc!=4){
	printf("./a.out | N |steps/frame | #frames \n");
	exit(0);
	}

nT=atoi(argv[1]);
frame=atoi(argv[2]);
steps=(int)(frame*atoi(argv[3]));
activeSpeed=0.0;

density[0]=0.944;                                  //Phic=0.844, density[0] shoudl be twice this acording to PRL.
initialize_random(density[0]);
print(0);

        for (int i=0; i<=nT; i++){
	totalForce(i);
	}
							   
//Calculate U0
printf("k | pe | density | density_c | sqrt(U[k]/U[k-1]) | pressure | delta\n");
int k=0;
pe=1e-10;
minimizePE(k);
density_c[0]=0.904;
printf("%d %.10Le %.10Le %.10Le %.10Le %.10Le %.f\n",k,pe,density[k],density_c[k],sqrt(U[k]/U[k-1]),virial_pressure,delta);

damp=0.01;
deltaT=0.02;					    //You can use larger deltaT here since the first minimization which involved large forces is done.
						    
	while (pe>1e-20){
	k++;
	calculate_density(k);
	dilate(k);

	minimizePE(k);
	calculate_critical_density(k);
	printf("%d %.10Le %.10Le %.10Le %.10Le %.10Le %.10Le\n",k,pe,density[k],density_c[k],sqrt(U[k]/U[k-1]),virial_pressure,delta);
	}

//Rescale
rB /= rA;
side /= rA;
rA = 1.0;

remove_rattlers();

printf("%d %d %d\n",contacts_noRattlers,nT_noRattlers,(nT_noRattlers-1)*2+1 );

return 0;
}

void minimizePE(int k){
long double force_cutoff=4e-14;
force_cutoff=1e-10;
totF=1.2*force_cutoff;				   //Just for while loop to run for the first step
int i=0;
        while (totF>force_cutoff){
	i++;     					  //for while loop only
	FIRE();
		if (i%frame==0){
		energy();
		averageTotalForce();
		//printf("%d %.30f %.30f %.30f %.30f %.30f\n",i,pe,ke,ke_COM,totF,virial_pressure);
		print(i);		   //Print movie
		}
	}
U[k]=pe;

}
void calculate_density(int k){
density[k] = density_c[k-1] + (density[k-1]-density_c[k-1])*pow(10,-1./10);
}

void dilate(int k){
long double packFrac=density[k];
delta = sqrt((2*packFrac*side*side)/(nT*M_PI*(rA*rA+rB*rB)));
rA = delta*rA;
rB = delta*rB;
}

void calculate_critical_density(int k){
//density_c[k]=(density[k]-density[k-1]*sqrt(U[k-1]/U[k-2]))/(1.-sqrt(U[k-1]/U[k-2]));  //This expression in PRL has a typo. Shastry's paper has correct expression
density_c[k]=(density[k]-density[k-1]*sqrt(U[k]/U[k-1]))/(1.-sqrt(U[k]/U[k-1]));
}

void remove_rattlers(){
int no_neighbors[nT+1];
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

                if (dx>side/2.) dx -= side;
                if (dy>side/2.) dy -= side;
                separation = dx*dx + dy*dy;
                separation = sqrt(separation);
    
                        if (separation<1.00*d0) {
                        no_neighbors[i]++; 
                        no_neighbors[j]++;
                        }
                }
        }

float z_average_no_rattlers=0;
double diameter;
vec position_noRattlers[nT];
n=0;
        for (int i=1; i<=nT; i++){
                if (no_neighbors[i]>2){
                n++;
		position_noRattlers[n].t = position[i].t;
		position_noRattlers[n].x = position[i].x;
		position_noRattlers[n].y = position[i].y;
                z_average_no_rattlers += no_neighbors[i];
                }
        }
nT_noRattlers=n;
contacts_noRattlers=0;
        for (int i=1; i<=nT_noRattlers; i++){
        ti=position_noRattlers[i].t;
                for (int j=i+1; j<=nT_noRattlers; j++){
                tj=position_noRattlers[j].t;
                        if (ti+tj==2) d0=2.*rA;
                        else if (ti+tj==3) d0=rA+rB;
                        else if (ti+tj==4) d0=2.*rB;

                dx=fabs(position_noRattlers[i].x-position_noRattlers[j].x);
                dy=fabs(position_noRattlers[i].y-position_noRattlers[j].y);

                if (dx>side/2.) dx -= side;
                if (dy>side/2.) dy -= side;
                separation = dx*dx + dy*dy;
                separation = sqrt(separation);
    
                        if (abs(separation-d0)<(1e-11)*d0) {
			contacts_noRattlers ++;
                        }
                }
        }
}
