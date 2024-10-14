//Energy minimization of a jammed active particle distribution in 2D

int max_iterations=10000;
long double density_c[10000],density[10000];
long double U[10000];
#include "./longDoubleHeader.h"

void calculate_density(int k);
void dilate(int k);
void minimizePE(int k);
void calculate_critical_density(int k);

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
//density[0]=1.344;                                  //Phic=0.844, density[0] shoudl be twice this acording to PRL.
initialize_random();

        for (int i=0; i<=nT; i++){
	totalForce(i);
	}
							   
//Calculate U0


int k=0;
pe=1e-10;
minimizePE(k);
density_c[0]=0.844;
printf("%d %.10Le %.10Le %.10Le %.10Le %.10Le \n",k,pe,density[k],density_c[k],sqrt(U[k]/U[k-1]),virial_pressure);

k=1;
calculate_density(k);
//printf("%lf\n",density[1]);
dilate(k);
//printf("%lf %lf\n",rA,rB);
minimizePE(k);
density_c[1]=0.844;				   //Cannot use iteration formula since you need U[k-2]
printf("%d %.10Le %.10Le %.10Le %.10Le %.10Le \n",k,pe,density[k],density_c[k],sqrt(U[k]/U[k-1]),virial_pressure);

	//for (k=2;k<=max_iterations;k++){
	while (pe>1e-24){
	k++;
	calculate_density(k);
	dilate(k);
	minimizePE(k);
	calculate_critical_density(k);
	printf("%d %.10Le %.10Le %.10Le %.10Le %.10Le \n",k,pe,density[k],density_c[k],sqrt(U[k]/U[k-1]),virial_pressure);
	}

return 0;
}

void minimizePE(int k){
long double force_cutoff=4e-14;
force_cutoff=4e-12;
totF=1.2*force_cutoff;				   //Just for while loop to run for the first step
int i=0;
        while (totF>force_cutoff){
		if (i%frame==0){
		energy();
		averageTotalForce();
		//printf("%d %.30f %.30f %.30f %.30f %.30f\n",i,pe,ke,ke_COM,totF,virial_pressure);
		print(i);		   //Print movie
		}
	FIRE();
	i++;     					  //for while loop only
	}
U[k]=pe;
}
void calculate_density(int k){
density[k] = density_c[k-1] + (density[k-1]-density_c[k-1])*pow(10,-1./10);
}

void dilate(int k){
long double packFrac=density[k];
long double alpha = sqrt((2*packFrac*side*side)/(nT*M_PI*(rA*rA+rB*rB)));
rA = alpha*rA;
rB = alpha*rB;
}

void calculate_critical_density(int k){
density_c[k]=(density[k]-density[k-1]*sqrt(U[k-1]/U[k-2]))/(1.-sqrt(U[k-1]/U[k-2]));
}
