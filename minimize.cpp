//Energy minimization of a jammed active particle distribution in 2D

#include "./header.h"

int main(int argc, char *argv[]){
	if (argc!=5){
	printf("./a.out | config.in | activeSpeed | steps/frame | #frames \n");
	exit(0);
	}

config=argv[1];
activeSpeed=atof(argv[2]);
frame=atoi(argv[3]);
steps=(int)(frame*atoi(argv[4]));

initialize_LAMMPS_config(); 					           //store positions, active velocity directors
//initialize();

print(0);

int i=0;
//totF=1.2*totF_cutoff;					   //This is just for the while loop to work on the first run
        for (i=1; i<=nT; i++){
	forceInt(i);
	}
							   
        for (i=0; i<=steps; i++){
        //while (totF>totF_cutoff){
		if (i%frame==0){
		energy();
		averageTotalForce();
		printf("%d %.30f %.30f %.30f %.30f %.30f %.30f\n",i,pe,pe_Eff,ke,ke_COM,totF,virial_pressure);
		print(i);		   //Print movie
		}
	brownianRun();
	//mdrun();
	//FIRE();
	//i++;     					  //for while loop only
	}


return 0;
}
