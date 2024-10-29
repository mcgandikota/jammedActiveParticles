#include "./header.h"

int main(int argc, char *argv[]){
	if (argc!=3){
	printf("./a.out | config.in | activeSpeed |\n");
	exit(0);
	}

config=argv[1];
activeSpeed=atof(argv[2]);
frame=1000;

initialize_no_angles(); 					           //store positions, active velocity directors
write_config();

double deltaZ;
int i=0;
energy();
averageTotalForce();
deltaZ=calculate_deltaZ();
totF=1.2*totF_cutoff;
			   
        while (totF > totF_cutoff && i<1e9){
	i++;     					  //for while loop only
		if (i%frame==0){
		//energy();   				  //This is time consuming
		averageTotalForce();
		deltaZ=calculate_deltaZ();		  //This is time consuming
		printf("%d %.30f %.30f %.30f %.30f %.30f %.30f\n",i,deltaZ,sqrt(2.*ke),pe,pe_Eff,totF,virial_pressure);
		//print(i);		   //Print movie
		}
	//brownianRun();
	//mdrun();
	FIRE();
	}

printf("%.16f\n",deltaZ);

print(i);
return 0;
}
