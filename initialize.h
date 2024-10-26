//Function declarations
void initialize();
void initialize_random(double density);
void initialize_LAMMPS_config();

//Functions
void initialize(){
FILE *inp;
char lines[1000];
char st;
double v;
//get nT and side
inp=fopen(config,"r");  
fgets(lines,1000,inp);
sscanf(lines, "%d %lf ",&nT,&side);
fgets(lines,1000,inp);

//Initialize Positions
int q,t;
double x,y,z;
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

void initialize_random(double density){
side = sqrt(nT*M_PI*(rA*rA+rB*rB)/(2*density));
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
		dr=sqrt(dx*dx+dy*dy);
			if (dr<0.8) goto here;
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
side_r=1./side;
fclose(inp);

//Initialize Positions
int q,t;
double x,y,z;
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
