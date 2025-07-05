//   Shallow Water Equations 2D Simulator (SWE2D)
//  ===========================================================================================
//   Square grid
//  ===========================================================================================
//   Version 1.0 - Mar 2025
//  ===========================================================================================
//   Computational Hydraulics Group - University of Zaragoza
//  ===========================================================================================

#include "define.h"
#include "lib/shallow_water.h"


// Inputs and storage
char filename[1024];
FILE *fp;
FILE *logFile;
FILE *qFile;

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


// MAIN CODE FUNTION
int main(){

    ///CASE NAME
    char case_name[] = "Abou_Cisse_2D_explosion_problem_sinB";
    char case_name_slashed[1024] = "";
    if (case_name != NULL){
        sprintf(case_name_slashed, "\\");
        strcat(case_name_slashed, case_name);
    }


// OPENING DATA STORAGE

	char tempfile[1024];

	//i=system ("rm outputFiles/*.out");	// Remove old files (linux only)


    directory_founder(tempfile, "outputFiles", case_name_slashed, "\\log", ".out");
	logFile=fopen(tempfile,"w");

	directory_founder(tempfile, "outputFiles", case_name_slashed, "\\discharge", ".out");
	qFile=fopen(tempfile,"w");

	#if DISPLAY
	Gnuplot gp;			// Graphical output pipe (linux only)
	#endif



// SIMULATOR HEADERS
	printf("\n   SHALLOW WATER EQUATIONS  <<<>>> 2D SIMULATOR");
	printf("\n   Version 1.0 - Mar 2025");
	printf("\n   Computational Hydraulics Group - University of Zaragoza");
	printf("\n   --------------------------------------------------");

	printf("\n\n>> Case:	%s", case_name);
	printf("\n");

	fprintf(logFile,"\n   SHALLOW WATER EQUATIONS <<<>>>2D SIMULATOR");
	fprintf(logFile,"\n   Version 1.0 - Mar 2025");
	fprintf(logFile,"\n   Computational Hydraulics Group - University of Zaragoza");
	fprintf(logFile,"\n   --------------------------------------------------");

	fprintf(logFile,"\n\n>> Case:   %s", case_name);
	fprintf(logFile,"\n");



// GEOMETRY DATA
	// Load geometry data
	int nX, nY; // Number of cells in X and Y directions
	int nCells;
	double dx; // Grid spacing
	double Xmin, Ymin; // Domain limits
	double Length, Width; // Domain size

	directory_founder(tempfile, "input", case_name_slashed, "\\elevation", ".input");
	read_raster_header(tempfile, &nX, &nY, &Xmin, &Ymin, &dx);

	nCells = nX*nY;
	Length = dx*nX;
	Width = dx*nY;

	printf("\n\n>> Geometry: ");
	printf("\n     nCells %d ", nCells);
	printf("\n     DeltaX %lf - (Xmin,Ymin) = (%lf,%lf) - Width %lf - Length %lf", dx, Xmin, Ymin, Width, Length);
	printf("\n");

	fprintf(logFile,"\n\n>> Geometry: ");
	fprintf(logFile,"\n     nCells %d ", nCells);
	fprintf(logFile,"\n     DeltaX %lf - (Xmin,Ymin) = (%lf,%lf) - Width %lf - Length %lf", dx, Xmin, Ymin, Width, Length);
	fprintf(logFile,"\n");


	// Change array declarations to 2D
	double *x, *y; // Cell center coordinates
	double area = dx*dx;

	x = (double*) malloc( (nX+1)*(nY+1)*sizeof(double) );
	y = (double*) malloc( (nX+1)*(nY+1)*sizeof(double) );

	// Initialize grid coordinates
	for(int j = 0; j <= nY; j++) { //rows
		for(int i = 0; i <= nX; i++) { //cols
			x[IDX(i,j,nX+1)] = Xmin + i*dx; // i for x-direction
			y[IDX(i,j,nX+1)] = Ymin + (nY-j)*dx; // j for y-direction
		}
	}

// MAIN VARIABLE DEFINITION

	// Declare variables
	clock_t CPUtime;

	int iter;										// Iteraction index
	int nout;										// Output index

	double *h, *qx, *qy, *qBx, *qBy;				// Conserved variables
	double *ux, *uy, *Bx, *By;						// Primitive variables
	double *n;										// 2D averaged Manning Coefficient
	double *zb;										// Bed surface level

	double *DU1;									// Variation conserved variable h in Dt
	double *DU2;									// Variation conserved variable qx in Dt
	double *DU3;									// Variation Conserved variable qy in Dt
	double *DU4;									// Variation conserved variable Bx in Dt
	double *DU5;									// Variation Conserved variable By in Dt

	double t, dt;
	double dtSW;

	// Mass error monitor
	double massWt0;									// Initial mass t = t0
	double massWtn;									// Mass at the time t = tn
	double divergence_loss;                         // Divergence count: div(hB); if not 0, then monopoles, and monopoles kinda bad
	double Qin;										// Discharge at intlet for the time t = t0
	double Qout;									// Discharge at outlet for the time t = t0
	double massWerror;
	double Qbalance;

	// Allocate memory

	h = (double*) malloc( nCells* sizeof(double) );
	qx = (double*) malloc( nCells* sizeof(double) );
	qy = (double*) malloc( nCells* sizeof(double) );
	qBx = (double*) malloc( nCells* sizeof(double) );
	qBy = (double*) malloc( nCells* sizeof(double) );

	ux = (double*) malloc( nCells* sizeof(double) );
	uy = (double*) malloc( nCells* sizeof(double) );
	Bx = (double*) malloc( nCells* sizeof(double) );
	By = (double*) malloc( nCells* sizeof(double) );

	n = (double*) malloc( nCells* sizeof(double) );
	zb = (double*) malloc( nCells* sizeof(double) );

	DU1 = (double*) malloc( nCells* sizeof(double) );
	DU2 = (double*) malloc( nCells* sizeof(double) );
	DU3 = (double*) malloc( nCells* sizeof(double) );
	DU4 = (double*) malloc( nCells* sizeof(double) );
	DU5 = (double*) malloc( nCells* sizeof(double) );

	//read raster elevation
	directory_founder(tempfile, "input", case_name_slashed, "\\elevation", ".input");
	read_raster(tempfile, nX, nY, zb);

	printf("\n\n>> Geometry loaded");
	printf("\n");
	fprintf(logFile,"\n\n>> Geometry loaded");
	fprintf(logFile,"\n");



// SIMULATION SETUP
	double simTime;
	double Toutput;
	double CFL;
	double n1;

	directory_founder(tempfile, "input", case_name_slashed, "\\simulation", ".input");
	//sprintf(tempfile,"input/simulation.input");
	read_simulation_setup(tempfile, &(simTime),
		&(Toutput), &(CFL), &(n1));


	double QIN, HIN;
	double HOUT, ZSOUT;
	QIN = -1.0;
	HIN = -1.0;
	HOUT = -1.0;
	ZSOUT = -1.0;

	directory_founder(tempfile, "input", case_name_slashed, "\\simulation", ".input");
	read_boundary_configuration(tempfile,
		&(QIN), &(HIN),
		&(HOUT), &(ZSOUT));


	printf("\n   --------------------------------------------------");
	printf("\n>> Simulation setup loaded:");
	printf("\n     Simulation time: %.1lf s",simTime);
	printf("\n     CFL: %.2lf",CFL);
	printf("\n     Friction term activated: %d",friction);
	printf("\n     Friction model: %d",Fmodel);
	printf("\n     Number of cells: %d %d Total: %d",nX, nY, nCells);
	printf("\n     dx: %.2lf m", dx);
	printf("\n     Data saved each %d iterations",Niter);
	printf("\n");
	printf("\n     Fluid density: %04.0lf kg/m3",rhow);
	printf("\n     Dynamic viscosity: %.6lf Pa.s",mu);
	printf("\n     Manning coefficient: %.3lf",n1);



	fprintf(logFile,"\n\n   --------------------------------------------------");
	fprintf(logFile,"\n>> Simulation setup loaded:");
	fprintf(logFile,"\n     Simulation time: %.1lf s",simTime);
	fprintf(logFile,"\n     CFL: %.2lf",CFL);
	fprintf(logFile,"\n     Friction term activated: %d",friction);
	fprintf(logFile,"\n     Friction model: %d",Fmodel);
	fprintf(logFile,"\n     Number of cells: %d %d Total: %d",nX, nY, nCells);
	fprintf(logFile,"\n     dx: %.2lf m",dx);
	fprintf(logFile,"\n     Data saved each %d iterations",Niter);
	fprintf(logFile,"\n");
	fprintf(logFile,"\n     Fluid density: %04.0lf kg/m3",rhow);
	fprintf(logFile,"\n     Dynamic viscosity: %.6lf Pa.s",mu);
	fprintf(logFile,"\n     Manning coefficient: %.3lf",n1);


	printf("\n\n>> Simulation setup loaded");
	printf("\n");
	fprintf(logFile,"\n\n>> Simulation setup loaded");
	fprintf(logFile,"\n");


// VARIABLE INITIALIZATION
	iter = 1;					// Iteration number
	t = 0.0;					// Initial time

	h_initialize_variables( nCells, h,  qx,  qy,  ux,  uy,
        n, qBx,  qBy, Bx, By,
		DU1, DU2, DU3, DU4, DU5);

	printf("\n\n>> Flow initialization completed");
	fprintf(logFile,"\n\n>> Flow initialization completed");



/////////////////////////////////////////////////////////////////////////
/// INITIAL CONDITIONS
/////////////////////////////////////////////////////////////////////////

	//read hini
    directory_founder(tempfile, "input", case_name_slashed, "\\hini", ".input");
	read_raster(tempfile, nX, nY, h);

	//read qvxini
    directory_founder(tempfile, "input", case_name_slashed, "\\qxini", ".input");
	read_raster(tempfile, nX, nY, qx);

	//read qvyini
    directory_founder(tempfile, "input", case_name_slashed, "\\qyini", ".input");
	read_raster(tempfile, nX, nY, qy);

	//read qBxini
    directory_founder(tempfile, "input", case_name_slashed, "\\qBxini", ".input");
	read_raster(tempfile, nX, nY, qBx);

	//read qByini
	directory_founder(tempfile, "input", case_name_slashed, "\\qByini", ".input");
	read_raster(tempfile, nX, nY, qBy);


	// Flow variables calculation
	h_compute_initial_flow_variables( nCells,
		h,  qx,  qy, ux,  uy, qBx, qBy, Bx, By, n,  n1);


	printf("\n\n>> Initial conditions loaded");
	fprintf(logFile,"\n\n>> Initial conditions loaded");


	// Initial mass balance
	h_compute_water_mass( nCells, h,  area,  &(massWtn));
	h_compute_divergence_condition(nX, nY, qBx, qBy, &(divergence_loss), dx);

	printf("mass:%lf\n",massWtn);
	printf("magnetic divergence: %.3e\n", divergence_loss);
	Qin = 0.0;
	Qout = 0.0;
	Qbalance = 0.0;


	fprintf(logFile,"\n>> INITIAL FLOW MASS %.6lf m3",massWtn);
	fprintf(logFile,"\n");

	// Store initial conditions
	//Cell data
	nout = 0;

	char celldata[1024];
	sprintf(celldata, "\\archivos_vtk\\celldata%d", nout);
    directory_founder(tempfile, "outputFiles", case_name_slashed, celldata, ".vtk");
	//sprintf(filename, "outputFiles/celldata%d.vtk",nout);
	write_vtk_cells(tempfile, nX, nY, x, y,
			zb, h, ux, uy, qBx, qBy, Bx, By);
	nout += 1;



	// Discharge data
	fprintf(qFile,"%.3lf\t %.6lf\t %.6lf\t %.6lf\n",
		t,
		Qin,
		Qout,
		Qbalance);

	printf("\n\nPress INTRO key to start ...");
	printf("\n\n");
	getchar();



CPUtime = clock();
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////     TIME LOOP     /////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
while(t <= simTime) {

// TIME STEP COMPUTATION
	h_compute_flow_time_step_2D(nX, nY,
			h, ux, uy, dx,
			qBx, qBy, Bx, By,
			&dtSW);
	dt = CFL * dtSW;




// DISPLAY TIME AND TIME STEP
	if(iter%Niter == 0){
		printf("\n--------------------------------------------------------------------------------------------------------");
		printf("\n>> Time: %.3lf seg  -  Iteration %d  -  Time step %.6lf",t,iter,dt);
	}

// FLUXES CALCULATION
	// X-direction fluxes
	h_compute_x_fluxes(nX, nY,
    	h, qx, qy, ux, uy,
    	zb, n, qBx, qBy, Bx, By,
    	DU1, DU2, DU3, DU4, DU5, dx, QIN, HIN);


	// Y-direction fluxes
	h_compute_y_fluxes(nX, nY,
    	h, qx, qy, ux, uy,
    	zb, n, qBx, qBy, Bx, By,
    	DU1, DU2, DU3, DU4, DU5, dx, QIN, HIN);



// NEGATIVE WATER DEPTH CHECKING
	h_check_depth_positivity( nCells,
		h,  DU1,  dx,  &(dt));



// SHALLOW WATER CELL UPDATE
	h_update_cells_2D(nCells,
			h, qx, qy, ux, uy,
			qBx, qBy, Bx, By,
			DU1, DU2, DU3, DU4, DU5,
			dx, dt);

    #ifdef debug
        printf("\nvariables justo despues de la update");
        print_variables(nX, nY,
            h, qx, qy, ux, uy,
            qBx, qBy, Bx, By,
            DU1, DU2, DU3, DU4, DU5,
            dx, dt);
    #endif // debug


	h_wet_dry_x(nX, nY,
		h, qx, ux, zb);

	h_wet_dry_y(nX, nY,
		h, qy, uy, zb);
   // West boundary
   h_set_west_boundary(nX, nY, h, qx, qy, ux, uy, qBx, qBy, Bx, By, QIN, HIN);
   // East boundary
   h_set_east_boundary(nX, nY, h, qx, qy, ux, uy, zb, qBx, qBy, Bx, By, QIN, HIN);
   // North boundary
   h_set_north_boundary(nX, nY, h, qBx, qBy, Bx, By, qx, qy, ux, uy, QIN, HIN);
   // South boundary
   h_set_south_boundary(nX, nY, h, qBx, qBy, Bx, By, qx, qy, ux, uy, QIN, HIN);

// SIMULATION MONITORS
	massWt0 = massWtn;
	h_compute_water_mass( nCells,
		h,  dx,  &(massWtn));
	h_compute_divergence_condition(nX, nY, qBx, qBy, &(divergence_loss), dx);

	//printf("mass:%lf\n",massWtn);


	//Compute mass error
	if(massWt0 != 0.0) {
		massWerror = (massWtn - (massWt0-(Qout-Qin)*dt)) / massWt0;
		if(fabs(massWerror) < 1e-16) massWerror = 1e-16;
	} else {
		massWerror = 0.0;
	}

	Qbalance += (Qin-Qout)*dt;


	if(iter%Niter == 0){
		printf("\n\tMass error %.3e\t",massWerror);
        printf("Magnetic divergence: %.3e\t", divergence_loss);
		printf("\n\tWater discharge IN %.6lf OUT %.6lf  [m3/s]\n",Qin,Qout);
	}


	Qin = 0.0; // West inflow
    Qout = 0.0; // East outflow


// UPDATE TIME
	iter++;
	t = t+dt;

	if(dt <= 1e-9) { break;}


	//write files!!!

// DATA OUTPUT
	if(t >= nout*Toutput){
		//Cell data
		//output vtk
        sprintf(celldata, "\\archivos_vtk\\celldata%d", nout);
        directory_founder(tempfile, "outputFiles", case_name_slashed,celldata, ".vtk");
		//sprintf(filename, "outputFiles/celldata%d.vtk",nout);
		write_vtk_cells(tempfile, nX, nY, x, y,
			zb, h, ux, uy, qBx, qBy, Bx, By);
		// Discharge data
		fprintf(qFile,"%.3lf\t %.6lf\t %.6lf\t %.6lf\n",
			t,
			Qin,
			Qout,
			Qbalance);

		nout += 1;
	}

	#ifdef debug/*
        for (int j = 0; j < nY; j++) {     // Rows
            for (int i = 0; i < nX; i++) { // Columns
                int ic = IDX(i,j,nX);
                if(fabs(Bx[ic]) > tol12){
                    printf("row=%d, column=%d, t=%f, iter=%d, Bx=%.12f\n", j, i, t, iter, Bx[ic]);
                    printf("h=%f\n", h[ic]);
                    getchar();
                }
                if(fabs(By[ic]) > tol12){
                    printf("row=%d, column=%d, t=%f, iter=%d, By=%.12f\n", j, i, t, iter, By[ic]);
                    printf("h=%f\n", h[ic]);
                    getchar();
                }
            }
        }*/

        printf("\nvariables final\n");
        print_variables(nX, nY,
            h, qx, qy, ux, uy,
            qBx, qBy, Bx, By,
            DU1, DU2, DU3, DU4, DU5,
            dx, dt);
        getchar();
    #endif // debug
}
//////////////////////////////////     END TIME LOOP     //////////////////////////////////////////////////////////

sprintf(celldata, "\\archivos_vtk\\celldata%d", nout);
directory_founder(tempfile, "outputFiles", case_name_slashed, celldata, ".vtk");
write_vtk_cells(tempfile, nX, nY, x, y,
    zb, h, ux, uy, qBx, qBy, Bx, By);



CPUtime = clock() - CPUtime ;
double time = ((float)CPUtime)/CLOCKS_PER_SEC;








// DISPLAY FINAL INFORMATION
printf("\n\n>> Final Time: %.3lf seg",t);
printf("\n\n>> Computation time %.3lf seg",time);
printf("\n\n>> Simulation completed!");
printf("\n\n ");

fprintf(logFile,"\n>> Final Time: %.3lf seg",t);
fprintf(logFile,"\n\n>> Computation time %.3lf seg",time);
fprintf(logFile,"\n\n>> Simulation completed!");
fprintf(logFile,"\n\n ");






// CLOSING DATA STORAGE
fclose(logFile);
fclose(qFile);



//Free memory
/*
liberate_memory(x, y, zb, h, ux, uy, qBx, qBy,
    Bx, By, n, DU1, DU2, DU3, DU4, DU5, qx, qy);
*/

} // End of main function
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







