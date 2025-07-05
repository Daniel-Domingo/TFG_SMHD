//   Shallow Water Equations 2D Simulator (SWE2D)
//  ===========================================================================================
//   Square grid
//  ===========================================================================================
//   Version 1.0 - Mar 2025
//  ===========================================================================================
//   Computational Hydraulics Group - University of Zaragoza
//  ===========================================================================================


#include "shallow_water.h"
void directory_founder(char *tempfile, const char *header, const char *case_name_slashed, const char *file_name, const char *extension){
    sprintf(tempfile, header);
    strcat(tempfile, case_name_slashed);
    strcat(tempfile, file_name);
    strcat(tempfile, extension);
}

int read_raster_header(const char *filename, int *nX, int *nY,
	double *Xmin, double *Ymin, double *dx){

	double nodata;

	FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return -1;
    }

	// Read header
    fscanf(file, "NCOLS %d\n", nX);
    fscanf(file, "NROWS %d\n", nY);
    fscanf(file, "XLLCORNER %lf\n", Xmin);
    fscanf(file, "YLLCORNER %lf\n", Ymin);
    fscanf(file, "CELLSIZE %lf\n", dx);
    fscanf(file, "NODATA_VALUE %lf\n", &nodata);

    fclose(file);
    return 0;

}

int read_raster(const char *filename, int nX, int nY,
	double *data){

	char buffer[256];

	FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return -1;
    }

	// Skip header lines
	for (int i = 0; i < 6; i++) {
        if (!fgets(buffer, sizeof(buffer), file)) {
            perror("Error reading header");
            fclose(file);
            return -1;
        }
    }

	// Read matrix
	for (int j = 0; j < nY; j++) {     // Rows
		for (int i = 0; i < nX; i++) { // Columns
			if (fscanf(file, "%lf", &data[IDX(i, j, nX)]) != 1) {
				perror("Error reading raster data");
				fclose(file);
				return -1;
			}
		}
}


    fclose(file);
    return 0;

}


int read_simulation_setup(char *tempfile, double *simTime,
	double *Toutput, double *CFL, double *n1 ){

	FILE *fp;
	char line[1024];
	const char* item;
	char* position;
	char* endPtr;

	//Check file
	fp = fopen(tempfile,"r");
	if(!fp){
		printf("File %s not found",tempfile);
		return 0;
	}

	//Initial condition for bed surface
    while (fgets(line, sizeof(line), fp)) {

		item = "simTime:";
        position = strstr(line, item);
        if(position){
			position += strlen(item);
			(*simTime) = (double)strtod(position, &endPtr); //double
        }

		item = "Toutput:";
        position = strstr(line, item);
        if(position){
			position += strlen(item);
			(*Toutput) = (double)strtod(position, &endPtr); //double
        }

		item = "CFL:";
        position = strstr(line, item);
        if(position){
			position += strlen(item);
			(*CFL) = (double)strtod(position, &endPtr); //double
        }

		item = "n1:";
        position = strstr(line, item);
        if(position){
			position += strlen(item);
			(*n1) = (double)strtod(position, &endPtr); //double
        }

    }
	fclose(fp);

	return 1;
}



int read_boundary_configuration(char *tempfile,
	double *QIN, double *HIN,
	double *HOUT, double *ZSOUT){

	FILE *fp;
	char line[1024];
	const char* item;
	char* position;
	char* endPtr;

	//Check file
	fp = fopen(tempfile,"r");
	if(!fp){
		printf("File %s not found",tempfile);
		return 0;
	}

	//Initial condition for bed surface
    while (fgets(line, sizeof(line), fp)) {

		item = "QIN:";
        position = strstr(line, item);
        if(position){
			position += strlen(item);
			(*QIN) = (double)strtod(position, &endPtr); //double
        }

		item = "HIN:";
        position = strstr(line, item);
        if(position){
			position += strlen(item);
			(*HIN) = (double)strtod(position, &endPtr); //double
        }

		item = "HOUT:";
        position = strstr(line, item);
        if(position){
			position += strlen(item);
			(*HOUT) = (double)strtod(position, &endPtr); //double
        }

		item = "ZSOUT:";
        position = strstr(line, item);
        if(position){
			position += strlen(item);
			(*ZSOUT) = (double)strtod(position, &endPtr); //double
        }

    }
	fclose(fp);

	return 1;
}



void h_initialize_variables(int nCells,
	double *h, double *qx, double *qy, double *ux, double *uy,
	double *n, double *qBx, double *qBy, double *Bx, double *By,
	double *DU1, double *DU2, double *DU3, double *DU4, double *DU5){

    	for(int ic = 0; ic < nCells; ic++){
		h[ic] = 0.0;
		qx[ic] = tol12;
		qy[ic] = tol12;
		ux[ic] = tol12;
		uy[ic] = tol12;
        qBx[ic] = tol12;
		qBy[ic] = tol12;
		Bx[ic] = tol12;
		By[ic] = tol12;

		n[ic] = tol12;

		DU1[ic] = 0.0;
		DU2[ic] = 0.0;
		DU3[ic] = 0.0;
		DU4[ic] = 0.0;
		DU5[ic] = 0.0;

	}

}


int h_compute_initial_flow_variables(int nCells,
	double *h, double *qx, double *qy, double *ux, double *uy,
	double *qBx, double *qBy, double *Bx, double *By,
	double *n, double n1){

    int ic;
	for(ic=0; ic<nCells; ic++){
		// Minimum depth control
		if(h[ic] >= hmin){
			// Minimum depth control
			ux[ic]=qx[ic]/h[ic];
			uy[ic]=qy[ic]/h[ic];
			Bx[ic]=qBx[ic]/h[ic];
			By[ic]=qBy[ic]/h[ic];

		}else{
			qx[ic]=0.0;
			qy[ic]=0.0;
			ux[ic]=0.0;
			uy[ic]=0.0;
			qBx[ic]=0.0;
			qBy[ic]=0.0;
			Bx[ic]=0.0;
			By[ic]=0.0;

		}
		// Cell Manning Coefficient
        n[ic] = n1;

	}

    return 1;

}



void h_compute_water_mass(int nCells,
	double *h, double area, double *massWt){

	(*massWt) = 0.0;
	for(int ic=0; ic<nCells; ic++){
		(*massWt) += h[ic]*area;
	}

}

///Ahora no va a funcionar (bien) si se ponen condiciones periódicas
void h_compute_divergence_condition(int nX, int nY, double *qBx, double *qBy,
    double *divergence_loss, double dx){
    double dif_aux;
    (*divergence_loss) = 0;
	for (int j = 0; j < nY-1; j++) {     // Rows
		for (int i = 0; i < nX-1; i++) { // Columns
            int origin = IDX(i,j,nX);   //tomamos una celda de referencia y aproximamos la condición de diverngencia con respecto de la que tiene a la derecha y encima
            int right = IDX(i+1,j,nX);
			int top = IDX(i,j+1,nX);
/*
            if (qBx[origin] != 0)
                printf("Row: %d, Column: %d, qBx: %.16f\n", j, i, qBx[origin]);
            if (qBy[origin] != 0)
                printf("Row: %d, Column: %d, qBy: %.16f\n", j, i, qBy[origin]);*/
			//Se pone el fabs para evitar que se anulen contribuciones negativas con positivas
			dif_aux = fabs((qBx[right]-qBx[origin] + qBy[top]-qBy[origin])/dx);
            (*divergence_loss) += dif_aux;
            //Esto es para hacer pruebas
            /*
            if(dif_aux != 0.0){
                printf("Row: %d, Column: %d, Divergence: %f\n", j, i, dif_aux);
            }
            */
		}
    }
}


///El tema de las normales no va a funcionar, aka, por ahora hay que usar mallas cartesianas
///No tengo del todo claro como acaba de furrular el tema de las velocidades de onda. Creo que está bien pero habrá que revisarlo
 void h_compute_flow_time_step_2D(int nX, int nY,
		double *h, double *ux, double *uy, double dx,
		double *qBx, double *qBy, double *Bx, double *By,
double *dtSW) {

    double uxROE, uyROE, cROE, BxROE, ByROE;
    double dt = 1e6;
    //Dejo por aquí una versión que funcione en el caso de que se quiera poner una malla no cartesiana
    /*
    double normalX, normalY;
    normalX = 1.0;
    normalY = 0.0;
    */

    // Process horizontal interfaces
	for (int j = 0; j < nY; j++) {     // Rows
		for (int i = 0; i < nX-1; i++) { // Columns
            int left = IDX(i,j,nX);
            int right = IDX(i+1,j,nX);
			if (h[left] > tol9 || h[right] > tol9){
                BxROE = (sqrt(h[left])*Bx[left]+sqrt(h[right])*Bx[right]) / (sqrt(h[left])+sqrt(h[right]));
                //Dejo por aquí una versión que funcione en el caso de que se quiera poner una malla no cartesiana
                //ByROE = (sqrt(h[left])*By[left]+sqrt(h[right])*By[right]) / (sqrt(h[left])+sqrt(h[right]));
				//cROE = sqrt(g*0.5*(h[left]+h[right])+(BxROE*normalX+ByROE*normalY)^2);
				cROE = sqrt(g*0.5*(h[left]+h[right])+BxROE*BxROE);
				uxROE = (sqrt(h[left])*ux[left]+sqrt(h[right])*ux[right]) / (sqrt(h[left])+sqrt(h[right]));
				if(fabs(uxROE)+cROE > tol9) {
					dt = MIN(dt, dx/(fabs(uxROE)+cROE));
				}
			}
        }
    }

    /*
    normalX = 0.0;
    normalY = 1.0;
    */

    // Process vertical interfaces
	for (int j = 0; j < nY-1; j++) {     // Rows
		for (int i = 0; i < nX; i++) { // Columns
            int bottom =IDX(i,j,nX);
            int top = IDX(i,j+1,nX);
			if (h[bottom] > tol9 || h[top] > tol9){
                ByROE = (sqrt(h[bottom])*By[bottom]+sqrt(h[top])*By[top]) / (sqrt(h[bottom])+sqrt(h[top]));
                //BxROE = (sqrt(h[left])*Bx[left]+sqrt(h[right])*Bx[right]) / (sqrt(h[left])+sqrt(h[right]));
                //Dejo por aquí una versión que funcione en el caso de que se quiera poner una malla no cartesiana
				//cROE = sqrt(g*0.5*(h[left]+h[right])+(BxROE*normalX+ByROE*normalY)^2);
				cROE = sqrt(g*0.5*(h[bottom]+h[top])+ByROE*ByROE);
				uyROE = (sqrt(h[top])*uy[top]+sqrt(h[bottom])*uy[bottom]) / (sqrt(h[top])+sqrt(h[bottom]));
            	if(fabs(uyROE)+cROE > tol9) {
             		dt = MIN(dt, dx/(fabs(uyROE)+cROE));
   		        }
			}
        }
    }

    (*dtSW) = dt;
}


//---------------Función: h_compute_x_fluxes------------------
///No van a estar bien las condiciones de positividad, cuidado con eso
//Estas funciones se pueden limpiar mucho para que rulen mejor; p.ej: todo lo relacionado a la quinta onda (al menos por ahora) no sirve de nada
void h_compute_x_step(int nX, int nY,
    double *h, double *qx, double *qy, double *ux, double *uy,
	double *zb, double *n, double *qBx, double *qBy, double *Bx, double *By,
	double *DU1, double *DU2, double *DU3, double *DU4, double *DU5,
    double dx, int left, int right){

	double uxROE, uyROE, cROE, BxROE, ByROE;       		// ROE averaged variables
	double lambda1, lambda2, lambda3, lambda4, lambda5;				// Eigenvalues
	double e1[5], e2[5], e3[5], e4[5], e5[5];					// Eigenvectors

	double deltah0, deltaqx0, deltaqy0, deltaqBx0, deltaqBy0;			// Spatial increments of the conserved variables
	double alpha1, alpha2, alpha3, alpha4, alpha5;					// Wave-strenghts

	double So, Fp, b1P;			// Hydrostatic pressure force at bed
	double nav;								// 2D averaged Manning Coefficient
	double Sf, Ff, b1F;				// Fupwind force at bed
	double beta1, beta2, beta3, beta4, beta5;		// Source term coefficients

	double hstar,hx,hxx,Qx,Qstart;			// Non-negative area values control - Beta limitation

	double lb1l, lb1r, lb3l, lb3r, lambda1plus, lambda1minus, lambda3plus, lambda3minus;		// Entropy correction

	double hav;

	double un, normalX, normalY;

    //El quinto autovector es super chungo, vamos a hacer un par de apaños para simplificar:
    double B, W, a2, u, v;
    double B2, B3, B4, u2, u3, u4;
    double denominator_e5, u_minus_c, u_plus_c, B_minus_u, B_plus_u, e5_3_aux;

	//normal direction in X
	normalX=1.0;
	normalY=0.0;
    if(h[left] >= tol9 || h[right] >= tol9){			// Wet wall

        // Averaged quantities at the walls
        hav = (h[left]+h[right])/2.;

        // ROE averaged variables at the walls
        // Note about the notation:
        // ux=u (x-axis component of the velocity)
        // uy=v (y-axis component of the velocity)
        uxROE = (sqrt(h[left])*ux[left]+sqrt(h[right])*ux[right]) / (sqrt(h[left])+sqrt(h[right]));
        uyROE = (sqrt(h[left])*uy[left]+sqrt(h[right])*uy[right]) / (sqrt(h[left])+sqrt(h[right]));
        BxROE = (sqrt(h[left])*Bx[left]+sqrt(h[right])*Bx[right]) / (sqrt(h[left])+sqrt(h[right]));
        ByROE = (sqrt(h[left])*By[left]+sqrt(h[right])*By[right]) / (sqrt(h[left])+sqrt(h[right]));
        cROE = sqrt(g*hav + pow(BxROE*normalX + ByROE*normalY, 2));   //In magnetohydrodynamics, cROE is modified


        //normal projection of the velocity
        un = uxROE*normalX + uyROE*normalY;

        // Eigenvalues
        lb1l = ux[left]*normalX + uy[left]*normalY-sqrt(g*h[left]);
        lb1r = ux[right]*normalX + uy[right]*normalY-sqrt(g*h[right]);
        lambda1 = un - cROE;
        lambda2 = (uxROE-BxROE)*normalX + (uyROE-ByROE)*normalY;
        lb3l = ux[left]*normalX + uy[left]*normalY+sqrt(g*h[left]);
        lb3r = ux[right]*normalX + uy[right]*normalY+sqrt(g*h[right]);
        lambda3 = (uxROE+BxROE)*normalX + (uyROE+ByROE)*normalY;
        lambda4 = un + cROE;
        lambda5 = 0.;

        // Eigenvectors
        e1[0] = 1;    e1[1] = (uxROE-cROE);    e1[2] = uyROE;     e1[3] = 0.;     e1[4] = ByROE;
        e2[0] = 0.;   e2[1] = 0.;              e2[2] = cROE;      e2[3] = 0.;     e2[4] = cROE;
        e3[0] = 0.;   e3[1] = 0.;              e3[2] = cROE;      e3[3] = 0.;     e3[4] = -cROE;
        e4[0] = 1;    e4[1] = (uxROE+cROE);    e4[2] = uyROE;     e4[3] = 0.;     e4[4] = ByROE;
        //El quinto autovector es super chungo, vamos a hacer un par de apaños para simplificar:
        B=BxROE; W=ByROE; a2=g*hav; u=uxROE; v=uyROE;
        B2=B*B; B3=B*B*B; B4=B*B*B*B;
        u2=u*u; u3=u*u*u; u4=u*u*u*u;
        denominator_e5 = B3*W-B*W*a2-B*W*u2+u3*v-B2*u*v-a2*u*v;
        e5[0] = (2*B3-2*B*u2);
        e5[1] = 0.;
        e5[2] = (W*u3-B2*W*u-W*a2*u+B3*v-B*a2*v-B*u2*v);
        e5[3] = (B4+B2*a2+u4-2*B2*u2-a2*u2);
        e5[4] = denominator_e5;


        // Wave-strenght coefficients
        deltah0 = h[right]-h[left];
        deltaqx0 = qx[right]-qx[left];
        deltaqy0 = qy[right]-qy[left];
        deltaqBx0 = qBx[right]-qBx[left];
        deltaqBy0 = qBy[right]-qBy[left];

        //Controles de división por 0
        B_minus_u = (B-u);
        B_plus_u = (B+u);
        u_minus_c = u-cROE;
        u_plus_c = u+cROE;
        e5_3_aux = e5[3];
        if(fabs(B_minus_u)<tol9){B_minus_u = tol9;}
        if(fabs(B_plus_u)<tol9){B_plus_u = tol9;}
        if(fabs(u_minus_c)<tol9){u_minus_c = tol9;}
        if(fabs(u_plus_c)<tol9){u_plus_c = tol9;}
        if(fabs(cROE)<tol9){cROE = tol9;}// printf("cROE=%f", cROE);}
        if(fabs(e5_3_aux)<tol9){e5_3_aux = tol9;}

        //if(fabs(e5[3])>tol9){printf("e5[3]=%.9f\n", e5[3]);}

        alpha1 = (deltah0*(u+cROE)/2-deltaqx0/2+deltaqBx0*B/(u_minus_c))/cROE;
        alpha2 = (-deltah0*(W+v)+deltaqy0+deltaqBx0*(W+v)/(B_minus_u)+deltaqBy0)/(2*cROE);
        alpha3 = (deltah0*(W-v)+deltaqy0+deltaqBx0*(W-v)/(B_plus_u)-deltaqBy0)/(2*cROE);
        alpha4 = (deltah0*(cROE-u)/2+deltaqx0/2-deltaqBx0*B/(u_plus_c))/cROE;
        alpha5 = deltaqBx0/e5_3_aux;

        #ifdef debug
            print_eigen(left, right, e1, e2, e3, e4, e5, alpha1, alpha2, alpha3, alpha4, alpha5);
        #endif // debug




        // Source term coefficients
        //Bed slope differential formulation

        #if bed_slope_integral
            double deltaz = zb[right] - zb[left];
            double l1 = zb[left] + h[left];
            double l2 = zb[right] + h[right];
            double dzp = deltaz;
            double hp;

            if (deltaz >= 0.0) {
                hp = h[left];
                if (l1 < zb[right]) {
                    dzp = h[left];
                }
            } else {
                hp = h[right];
                if (l2 < zb[left]) {
                    dzp = -h[right];
                }
            }

            // Calculate the bed slope source term
            b1P= (1./(2.*cROE)) * g * (hp - 0.5 * fabs(dzp)) * dzp;
        #else
            So = -(zb[right]-zb[left])/dx;
            Fp = g*hav*So*dx;
            b1P = (-1./(2.*cROE)) * Fp;
        #endif

        //friction
        nav = 0.5*(n[left]+n[right]);
        if(h[left] < tol9 || h[right] < tol9){
            nav=0.0;
        }
        if(friction == 1) {
            // Manning
            if(Fmodel == 1) { Sf = nav*nav*un*sqrt(uxROE*uxROE + uyROE*uyROE)/pow(hav,4./3.); } // Unit formulation
            //Viscous
            if(Fmodel == 2) { Sf = 3.*mu*un / (rhow*g*hav); } // Unit formulation
        } else {
            Sf = 0.0;
        }

        Ff = -g*hav*Sf*dx;

        b1F = (-1./(2.*cROE)) * Ff;

        // Friction fix
        Qstart = (qx[left]+alpha1*e1[1])*normalX + (qy[left]+alpha1*e1[2])*normalY - b1P;
        Qx = Qstart - b1F;
        if(fabs(Qstart)<tol9) Qstart=0.0;
        if(fabs(Qx)<tol9) Qx=0.0;

        if(Qx*Qstart < 0.0){
            b1F = Qstart;
        }

        // beta coefficient
        beta1 = b1P + b1F;
        beta2= 0.0;
        beta3 = -beta1;
        beta4 = 0.0;
        beta5 = 0.0;

        // Positivity fix
        hstar = h[left]+alpha1;
        if(hstar<tol9) hstar = 0.0;

        if( (h[left]>0.0 && h[right]>0.0) &&
            (hstar>0.0) &&
            (lambda1*lambda3 < 0.0) ){

            hx = h[left]+alpha1-beta1/lambda1;
            hxx = h[right]-alpha3+beta3/lambda3;
            if(fabs(hx)<tol9) hx = 0.0;
            if(fabs(hxx)<tol9) hxx = 0.0;

            if(hx < 0.0){
                beta1 = hstar*lambda1;
                beta3 = -beta1;
            }

            if(hxx < 0.0){
                beta1 = hstar*lambda3;
                beta3 = -beta1;
            }
        }


        // Update contributions
        hx = h[left]+alpha1-beta1/lambda1;
        hxx = h[right]-alpha3+beta3/lambda3;
        if(fabs(hx)<tol9) hx = 0.0;
        if(fabs(hxx)<tol9) hxx = 0.0;

        // First wave
        if(h[left]<tol9 && hx < 0.0){	// dry-wet wall
            DU1[right] += (lambda1*alpha1-beta1)*e1[0];
            DU2[right] += 0.0;
            DU3[right] += 0.0;
            DU4[right] += 0.0;
            DU5[right] += 0.0;

        } else if(h[right]<tol9 && hxx < 0.0){ // wet-dry wall
            DU1[left] += (lambda1*alpha1-beta1)*e1[0];
            DU2[left] += 0.0;
            DU3[left] += 0.0;
            DU4[left] += 0.0;
            DU5[left] += 0.0;

        } else if(lb1l < 0.0 && lb1r > 0.0){
            lambda1minus = lb1l*(lb1r-lambda1)/(lb1r-lb1l);
            DU1[left] += (lambda1minus*alpha1-beta1)*e1[0];
            DU2[left] += (lambda1minus*alpha1-beta1)*e1[1];
            DU3[left] += (lambda1minus*alpha1-beta1)*e1[2];
            DU4[left] += (lambda1minus*alpha1-beta1)*e1[3];
            DU5[left] += (lambda1minus*alpha1-beta1)*e1[4];
            lambda1plus = lb1r*(lambda1-lb1l)/(lb1r-lb1l);
            DU1[right] += (lambda1plus*alpha1)*e1[0];
            DU2[right] += (lambda1plus*alpha1)*e1[1];
            DU3[right] += (lambda1plus*alpha1)*e1[2];
            DU4[right] += (lambda1plus*alpha1)*e1[3];
            DU5[right] += (lambda1plus*alpha1)*e1[4];
            //ni flower que hacer con los DU4 y DU5

        } else if(lambda1 < 0.0){
            DU1[left] += (lambda1*alpha1-beta1)*e1[0];
            DU2[left] += (lambda1*alpha1-beta1)*e1[1];
            DU3[left] += (lambda1*alpha1-beta1)*e1[2];
            DU4[left] += (lambda1*alpha1-beta1)*e1[3];
            DU5[left] += (lambda1*alpha1-beta1)*e1[4];

        } else if(lambda1 >= 0.0){
            DU1[right] += (lambda1*alpha1-beta1)*e1[0];
            DU2[right] += (lambda1*alpha1-beta1)*e1[1];
            DU3[right] += (lambda1*alpha1-beta1)*e1[2];
            DU4[right] += (lambda1*alpha1-beta1)*e1[3];
            DU5[right] += (lambda1*alpha1-beta1)*e1[4];
        }

        // Second wave
        if(h[left]<tol9 && hx < 0.0){ 		// dry-wet wall
            DU1[right] += (lambda2*alpha2-beta2)*e2[0];
            DU2[right] += 0.0;
            DU3[right] += 0.0;
            DU4[right] += 0.0;
            DU5[right] += 0.0;

        } else if(h[right]<tol9 && hxx < 0.0){	// wet-dry wall
            DU1[left] += (lambda2*alpha2-beta2)*e2[0];
            DU2[left] += 0.0;
            DU3[left] += 0.0;
            DU4[left] += 0.0;
            DU5[left] += 0.0;

        } else if(lambda2 < 0.0){
            DU1[left] += (lambda2*alpha2-beta2)*e2[0];
            DU2[left] += (lambda2*alpha2-beta2)*e2[1];
            DU3[left] += (lambda2*alpha2-beta2)*e2[2];
            DU4[left] += (lambda2*alpha2-beta2)*e2[3];
            DU5[left] += (lambda2*alpha2-beta2)*e2[4];

        } else if(lambda2 >= 0.0){
            DU1[right] += (lambda2*alpha2-beta2)*e2[0];
            DU2[right] += (lambda2*alpha2-beta2)*e2[1];
            DU3[right] += (lambda2*alpha2-beta2)*e2[2];
            DU4[right] += (lambda2*alpha2-beta2)*e2[3];
            DU5[right] += (lambda2*alpha2-beta2)*e2[4];
        }

        // Third wave
        if(h[left]<tol9 && hx < 0.0){ 		// dry-wet wall
            DU1[right] += (lambda3*alpha3-beta3)*e3[0];
            DU2[right] += 0.0;
            DU3[right] += 0.0;
            DU4[right] += 0.0;
            DU5[right] += 0.0;

        } else if(h[right]<tol9 && hxx < 0.0){	// wet-dry wall
            DU1[left] += (lambda3*alpha3-beta3)*e3[0];
            DU2[left] += 0.0;
            DU3[left] += 0.0;
            DU4[left] += 0.0;
            DU5[left] += 0.0;

        } else if(lb3l < 0.0 && lb3r > 0.0){
            lambda3minus=lb3l*(lb3r-lambda3)/(lb3r-lb3l);
            DU1[left] += (lambda3minus*alpha3)*e3[0];
            DU2[left] += (lambda3minus*alpha3)*e3[1];
            DU3[left] += (lambda3minus*alpha3)*e3[2];
            DU4[left] += (lambda3minus*alpha3)*e3[3];
            DU5[left] += (lambda3minus*alpha3)*e3[4];
            lambda3plus=lb3r*(lambda3-lb3l)/(lb3r-lb3l);
            DU1[right] += (lambda3plus*alpha3-beta3)*e3[0];
            DU2[right] += (lambda3plus*alpha3-beta3)*e3[1];
            DU3[right] += (lambda3plus*alpha3-beta3)*e3[2];
            DU4[right] += (lambda3plus*alpha3-beta3)*e3[3];
            DU5[right] += (lambda3plus*alpha3-beta3)*e3[4];
            //ni flower que hacer con los DU4 y DU5

        } else if(lambda3 < 0.0){
            DU1[left] += (lambda3*alpha3-beta3)*e3[0];
            DU2[left] += (lambda3*alpha3-beta3)*e3[1];
            DU3[left] += (lambda3*alpha3-beta3)*e3[2];
            DU4[left] += (lambda3*alpha3-beta3)*e3[3];
            DU5[left] += (lambda3*alpha3-beta3)*e3[4];

        } else if(lambda3 >= 0.0){
            DU1[right] += (lambda3*alpha3-beta3)*e3[0];
            DU2[right] += (lambda3*alpha3-beta3)*e3[1];
            DU3[right] += (lambda3*alpha3-beta3)*e3[2];
            DU4[right] += (lambda3*alpha3-beta3)*e3[3];
            DU5[right] += (lambda3*alpha3-beta3)*e3[4];
        }

        // Fourth wave
        if(h[left]<tol9 && hx < 0.0){ 		// dry-wet wall
            DU1[right] += (lambda4*alpha4-beta4)*e4[0];
            DU2[right] += 0.0;
            DU3[right] += 0.0;
            DU4[right] += 0.0;
            DU5[right] += 0.0;

        } else if(h[right]<tol9 && hxx < 0.0){	// wet-dry wall
            DU1[left] += (lambda4*alpha4-beta4)*e4[0];
            DU2[left] += 0.0;
            DU3[left] += 0.0;
            DU4[left] += 0.0;
            DU5[left] += 0.0;

        } else if(lambda4 < 0.0){
            DU1[left] += (lambda4*alpha4-beta4)*e4[0];
            DU2[left] += (lambda4*alpha4-beta4)*e4[1];
            DU3[left] += (lambda4*alpha4-beta4)*e4[2];
            DU4[left] += (lambda4*alpha4-beta4)*e4[3];
            DU5[left] += (lambda4*alpha4-beta4)*e4[4];

        } else if(lambda4 >= 0.0){
            DU1[right] += (lambda4*alpha4-beta4)*e4[0];
            DU2[right] += (lambda4*alpha4-beta4)*e4[1];
            DU3[right] += (lambda4*alpha4-beta4)*e4[2];
            DU4[right] += (lambda4*alpha4-beta4)*e4[3];
            DU5[right] += (lambda4*alpha4-beta4)*e4[4];
        }

        // Fifth wave
        if(h[left]<tol9 && hx < 0.0){ 		// dry-wet wall
            DU1[right] += (lambda5*alpha5-beta5)*e5[0];
            DU2[right] += 0.0;
            DU3[right] += 0.0;
            DU4[right] += 0.0;
            DU5[right] += 0.0;

        } else if(h[right]<tol9 && hxx < 0.0){	// wet-dry wall
            DU1[left] += (lambda5*alpha5-beta5)*e5[0];
            DU2[left] += 0.0;
            DU3[left] += 0.0;
            DU4[left] += 0.0;
            DU5[left] += 0.0;

        } else if(lambda5 < 0.0){
            DU1[left] += (lambda5*alpha5-beta5)*e5[0];
            DU2[left] += (lambda5*alpha5-beta5)*e5[1];
            DU3[left] += (lambda5*alpha5-beta5)*e5[2];
            DU4[left] += (lambda5*alpha5-beta5)*e5[3];
            DU5[left] += (lambda5*alpha5-beta5)*e5[4];

        } else if(lambda5 >= 0.0){
            DU1[right] += (lambda5*alpha5-beta5)*e5[0];
            DU2[right] += (lambda5*alpha5-beta5)*e5[1];
            DU3[right] += (lambda5*alpha5-beta5)*e5[2];
            DU4[right] += (lambda5*alpha5-beta5)*e5[3];
            DU5[right] += (lambda5*alpha5-beta5)*e5[4];
        }

    } // End wet walls
}

void h_compute_y_step(int nX, int nY,
	double *h, double *qx, double *qy, double *ux, double *uy,
	double *zb, double *n, double *qBx, double *qBy, double *Bx, double *By,
	double *DU1, double *DU2, double *DU3, double *DU4, double *DU5,
	double dx, int top, int bottom){


	double uxROE, uyROE, cROE, BxROE, ByROE;       		// ROE averaged variables
	double lambda1, lambda2, lambda3, lambda4, lambda5;				// Eigenvalues
	double e1[5], e2[5], e3[5], e4[5], e5[5];					// Eigenvectors

	double deltah0, deltaqx0, deltaqy0, deltaqBx0, deltaqBy0;			// Spatial increments of the conserved variables
	double alpha1, alpha2, alpha3, alpha4, alpha5;					// Wave-strenghts

	double So, Fp, b1P;			// Hydrostatic pressure force at bed
	double nav;								// 2D averaged Manning Coefficient
	double Sf, Ff, b1F;				// Fupwind force at bed
	double beta1, beta2, beta3, beta4, beta5;		// Source term coefficients

	double hstar,hx,hxx,Qx,Qstart;			// Non-negative area values control - Beta limitation

	double lb1l, lb1r, lb3l, lb3r, lambda1plus, lambda1minus, lambda3plus, lambda3minus;		// Entropy correction

	double hav;

	double un, normalX, normalY;

    //El quinto autovector es super chungo, vamos a hacer un par de apaños para simplificar:
    double B, W, a2, u, v;
    double W2, W3, W4, v2, v3, v4;
    double denominator_e5, v_minus_c, v_plus_c, W_minus_v, W_plus_v;

	//normal direction in Y
	normalX=0.0;
	normalY=1.0;


    if(h[bottom] >= tol9 || h[top] >= tol9){			// Wet wall

        // Averaged quantities at the walls
        hav = (h[bottom]+h[top])/2.;

        // ROE averaged variables at the walls
        uxROE = (sqrt(h[bottom])*ux[bottom]+sqrt(h[top])*ux[top]) / (sqrt(h[bottom])+sqrt(h[top]));
        uyROE = (sqrt(h[bottom])*uy[bottom]+sqrt(h[top])*uy[top]) / (sqrt(h[bottom])+sqrt(h[top]));
        BxROE = (sqrt(h[bottom])*Bx[bottom]+sqrt(h[top])*Bx[top]) / (sqrt(h[bottom])+sqrt(h[top]));
        ByROE = (sqrt(h[bottom])*By[bottom]+sqrt(h[top])*By[top]) / (sqrt(h[bottom])+sqrt(h[top]));
        cROE = sqrt(g*hav + pow(BxROE*normalX + ByROE*normalY, 2));   //In magnetohydrodynamics, cROE is modified


        //normal projection of the velocity
        un = uxROE*normalX + uyROE*normalY;

        // Eigenvalues
        lb1l = ux[bottom]*normalX + uy[bottom]*normalY-sqrt(g*h[bottom]);
        lb1r = ux[top]*normalX + uy[top]*normalY-sqrt(g*h[top]);
        lambda1 = un - cROE;
        lambda2 = (uxROE-BxROE)*normalX + (uyROE-ByROE)*normalY;
        lb3l = ux[bottom]*normalX + uy[bottom]*normalY+sqrt(g*h[bottom]);
        lb3r = ux[top]*normalX + uy[top]*normalY+sqrt(g*h[top]);
        lambda3 = (uxROE+BxROE)*normalX + (uyROE+ByROE)*normalY;
        lambda4 = un + cROE;
        lambda5 = 0.;




        // Eigenvectors
        e1[0] = 1;    e1[1] = uxROE;     e1[2] = (uyROE-cROE);     e1[3] = BxROE;    e1[4] = 0.;
        e2[0] = 0.;   e2[1] = -cROE;     e2[2] = 0.;               e2[3] = -cROE;    e2[4] = 0.;
        e3[0] = 0.;   e3[1] = -cROE;     e3[2] = 0.;               e3[3] = cROE;     e3[4] = 0.;
        e4[0] = 1;    e4[1] = uxROE;     e4[2] = (uyROE+cROE);     e4[3] = BxROE;    e4[4] = 0.;
        //El quinto autovector es super chungo, vamos a hacer un par de apaños para simplificar:
        B=BxROE; W=ByROE; a2=g*hav; u=uxROE; v=uyROE;
        W2=W*W; W3=W*W*W; W4=W*W*W*W;
        v2=v*v; v3=v*v*v; v4=v*v*v*v;
        denominator_e5 = W4+W2*a2+v4-2*W2*v2-a2*v2;
        e5[0] = 2*W*denominator_e5/(W2+a2-v2);
        e5[1] = (3*W3*u+W*a2*u+B*v3-3*W*u*v2-B*W2*v2-B*a2*v);
        e5[2] = 0.;
        e5[3] = B*W3-B*W*a2-u*v3-B*W*v2+W2*u*v+a2*u*v;
        e5[4] = denominator_e5;


        // Wave-strenght coefficients
        deltah0 = h[top]-h[bottom];
        deltaqx0 = qx[top]-qx[bottom];
        deltaqy0 = qy[top]-qy[bottom];
        deltaqBx0 = qBx[top]-qBx[bottom];
        deltaqBy0 = qBy[top]-qBy[bottom];

        //Controles de división por 0
        W_minus_v = (W-v);
        W_plus_v = (W+v);
        v_minus_c = v-cROE;
        v_plus_c = v+cROE;
        if(fabs(W_minus_v) < tol9) W_minus_v = tol9;
        if(fabs(W_plus_v) < tol9) W_plus_v = tol9;
        if(fabs(v_minus_c) < tol9) v_minus_c = tol9;
        if(fabs(v_plus_c) < tol9) v_plus_c = tol9;
        if(fabs(cROE) < tol9) cROE = tol9;
        if(fabs(denominator_e5) < tol9) denominator_e5 = tol9;

        alpha1 = (deltah0*(v+cROE)/2-deltaqy0/2+deltaqBy0*W/(v_minus_c))/cROE;
        alpha2 = (deltah0*(B+u)-deltaqx0-deltaqBx0-deltaqBy0*(B+u)/(W_minus_v))/(2*cROE);
        alpha3 = (-deltah0*(B-u)-deltaqx0+deltaqBx0-deltaqBy0*(B-u)/(W_plus_v))/(2*cROE);
        alpha4 = (deltah0*(-v+cROE)/2+deltaqy0/2-deltaqBy0*W/(v_plus_c))/cROE;
        alpha5 = deltaqBy0/denominator_e5;

        #ifdef debug
            print_eigen(bottom, top, e1, e2, e3, e4, e5, alpha1, alpha2, alpha3, alpha4, alpha5);
        #endif // debug




        // Source term coefficients

        #if bed_slope_integral
            double deltaz = zb[top] - zb[bottom];
            double l1 = zb[bottom] + h[bottom];
            double l2 = zb[top] + h[top];
            double dzp = deltaz;
            double hp;

            if (deltaz >= 0.0) {
                hp = h[bottom];
                if (l1 < zb[top]) {
                    dzp = h[bottom];
                }
            } else {
                hp = h[top];
                if (l2 < zb[bottom]) {
                    dzp = -h[top];
                }
            }

            // Calculate the bed slope source term
            b1P= (1./(2.*cROE)) * g * (hp - 0.5 * fabs(dzp)) * dzp;
        #else
            So = -(zb[top]-zb[bottom])/dx;
            Fp = g*hav*So*dx;
            b1P = (-1./(2.*cROE)) * Fp;
        #endif


        nav = 0.5*(n[bottom]+n[top]);
        if(h[bottom] < tol9 || h[top] < tol9){
            nav=0.0;
        }
        if(friction == 1) {
            // Manning
            if(Fmodel == 1) { Sf = nav*nav*un*sqrt(uxROE*uxROE + uyROE*uyROE)/pow(hav,4./3.); } // Unit formulation
            //Viscous
            if(Fmodel == 2) { Sf = 3.*mu*un / (rhow*g*hav); } // Unit formulation
        } else {
            Sf = 0.0;
        }

        Ff = -g*hav*Sf*dx;

        b1F = (-1./(2.*cROE)) * Ff;

        // Friction fix
        Qstart = (qx[bottom]+alpha1*e1[1])*normalX + (qy[bottom]+alpha1*e1[2])*normalY - b1P;
        Qx = Qstart - b1F;
        if(fabs(Qstart)<tol9) Qstart=0.0;
        if(fabs(Qx)<tol9) Qx=0.0;

        if(Qx*Qstart < 0.0){
            b1F = Qstart;
        }

        // beta coefficients
        beta1 = b1P + b1F;
        beta2 = 0.0;
        beta3 = -beta1;
        beta4 = 0.0;
        beta5 = 0.0;

        // Positivity fix
        hstar = h[bottom]+alpha1;
        if(hstar<tol9) hstar = 0.0;

        if( (h[bottom]>0.0 && h[top]>0.0) &&
            (hstar>0.0) &&
            (lambda1*lambda3 < 0.0) ){

            hx = h[bottom]+alpha1-beta1/lambda1;
            hxx = h[top]-alpha3+beta3/lambda3;
            if(fabs(hx)<tol9) hx = 0.0;
            if(fabs(hxx)<tol9) hxx = 0.0;

            if(hx < 0.0){
                beta1 = hstar*lambda1;
                beta3 = -beta1;
            }

            if(hxx < 0.0){
                beta1 = hstar*lambda3;
                beta3 = -beta1;
            }
        }

        // Update contributions
        hx = h[bottom]+alpha1-beta1/lambda1;
        hxx = h[top]-alpha3+beta3/lambda3;
        if(fabs(hx)<tol9) hx = 0.0;
        if(fabs(hxx)<tol9) hxx = 0.0;


        // First wave
        if(h[bottom]<tol9 && hx < 0.0){	// dry-wet wall
            DU1[top] += (lambda1*alpha1-beta1)*e1[0];
            DU2[top] += 0.0;
            DU3[top] += 0.0;
            DU4[top] += 0.0;
            DU5[top] += 0.0;

        } else if(h[top]<tol9 && hxx < 0.0){ // wet-dry wall
            DU1[bottom] += (lambda1*alpha1-beta1)*e1[0];
            DU2[bottom] += 0.0;
            DU3[bottom] += 0.0;
            DU4[bottom] += 0.0;
            DU5[bottom] += 0.0;

        } else if(lb1l < 0.0 && lb1r > 0.0){
            lambda1minus = lb1l*(lb1r-lambda1)/(lb1r-lb1l);
            DU1[bottom] += (lambda1minus*alpha1-beta1)*e1[0];
            DU2[bottom] += (lambda1minus*alpha1-beta1)*e1[1];
            DU3[bottom] += (lambda1minus*alpha1-beta1)*e1[2];
            DU4[bottom] += (lambda1minus*alpha1-beta1)*e1[3];
            DU5[bottom] += (lambda1minus*alpha1-beta1)*e1[4];
            lambda1plus = lb1r*(lambda1-lb1l)/(lb1r-lb1l);
            DU1[top] += (lambda1plus*alpha1)*e1[0];
            DU2[top] += (lambda1plus*alpha1)*e1[1];
            DU3[top] += (lambda1plus*alpha1)*e1[2];
            DU4[top] += (lambda1plus*alpha1)*e1[3];
            DU5[top] += (lambda1plus*alpha1)*e1[4];
            //ni flower que hacer con los DU4 y DU5

        } else if(lambda1 < 0.0){
            DU1[bottom] += (lambda1*alpha1-beta1)*e1[0];
            DU2[bottom] += (lambda1*alpha1-beta1)*e1[1];
            DU3[bottom] += (lambda1*alpha1-beta1)*e1[2];
            DU4[bottom] += (lambda1*alpha1-beta1)*e1[3];
            DU5[bottom] += (lambda1*alpha1-beta1)*e1[4];

        } else if(lambda1 >= 0.0){
            DU1[top] += (lambda1*alpha1-beta1)*e1[0];
            DU2[top] += (lambda1*alpha1-beta1)*e1[1];
            DU3[top] += (lambda1*alpha1-beta1)*e1[2];
            DU4[top] += (lambda1*alpha1-beta1)*e1[3];
            DU5[top] += (lambda1*alpha1-beta1)*e1[4];
        }

        // Second wave
        if(h[bottom]<tol9 && hx < 0.0){ 		// dry-wet wall
            DU1[top] += (lambda2*alpha2-beta2)*e2[0];
            DU2[top] += 0.0;
            DU3[top] += 0.0;
            DU4[top] += 0.0;
            DU5[top] += 0.0;

        } else if(h[top]<tol9 && hxx < 0.0){	// wet-dry wall
            DU1[bottom] += (lambda2*alpha2-beta2)*e2[0];
            DU2[bottom] += 0.0;
            DU3[bottom] += 0.0;
            DU4[bottom] += 0.0;
            DU5[bottom] += 0.0;

        } else if(lambda2 < 0.0){
            DU1[bottom] += (lambda2*alpha2-beta2)*e2[0];
            DU2[bottom] += (lambda2*alpha2-beta2)*e2[1];
            DU3[bottom] += (lambda2*alpha2-beta2)*e2[2];
            DU4[bottom] += (lambda2*alpha2-beta2)*e2[3];
            DU5[bottom] += (lambda2*alpha2-beta2)*e2[4];

        } else if(lambda2 >= 0.0){
            DU1[top] += (lambda2*alpha2-beta2)*e2[0];
            DU2[top] += (lambda2*alpha2-beta2)*e2[1];
            DU3[top] += (lambda2*alpha2-beta2)*e2[2];
            DU4[top] += (lambda2*alpha2-beta2)*e2[3];
            DU5[top] += (lambda2*alpha2-beta2)*e2[4];
        }

        // Third wave
        if(h[bottom]<tol9 && hx < 0.0){ 		// dry-wet wall
            DU1[top] += (lambda3*alpha3-beta3)*e3[0];
            DU2[top] += 0.0;
            DU3[top] += 0.0;
            DU4[top] += 0.0;
            DU5[top] += 0.0;

        } else if(h[top]<tol9 && hxx < 0.0){	// wet-dry wall
            DU1[bottom] += (lambda3*alpha3-beta3)*e3[0];
            DU2[bottom] += 0.0;
            DU3[bottom] += 0.0;
            DU4[bottom] += 0.0;
            DU5[bottom] += 0.0;

        } else if(lb3l < 0.0 && lb3r > 0.0){
            lambda3minus=lb3l*(lb3r-lambda3)/(lb3r-lb3l);
            DU1[bottom] += (lambda3minus*alpha3)*e3[0];
            DU2[bottom] += (lambda3minus*alpha3)*e3[1];
            DU3[bottom] += (lambda3minus*alpha3)*e3[2];
            DU4[bottom] += (lambda3minus*alpha3)*e3[3];
            DU5[bottom] += (lambda3minus*alpha3)*e3[4];
            lambda3plus=lb3r*(lambda3-lb3l)/(lb3r-lb3l);
            DU1[top] += (lambda3plus*alpha3-beta3)*e3[0];
            DU2[top] += (lambda3plus*alpha3-beta3)*e3[1];
            DU3[top] += (lambda3plus*alpha3-beta3)*e3[2];
            DU4[top] += (lambda3plus*alpha3-beta3)*e3[3];
            DU5[top] += (lambda3plus*alpha3-beta3)*e3[4];
            //ni flower que hacer con los DU4 y DU5

        } else if(lambda3 < 0.0){
            DU1[bottom] += (lambda3*alpha3-beta3)*e3[0];
            DU2[bottom] += (lambda3*alpha3-beta3)*e3[1];
            DU3[bottom] += (lambda3*alpha3-beta3)*e3[2];
            DU4[bottom] += (lambda3*alpha3-beta3)*e3[3];
            DU5[bottom] += (lambda3*alpha3-beta3)*e3[4];

        } else if(lambda3 >= 0.0){
            DU1[top] += (lambda3*alpha3-beta3)*e3[0];
            DU2[top] += (lambda3*alpha3-beta3)*e3[1];
            DU3[top] += (lambda3*alpha3-beta3)*e3[2];
            DU4[top] += (lambda3*alpha3-beta3)*e3[3];
            DU5[top] += (lambda3*alpha3-beta3)*e3[4];
        }

        // Fourth wave
        if(h[bottom]<tol9 && hx < 0.0){ 		// dry-wet wall
            DU1[top] += (lambda4*alpha4-beta4)*e4[0];
            DU2[top] += 0.0;
            DU3[top] += 0.0;
            DU4[top] += 0.0;
            DU5[top] += 0.0;

        } else if(h[top]<tol9 && hxx < 0.0){	// wet-dry wall
            DU1[bottom] += (lambda4*alpha4-beta4)*e4[0];
            DU2[bottom] += 0.0;
            DU3[bottom] += 0.0;
            DU4[bottom] += 0.0;
            DU5[bottom] += 0.0;

        } else if(lambda4 < 0.0){
            DU1[bottom] += (lambda4*alpha4-beta4)*e4[0];
            DU2[bottom] += (lambda4*alpha4-beta4)*e4[1];
            DU3[bottom] += (lambda4*alpha4-beta4)*e4[2];
            DU4[bottom] += (lambda4*alpha4-beta4)*e4[3];
            DU5[bottom] += (lambda4*alpha4-beta4)*e4[4];

        } else if(lambda4 >= 0.0){
            DU1[top] += (lambda4*alpha4-beta4)*e4[0];
            DU2[top] += (lambda4*alpha4-beta4)*e4[1];
            DU3[top] += (lambda4*alpha4-beta4)*e4[2];
            DU4[top] += (lambda4*alpha4-beta4)*e4[3];
            DU5[top] += (lambda4*alpha4-beta4)*e4[4];
        }

        // Fifth wave
        if(h[bottom]<tol9 && hx < 0.0){ 		// dry-wet wall
            DU1[top] += (lambda5*alpha5-beta5)*e5[0];
            DU2[top] += 0.0;
            DU3[top] += 0.0;
            DU4[top] += 0.0;
            DU5[top] += 0.0;

        } else if(h[top]<tol9 && hxx < 0.0){	// wet-dry wall
            DU1[bottom] += (lambda5*alpha5-beta5)*e5[0];
            DU2[bottom] += 0.0;
            DU3[bottom] += 0.0;
            DU4[bottom] += 0.0;
            DU5[bottom] += 0.0;

        } else if(lambda5 < 0.0){
            DU1[bottom] += (lambda5*alpha5-beta5)*e5[0];
            DU2[bottom] += (lambda5*alpha5-beta5)*e5[1];
            DU3[bottom] += (lambda5*alpha5-beta5)*e5[2];
            DU4[bottom] += (lambda5*alpha5-beta5)*e5[3];
            DU5[bottom] += (lambda5*alpha5-beta5)*e5[4];

        } else if(lambda5 >= 0.0){
            DU1[top] += (lambda5*alpha5-beta5)*e5[0];
            DU2[top] += (lambda5*alpha5-beta5)*e5[1];
            DU3[top] += (lambda5*alpha5-beta5)*e5[2];
            DU4[top] += (lambda5*alpha5-beta5)*e5[3];
            DU5[top] += (lambda5*alpha5-beta5)*e5[4];
        }

    } // End wet walls
}

//Estas solo hacen el bucle sobre todas las celdas
///OJO, que ahora también contienen la condición de frontera periódica
void h_compute_x_fluxes(int nX, int nY,
    double *h, double *qx, double *qy, double *ux, double *uy,
	double *zb, double *n, double *qBx, double *qBy, double *Bx, double *By,
	double *DU1, double *DU2, double *DU3, double *DU4, double *DU5,
    double dx, double QIN, double HIN){

   	// Process horizontal interfaces
   	if (QIN ==-2 && HIN < 0.0){
        int i_aux, i, j, right, left;
        for (j = 0; j < nY; j++) {     // Rows
            for (i = 0; i < nX; i++) { // Columns
                i_aux = i + 1;
                if(i==nX-1) {i_aux = 0;}
                left = IDX(i,j,nX);
                right = IDX(i_aux,j,nX);
                h_compute_x_step(nX, nY,
                    h, qx, qy, ux, uy,
                    zb, n, qBx, qBy, Bx, By,
                    DU1, DU2, DU3, DU4, DU5, dx,
                    left, right);
            }
        } // End of fluxes calculation - Wall loop
   	}

   	else{
        for (int j = 0; j < nY; j++) {     // Rows
            for (int i = 0; i < nX-1; i++) { // Columns
                int left = IDX(i,j,nX);
                int right = IDX(i+1,j,nX);
                h_compute_x_step(nX, nY,
                    h, qx, qy, ux, uy,
                    zb, n, qBx, qBy, Bx, By,
                    DU1, DU2, DU3, DU4, DU5, dx,
                    left, right);
            }
        }
   	}
   	#ifdef debug
        printf("x");
        print_DU(nX, nY,
            h, qx, qy, ux, uy,
            DU1, DU2, DU3, DU4, DU5,
            dx);
    #endif // debug
}

void h_compute_y_fluxes(int nX, int nY,
	double *h, double *qx, double *qy, double *ux, double *uy,
	double *zb, double *n, double *qBx, double *qBy, double *Bx, double *By,
	double *DU1, double *DU2, double *DU3, double *DU4, double *DU5,
	double dx, double QIN, double HIN){

    // Process vertical interfaces
    if (QIN <0.0 && HIN ==-2){
        int j_aux, j, i;
        for (j = 0; j < nY; j++) {     // Rows
            for (i = 0; i < nX; i++) { // Columns
                j_aux = j + 1;
                if(j==nY-1) {j_aux = 0;}
                int bottom =IDX(i,j,nX);
                int top = IDX(i,j_aux,nX);
                h_compute_y_step(nX, nY,
                    h, qx, qy, ux, uy,
                    zb, n, qBx, qBy, Bx, By,
                    DU1, DU2, DU3, DU4, DU5, dx,
                    top, bottom);
            }
        } // End of fluxes calculation - Wall loop
    }

    else{
        for (int j = 0; j < nY-1; j++) {     // Rows
            for (int i = 0; i < nX; i++) { // Columns
                int bottom =IDX(i,j,nX);
                int top = IDX(i,j+1,nX);
                h_compute_y_step(nX, nY,
                    h, qx, qy, ux, uy,
                    zb, n, qBx, qBy, Bx, By,
                    DU1, DU2, DU3, DU4, DU5, dx,
                    top, bottom);
            }
        }
    }
    #ifdef debug
        printf("DUy");
        print_DU(nX, nY,
            h, qx, qy, ux, uy,
            DU1, DU2, DU3, DU4, DU5,
            dx);
    #endif // debug
}





void h_check_depth_positivity(int nCells,
	double *h, double *DU1, double dx, double *dt){

    double aux1;

    for(int ic = 0; ic < nCells; ic++){

		aux1 = h[ic] - (*dt)*DU1[ic]/dx;


		if(fabs(aux1) < tol9) aux1 = 0.0;

		while(aux1 < 0.0) {
			(*dt) = 0.50 * (*dt);

			aux1 = h[ic] - (*dt)*DU1[ic]/dx;
			if(fabs(aux1) < tol9) aux1 = 0.0;
		}
	}

}



///Aquí hay alguna diferencia que no tengo clara, pues puede haber campo sin "agua" pero weno, por ahora asumimos que no
void h_update_cells_2D(int nCells,
	double *h, double *qx, double *qy, double *ux, double *uy,
    double *qBx, double *qBy, double *Bx, double *By,
	double *DU1, double *DU2, double *DU3, double *DU4, double *DU5,
	double dx, double dt){

        for(int ic = 0; ic < nCells; ic++){
		h[ic] = h[ic] - DU1[ic]*dt/dx;
		qx[ic] = qx[ic] - DU2[ic]*dt/dx;
		qy[ic] = qy[ic] - DU3[ic]*dt/dx;
		qBx[ic] = qBx[ic] - DU4[ic]*dt/dx;
        qBy[ic] = qBy[ic] - DU5[ic]*dt/dx;


		if(h[ic] >= tol9){
			// Minimum depth control
			if(h[ic] >= hmin){
				ux[ic] = qx[ic]/h[ic];
				uy[ic] = qy[ic]/h[ic];
				Bx[ic] = qBx[ic]/h[ic];
				By[ic] = qBy[ic]/h[ic];
			}else{
				qx[ic] = tol12;
				qy[ic] = tol12;
				ux[ic] = tol12;
				uy[ic] = tol12;
				qBx[ic] = tol12;
				qBy[ic] = tol12;
				Bx[ic] = tol12;
				By[ic] = tol12;
			}
		} else {
			h[ic] = tol12;
			qx[ic] = tol12;
            qy[ic] = tol12;
            ux[ic] = tol12;
            uy[ic] = tol12;
            qBx[ic] = tol12;
            qBy[ic] = tol12;
            Bx[ic] = tol12;
            By[ic] = tol12;
		}

		// Reset U fluxes differences
		DU1[ic]=0.0;
		DU2[ic]=0.0;
		DU3[ic]=0.0;
		DU4[ic]=0.0;
		DU5[ic]=0.0;

	}



}

void h_wet_dry_x(int nX, int nY,
	double *h, double *qx, double *ux,	double *zb){
	for (int j = 0; j < nY; j++) {     // Rows
		for (int i = 0; i < nX-1; i++) { // Columns
            int left = IDX(i,j,nX);
            int right = IDX(i+1,j,nX);
			if((h[right] < tol9) && (h[left] + zb[left] <zb[right])){
				qx[left]=0.0;
				ux[left]=0.0;
			}
			if((h[left] < tol9) && (h[right] + zb[right] <zb[left])){
				qx[right]=0.0;
				ux[right]=0.0;
			}
		}
	}
}

void h_wet_dry_y(int nX, int nY,
	double *h, double *qy, double *uy,	double *zb){

	// Process vertical interfaces
	for (int j = 0; j < nY-1; j++) {     // Rows
		for (int i = 0; i < nX; i++) { // Columns
			int bottom =IDX(i,j,nX);
			int top = IDX(i,j+1,nX);
			if((h[top] < tol9) && (h[bottom] + zb[bottom] <zb[top])){
				qy[bottom]=0.0;
				uy[bottom]=0.0;
			}
			if((h[bottom] < tol9) && (h[top] + zb[top] <zb[bottom])){
				qy[top]=0.0;
				uy[top]=0.0;
			}
		}
	}
}


///Ahora mismo, las condiciones de frontera funcionan solo en algunos casos super concretos, revisar
//Para poner salida libre, es tan sencillo como poner algo tipo -3 en todo
//Por razones que no explicare, las condiciones de Von-Karman (se activan con -2) están en las funciones
void h_set_west_boundary(int nX, int nY, double *h,
	double *qx, double *qy, double *ux, double *uy,
    double *qBx, double *qBy, double *Bx, double *By,
	double QIN, double HIN){
	for(int j = 0; j < nY; j++) { //
		int idx = IDX(0,j,nX);

		if(QIN > 0.0){
			qx[idx] = QIN;
			qy[idx] = 0.0;
		}
		if(HIN > 0.0) {
			h[idx] = HIN;
			qy[idx] = 0.0;
		}


		//Pared Rígida
		if(QIN ==-1 && HIN ==-1){
			qx[idx]=0.0;
		}


		if(h[idx]<0.0) h[idx]=0.0;

		if(h[idx]>=hmin){
			ux[idx] = qx[idx]/h[idx];
			uy[idx] = qy[idx]/h[idx];
			Bx[idx] = qBx[idx]/h[idx];
			By[idx] = qBy[idx]/h[idx];
		}else{
			ux[idx] = 0.0;
			uy[idx] = 0.0;
			Bx[idx] = 0.0;
			By[idx] = 0.0;
		}
	}
}

void h_set_east_boundary(int nX, int nY, double *h,
	double *qx, double *qy, double *ux, double *uy, double *zb,
    double *qBx, double *qBy, double *Bx, double *By,
	double HOUT, double ZSOUT){

    int idx;

	for(int j = 0; j < nY; j++) { //
		idx = IDX(nX-1,j,nX);

		if(HOUT > 0.0) {
			h[idx] = HOUT;
			qy[idx] = 0.0;
		}

		if(ZSOUT > 0.0) {
			h[idx] = ZSOUT-zb[idx];
			qy[idx] = 0.0;
		}

		if(HOUT ==-1 && ZSOUT ==-1){
			qx[idx]=0.0;
		}


		if(h[idx]<0.0) h[idx]=0.0;

		if(h[idx]>=hmin){
			ux[idx] = qx[idx]/h[idx];
			uy[idx] = qy[idx]/h[idx];
			Bx[idx] = qBx[idx]/h[idx];
			By[idx] = qBy[idx]/h[idx];
		}else{
			ux[idx] = 0.0;
			uy[idx] = 0.0;
			Bx[idx] = 0.0;
			By[idx] = 0.0;
		}
	}
}

void h_set_north_boundary(int nX, int nY, double *h,
    double *qBx, double *qBy, double *Bx, double *By,
	double *qx, double *qy, double *ux, double *uy,
	double HOUT, double ZSOUT){
	for(int i = 0; i < nX; i++) {
		int idx = IDX(i,nY-1,nX);


        if(HOUT == -1 && ZSOUT == -1){
            qy[idx] = 0.0;
        }

        if(h[idx]<0.0) h[idx]=0.0;

		if(h[idx]>=hmin){
			ux[idx] = qx[idx]/h[idx];
			uy[idx] = qy[idx]/h[idx];
			Bx[idx] = qBx[idx]/h[idx];
			By[idx] = qBy[idx]/h[idx];
		}else{
			ux[idx] = 0.0;
			uy[idx] = 0.0;
			Bx[idx] = 0.0;
			By[idx] = 0.0;
		}
    }

}

void h_set_south_boundary(int nX, int nY, double *h,
    double *qBx, double *qBy, double *Bx, double *By,
	double *qx, double *qy, double *ux, double *uy,
	double HOUT, double ZSOUT){

	for(int i = 0; i < nX; i++) {
		int idx = IDX(i,0,nX);


        if(HOUT == -1 && ZSOUT == -1){
            qy[idx] = 0.0;
        }


		if(h[idx]<0.0) h[idx]=0.0;

		if(h[idx]>=hmin){
			ux[idx] = qx[idx]/h[idx];
			uy[idx] = qy[idx]/h[idx];
			Bx[idx] = qBx[idx]/h[idx];
			By[idx] = qBy[idx]/h[idx];
		}else{
			ux[idx] = 0.0;
			uy[idx] = 0.0;
			Bx[idx] = 0.0;
			By[idx] = 0.0;
		}
	}

}



int write_vtk_cells(const char *filename, int nX, int nY, double *x, double *y,
	double *zb, double *h, double *ux, double *uy,
    double *qBx, double *qBy, double *Bx, double *By){

    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening output vtk file");
        return -1;
    }

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "VTK file containing raster data\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(file, "POINTS %d float\n", (nX + 1) * (nY + 1));

    // Write coordinates
    for (int j = 0; j <= nY; j++) {
        for (int i = 0; i <= nX; i++) {
            fprintf(file, "%.6f %.6f 0.0\n", x[IDX(i,j,nX+1)], y[IDX(i,j,nX+1)]);
        }
    }

    // Define the cells (connectivity)
    fprintf(file, "CELLS %d %d\n", nX*nY, 5*nX*nY);
    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            // For structured grid, points are ordered by their indices
            int p0 = IDX(i, j, nX+1);
            int p1 = IDX(i+1, j, nX+1);
            int p2 = IDX(i+1, j+1, nX+1);
            int p3 = IDX(i, j+1, nX+1);
            fprintf(file, "4 %d %d %d %d\n", p0, p1, p2, p3);
        }
    }

    // Define the cell types
    fprintf(file, "CELL_TYPES %d\n", nX*nY);
    for (int i = 0; i < nX*nY; i++) {
        fprintf(file, "9\n");  // 9 is the VTK code for quad cells
    }

    // Write elevation values as CELL_DATA
    fprintf(file, "CELL_DATA %d\n", nX*nY);
    fprintf(file, "SCALARS z double\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            fprintf(file, "%.6f\n", zb[IDX(i,j,nX)]);
        }
    }

    fprintf(file, "SCALARS h double\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            fprintf(file, "%.6f\n", h[IDX(i,j,nX)]);
        }
    }
	fprintf(file, "SCALARS h+z double\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            fprintf(file, "%.6f\n", h[IDX(i,j,nX)] + zb[IDX(i,j,nX)]);
        }
    }

    fprintf(file, "VECTORS v double\n");
    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            fprintf(file, "%.6f %.6f 0.0\n", ux[IDX(i,j,nX)], uy[IDX(i,j,nX)]);
        }
    }

    fprintf(file, "VECTORS B double\n");
    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            fprintf(file, "%.6f %.6f 0.0\n", Bx[IDX(i,j,nX)], By[IDX(i,j,nX)]);
        }
    }

    fclose(file);
    return 0;
}


void liberate_memory(double *x, double *y,
	double *zb, double *h, double *ux, double *uy,
    double *qBx, double *qBy, double *Bx, double *By,
    double *n, double *DU1, double *DU2, double *DU3,
    double *DU4, double *DU5, double *qx, double *qy){

	free(h);
	free(qx);
	free(qy);
	free(qBx);
	free(qBy);

	free(ux);
	free(uy);
	free(Bx);
	free(By);

	free(n);
	free(zb);

	free(DU1);
	free(DU2);
	free(DU3);
	free(DU4);
	free(DU5);
}


#ifdef debug
    void print_variables(int nX, int nY, double *h, double *qx, double *qy,
        double *ux, double *uy, double *qBx, double *qBy, double *Bx,
        double *By,	double *DU1, double *DU2, double *DU3, double *DU4,
        double *DU5, double dx, double dt){

        printf("\n h:\n");
        for (int j = 0; j < nX; j++){
            for(int i = 0; i < nY; i++){
                printf("%.14f\t", h[IDX(i,j,nX)]);
            }
            printf("\n");
        }
        printf("\n ux:\n");
        for (int j = 0; j < nX; j++){
            for(int i = 0; i < nY; i++){
                printf("%.14f\t", ux[IDX(i,j,nX)]);
            }
            printf("\n");
        }
        printf("\n uy:\n");
        for (int j = 0; j < nX; j++){
            for(int i = 0; i < nY; i++){
                printf("%.14f\t", uy[IDX(i,j,nX)]);
            }
            printf("\n");
        }
        printf("\n Bx:\n");
        for (int j = 0; j < nX; j++){
            for(int i = 0; i < nY; i++){
                printf("%.14f\t", Bx[IDX(i,j,nX)]);
            }
            printf("\n");
        }
        printf("\n By:\n");
        for (int j = 0; j < nX; j++){
            for(int i = 0; i < nY; i++){
                printf("%.14f\t", By[IDX(i,j,nX)]);
            }
            printf("\n");
        }
        printf("-----------------------------------------------------");
    }

    void print_DU(int nX, int nY, double *h, double *qx, double *qy,
    double *ux, double *uy, double *DU1, double *DU2, double *DU3, double *DU4, double *DU5, double dx){
    printf("\n DU1:\n");
    for (int j = 0; j < nX; j++){
        for(int i = 0; i < nY; i++){
           printf("%.14f\t", DU1[IDX(i,j,nX)]);
        }
       printf("\n");
    }
   printf("\n DU2:\n");
    for (int j = 0; j < nX; j++){
        for(int i = 0; i < nY; i++){
           printf("%.14f\t", DU2[IDX(i,j,nX)]);
        }
       printf("\n");
    }
   printf("\n DU3:\n");
    for (int j = 0; j < nX; j++){
        for(int i = 0; i < nY; i++){
           printf("%.14f\t", DU3[IDX(i,j,nX)]);
        }
       printf("\n");
    }
    printf("\n DU4:\n");
    for (int j = 0; j < nX; j++){
        for(int i = 0; i < nY; i++){
           printf("%.14f\t", DU4[IDX(i,j,nX)]);
        }
       printf("\n");
    }
    printf("\n DU5:\n");
    for (int j = 0; j < nX; j++){
        for(int i = 0; i < nY; i++){
           printf("%.14f\t", DU5[IDX(i,j,nX)]);
        }
       printf("\n");
    }
}

    void print_eigen(int row, int column, double *e1, double *e2, double *e3, double *e4, double *e5, double alpha1, double alpha2, double alpha3, double alpha4, double alpha5){
        printf("Pared: row=%d, column=%d\n", row, column);
        printf("e1:");
        for (int i=0; i<5; i++)
            printf("%.12f\t", e1[i]);
        printf("\n");
        printf("e2:");
        for (int i=0; i<5; i++)
            printf("%.12f\t", e2[i]);
        printf("\n");
        printf("e3:");
        for (int i=0; i<5; i++)
            printf("%.12f\t", e3[i]);
        printf("\n");
        printf("e4:");
        for (int i=0; i<5; i++)
            printf("%.12f\t", e4[i]);
        printf("\n");
        printf("e5:");
        for (int i=0; i<5; i++)
            printf("%.12f\t", e5[i]);
        printf("\n");

        printf("alpha1: %.12f\t", alpha1);
        printf("alpha2: %.12f\t", alpha2);
        printf("alpha3: %.12f\t", alpha3);
        printf("alpha4: %.12f\t", alpha4);
        printf("alpha5: %.12f\t", alpha5);
        printf("\n*********************************\n");

    }

#endif //debug
