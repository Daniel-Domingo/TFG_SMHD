//   Shallow Water Equations 2D Simulator (SWE2D)
//  ===========================================================================================
//   Square grid
//  ===========================================================================================
//   Version 1.0 - Mar 2025
//  ===========================================================================================
//   Computational Hydraulics Group - University of Zaragoza
//  ===========================================================================================


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <cstring>
#include "../define.h"
void directory_founder(char *tempfile, const char *header,
    const char *case_name_slashed, const char *file_name, const char *extension);

int read_raster_header(const char *filename, int *nX, int *nY,
	double *Xmin, double *Ymin, double *dx);

int read_raster(const char *filename, int nX, int nY,
	double *data);

int read_simulation_setup(char *tempfile, double *simTime,
	double *Toutput, double *CFL, double *n1);

int read_boundary_configuration(char *tempfile,
	double *QIN, double *HIN,
	double *HOUT, double *ZSOUT);

void h_initialize_variables(int nCells,
	double *h, double *qx, double *qy, double *ux, double *uy,
	double *n, double *qBx, double *qBy, double *Bx, double *By,
	double *DU1, double *DU2, double *DU3, double *DU4, double *DU5);

int h_set_inlet_initial_conditions(int nCells,
	double *h, double *u, double *M,
	double QIN, double HIN);

int h_set_outlet_initial_conditions(int nCells,
	double *h, double *u, double *M, double *zb,
	double HOUT, double ZSOUT);

int h_compute_initial_flow_variables(int nCells,
	double *h, double *qx, double *qy, double *ux, double *uy,
	double *qBx, double *qBy, double *Bx, double *By,
	double *n, double n1);


void h_compute_water_mass(int nCells,
	double *h, double area, double *massWt);

void h_compute_divergence_condition(int nX, int nY, double *qBx, double *qBy,
    double *divergence_loss, double dx);

void h_compute_flow_time_step_2D(int nX, int nY,
    double *h, double *ux, double *uy, double dx,
	double *qBx, double *qBy, double *Bx, double *By,
	double *dtSW);

void h_compute_x_fluxes(int nX, int nY,
	double *h, double *qx, double *qy, double *ux, double *uy,
	double *zb, double *n, double *qBx, double *qBy, double *Bx, double *By,
	double *DU1, double *DU2, double *DU3, double *DU4, double *DU5,
	double dx, double QIN, double HIN);

void h_compute_y_fluxes(int nX, int nY,
	double *h, double *qx, double *qy, double *ux, double *uy,
	double *zb, double *n, double *qBx, double *qBy, double *Bx, double *By,
	double *DU1, double *DU2, double *DU3, double *DU4, double *DU5,
	double dx, double QIN, double HIN);

void h_check_depth_positivity(int nCells,
	double *h, double *DU1, double dx, double *dt);

void h_update_cells_2D(int nCells,
	double *h, double *qx, double *qy, double *ux, double *uy,
    double *qBx, double *qBy, double *Bx, double *By,
	double *DU1, double *DU2, double *DU3, double *DU4, double *DU5,
	double dx, double dt);

void h_wet_dry_x(int nX, int nY,
	double *h, double *qx, double *ux,	double *zb);

void h_wet_dry_y(int nX, int nY,
	double *h, double *qy, double *uy,	double *zb);

void h_set_west_boundary(int nX, int nY, double *h,
	double *qx, double *qy, double *ux, double *uy,
    double *qBx, double *qBy, double *Bx, double *By,
	double QIN, double HIN);

void h_set_east_boundary(int nX, int nY, double *h,
	double *qx, double *qy, double *ux, double *uy, double *zb,
    double *qBx, double *qBy, double *Bx, double *By,
	double HOUT, double ZSOUT);

void h_set_north_boundary(int nX, int nY, double *h,
    double *qBx, double *qBy, double *Bx, double *By,
	double *qx, double *qy, double *ux, double *uy,
	double HOUT, double ZSOUT);

void h_set_south_boundary(int nX, int nY, double *h,
    double *qBx, double *qBy, double *Bx, double *By,
	double *qx, double *qy, double *ux, double *uy,
	double HOUT, double ZSOUT);

int write_vtk_cells(const char *filename, int nX, int nY, double *x, double *y,
	double *zb, double *h, double *ux, double *uy,
    double *qBx, double *qBy, double *Bx, double *By);

void liberate_memory(double *x, double *y,
	double *zb, double *h, double *ux, double *uy,
    double *qBx, double *qBy, double *Bx, double *By,
    double *n, double *DU1, double *DU2, double *DU3,
    double *DU4, double *DU5, double *qx, double *qy);

#ifdef debug
    void print_variables(int nX, int nY, double *h, double *qx, double *qy,
        double *ux, double *uy, double *qBx, double *qBy, double *Bx,
        double *By,	double *DU1, double *DU2, double *DU3, double *DU4,
        double *DU5, double dx, double dt);

    void print_DU(int nX, int nY, double *h, double *qx, double *qy,
        double *ux, double *uy, double *DU1, double *DU2, double *DU3, double *DU4, double *DU5, double dx);

    void print_eigen(int row, int column, double *e1, double *e2, double *e3, double *e4, double *e5,
        double alpha1, double alpha2, double alpha3, double alpha4, double alpha5);
#endif //debug
