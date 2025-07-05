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
#include <math.h>
//#include <iostream>
//#include <string>
//#include <vector>
#include <time.h>       			// clock_t, clock(), CLOCKS_PER_SEC
#include <unistd.h>					// Para usleep()



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Mesh discretization
#define ZRlevel 0.0			  // Reference bed elevation



////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Simulation setup
#define Niter 100       // Number of iterations for screen output
#define DISPLAY 0

#define friction 0			// Switch off friction term
#define Fmodel 1			// Basal resistance model (1: Turbulent model | 2: Viscous model)
#define bed_slope_integral 0			// swtich between integral and differential bed slope formulation

#define rhow 1000.0   // Water density
#define mu 0.000			// Water dynamic viscosity
#define hmin 0.001    // Minimum depth for advective flux



//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Phisical-mathematical parameters
#define PI 3.14159265358979323846
#define g 1
#define tol6 1e-6
#define tol9 1e-9
#define tol12 1e-12



// Macros
#define MIN(x,y) (((x)<(y)) ? (x) : (y))			// Macro to find the minimum of two numbers
#define MAX(a,b) (((a)>(b)) ? (a) : (b))            // Macro to find the maximum of two numbers
#define IDX(i, j, ncols) ((j) * (ncols) + (i))

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#define debug
