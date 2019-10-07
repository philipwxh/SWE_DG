#ifndef _SWE_FILEIO_HEADER
#define _SWE_FILEIO_HEADER

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <unistd.h>

#include "fem.h"

void SWE_fileio( double Lx, double Ly, int N, int K1D, int K, int Nsteps, 
     double curve, double CFL, double FinalTime, double dt, double L2_err, double Linf_err, 
     clock_t time_setup, clock_t time_iteration, clock_t time_error, clock_t time_total );

void SWE_error_output( Mesh *mesh, MatrixXd &Q, MatrixXd &Vq, MatrixXd &btm,
					   double g, double FinalTime, double CFL, int N, int K1D,
     void ( *SWE_solution )( MatrixXd, MatrixXd, double, double &, MatrixXd &, MatrixXd &, MatrixXd &, MatrixXd &));

#endif