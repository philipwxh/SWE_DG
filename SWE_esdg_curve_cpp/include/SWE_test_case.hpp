
#ifndef _SWE_TEST_CASE_HEADER
#define _SWE_TEST_CASE_HEADER

#include <iostream>
#include "fem.h"

void VortexSWESolution2d( MatrixXd x, MatrixXd y, double t, dfloat &g,
		                      MatrixXd &h, MatrixXd &u, MatrixXd &v, MatrixXd &b );

void H2restSWESolution2d( MatrixXd x, MatrixXd y, double t, dfloat &g,
                          MatrixXd &h, MatrixXd &u, MatrixXd &v, MatrixXd &b );

void H2BTMSWESolution2d( MatrixXd x, MatrixXd y, double t, dfloat &g,
                         MatrixXd &h, MatrixXd &u, MatrixXd &v, MatrixXd &b );

void SWE_sol_err( int N, int K1D, double g, double FinalTime, Mesh *mesh, MatrixXd &Q, MatrixXd &btm, 
  void ( *SWE_solution )( MatrixXd, MatrixXd, double, double &, MatrixXd &, MatrixXd &, MatrixXd &, MatrixXd &),
                   double &Linf_err, double &L2_err );

void SWE_sol_err_sbp( int N, int K1D, double g, double FinalTime, Mesh *mesh, MatrixXd &Q, MatrixXd &btm, 
  void ( *SWE_solution )( MatrixXd, MatrixXd, double, double &, MatrixXd &, MatrixXd &, MatrixXd &, MatrixXd &),
                   double &Linf_err, double &L2_err );
#endif