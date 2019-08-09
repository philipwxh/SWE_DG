#ifndef _SWE_TEST_CASE_HEADER
#define _SWE_TEST_CASE_HEADER


#include "fem.h"

void VortexSWESolution2d( MatrixXd x, MatrixXd y, double t, dfloat &g,
		                      MatrixXd &h, MatrixXd &u, MatrixXd &v, MatrixXd &b );

void H2restSWESolution2d( MatrixXd x, MatrixXd y, double t, dfloat &g,
                          MatrixXd &h, MatrixXd &u, MatrixXd &v, MatrixXd &b );

void H2BTMSWESolution2d( MatrixXd x, MatrixXd y, double t, dfloat &g,
                         MatrixXd &h, MatrixXd &u, MatrixXd &v, MatrixXd &b );

#endif