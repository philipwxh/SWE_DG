#include "SWE_test_case.hpp"

void VortexSWESolution2d( MatrixXd x, MatrixXd y, double t, dfloat &g,
		                      MatrixXd &h, MatrixXd &u, MatrixXd &v, MatrixXd &b ){
  
  /*  This function evalute exact solution for 2D shallow water equation with vortex initial
      condition
      input:
        MatrixXd x:  array of x coordinates of the mesh points
        MatrixXd y:  array of y coordinates of the mesh points
        double t:    time for the SWE equation 
        MatrixXd &h: variable to hold the computed result for h 
        MatrixXd &u: variable to hold the computed result for u 
        MatrixXd &v: variable to hold the computed result for v
  */

  MatrixXd H_inf( x.rows(), x.cols() );
  H_inf.fill( 1.0 );
  MatrixXd u_inf( x.rows(), x.cols() );
  u_inf.fill( 1.0 );
  MatrixXd v_inf = x * 0;
  MatrixXd xc = x * 0;
  MatrixXd yc = y * 0;

  MatrixXd xt = x - xc - u_inf * t;
  MatrixXd yt = y - yc - v_inf * t;
  MatrixXd r_sq = xt.array().square() + yt.array().square();
  double beta = 5.0;
  h = H_inf.array() - beta * beta / ( 32.0 * PI * PI ) * ( -2.0 * ( r_sq.array() - 1.0 ) ).exp();
  u = u_inf.array() - beta   / ( 2.0 * PI ) * ( -( r_sq.array() - 1.0 ) ).exp() * yt.array();
  v = v_inf.array() + beta   / ( 2.0 * PI ) * ( -( r_sq.array() - 1.0 ) ).exp() * xt.array();
  if( g < 0 ){
    g = 2.0;
    b = x * 0;
    printf( "Setting up the initial condition to run VortexSWESolution2d initial condition\n" );
  }
}

void H2restSWESolution2d( MatrixXd x, MatrixXd y, double t, dfloat &g,
                          MatrixXd &h, MatrixXd &u, MatrixXd &v, MatrixXd &b ){
  
  /*  This function evalute exact solution for 2D shallow water equation with vortex initial
      condition
      input:
        MatrixXd x:  array of x coordinates of the mesh points
        MatrixXd y:  array of y coordinates of the mesh points
        double t:    time for the SWE equation 
        MatrixXd &h: variable to hold the computed result for h 
        MatrixXd &u: variable to hold the computed result for u 
        MatrixXd &v: variable to hold the computed result for v
  */

  u = x * 0;
  v = y * 0;
  h = u.array() + 2;

  if( g < 0 ){
    g = 9.8;
    b = x * 0;
    printf( "Setting up the initial condition to run H2restSWESolution2d initial condition\n" );
  }
}

void H2BTMSWESolution2d( MatrixXd x, MatrixXd y, double t, dfloat &g,
                         MatrixXd &h, MatrixXd &u, MatrixXd &v, MatrixXd &b ){
  
  /*  This function evalute exact solution for 2D shallow water equation with vortex initial
      condition
      input:
        MatrixXd x:  array of x coordinates of the mesh points
        MatrixXd y:  array of y coordinates of the mesh points
        double t:    time for the SWE equation 
        MatrixXd &h: variable to hold the computed result for h 
        MatrixXd &u: variable to hold the computed result for u 
        MatrixXd &v: variable to hold the computed result for v
  */

  MatrixXd H( x.rows(), x.cols() );
  H.fill( 2.0 );

  u = x * 0;
  v = y * 0;

  if( g < 0 ){
    g = 9.8;
    b = 0.1 * ( 2 * PI * x.array() ).sin() * ( 2 * PI * y.array() ).cos() + 0.5;
    printf( "Setting up the initial condition to run H2BTMSWESolution2d initial condition\n" );    
  }
  h = H - b;
}