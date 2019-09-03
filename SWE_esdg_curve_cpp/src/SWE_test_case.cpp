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

void SWE_sol_err( int N, int K1D, double g, double FinalTime, Mesh *mesh, MatrixXd &Q, MatrixXd &btm, 
  void ( *SWE_solution )( MatrixXd, MatrixXd, double, double &, MatrixXd &, MatrixXd &, MatrixXd &, MatrixXd &),
                   double &Linf_err, double &L2_err ){
  // update the 2D SWE term h, hu and hv onto CPU
  MatrixXd h  = Q.middleRows( 0, mesh->Np );
  MatrixXd hu = Q.middleRows( mesh->Np, mesh->Np );
  MatrixXd hv = Q.middleRows( 2 * mesh->Np, mesh->Np );
  MatrixXd hex_inf, uex_inf, vex_inf;
  SWE_solution( mesh->x, mesh->y, FinalTime, g, hex_inf, uex_inf, vex_inf, btm );

  MatrixXd huex_inf = hex_inf.array() * uex_inf.array();
  MatrixXd hvex_inf = hex_inf.array() * vex_inf.array();
  // evaluate the L2 error
  MatrixXd SWE_hError  = ( hex_inf  - h  ).array().abs();
  MatrixXd SWE_huError = ( huex_inf - hu ).array().abs();
  MatrixXd SWE_hvError = ( hvex_inf - hv ).array().abs();
  Linf_err = max( SWE_hError.maxCoeff(), max( SWE_huError.maxCoeff(), SWE_hvError.maxCoeff() ) );
  printf( "L_inf for SWE = %g\n", Linf_err );

  // finer quadrature for error eval
  Mesh *mesh_err = new Mesh;
  TriMesh2d( mesh_err, 2 * K1D, K1D ); 
  InitRefTri( mesh_err, N + 2 );

  MatrixXd wq_err = mesh_err->wq;

  MatrixXd Vqtmp = Vandermonde2D( N, mesh_err->rq, mesh_err->sq );
  MatrixXd Vq = Vandermonde2D( N, mesh->rq, mesh->sq );  
  MatrixXd Vq2 = mrdivide( Vqtmp, Vq );
  MatrixXd wJq = wq_err.asDiagonal() * ( Vq2 * ( mesh->Vq * mesh->J ) );  

  MatrixXd xq2 = Vq2 * ( mesh->Vq * mesh->x );
  MatrixXd yq2 = Vq2 * ( mesh->Vq * mesh->y );
  MatrixXd btmq = Vq2 * ( mesh->Vq * btm );
  MatrixXd hex, uex, vex;

  SWE_solution( xq2, yq2, FinalTime, g, hex, uex, vex, btmq );

  MatrixXd huex = hex.array() * uex.array();
  MatrixXd hvex = hex.array() * vex.array();

  MatrixXd werr = wJq.array()*(
                  ( hex  - Vq2 * ( mesh->Vq * h  ) ).array().square() +
                  ( huex - Vq2 * ( mesh->Vq * hu ) ).array().square() +
                  ( hvex - Vq2 * ( mesh->Vq * hv ) ).array().square() );
  L2_err = sqrt( werr.sum() );
  printf( "L2 error for SWE = %g\n", L2_err );
  delete( mesh_err );
}

void SWE_sol_err_sbp( int N, int K1D, double g, double FinalTime, Mesh *mesh, MatrixXd &Q, MatrixXd &btm, 
  void ( *SWE_solution )( MatrixXd, MatrixXd, double, double &, MatrixXd &, MatrixXd &, MatrixXd &, MatrixXd &),
                   double &Linf_err, double &L2_err ){
  // update the 2D SWE term h, hu and hv onto CPU
  MatrixXd h  = Q.middleRows( 0, mesh->Nq );
  // cout << "h block: row: " << h.rows() << " . col: "<< h.cols() << ". " << endl << h << endl << endl;
  MatrixXd hu = Q.middleRows( mesh->Nq, mesh->Nq );
  MatrixXd hv = Q.middleRows( 2 * mesh->Nq, mesh->Nq );
  MatrixXd hex_inf, uex_inf, vex_inf;
  MatrixXd xq_sbp = mesh->Vq * mesh->x;
  MatrixXd yq_sbp = mesh->Vq * mesh->y;

  SWE_solution( xq_sbp, yq_sbp, FinalTime, g, hex_inf, uex_inf, vex_inf, btm );

  MatrixXd huex_inf = hex_inf.array() * uex_inf.array();
  MatrixXd hvex_inf = hex_inf.array() * vex_inf.array();
  // evaluate the L2 error
  MatrixXd SWE_hError  = ( hex_inf  - h  ).array().abs();
  MatrixXd SWE_huError = ( huex_inf - hu ).array().abs();
  MatrixXd SWE_hvError = ( hvex_inf - hv ).array().abs();
  Linf_err = max( SWE_hError.maxCoeff(), max( SWE_huError.maxCoeff(), SWE_hvError.maxCoeff() ) );
  printf( "L_inf for SWE = %g\n", Linf_err );

  // finer quadrature for error eval
  Mesh *mesh_err = new Mesh;
  TriMesh2d( mesh_err, 2 * K1D, K1D ); 
  InitRefTri( mesh_err, N + 2 );

  MatrixXd wq_err = mesh_err->wq;

  MatrixXd Vqtmp = Vandermonde2D( N, mesh_err->rq, mesh_err->sq );
  MatrixXd Vq = Vandermonde2D( N, mesh->rq, mesh->sq );  
  MatrixXd Vq2 = mrdivide( Vqtmp, Vq );
  MatrixXd wJq = wq_err.asDiagonal() * ( Vq2 * ( mesh->Vq * mesh->J ) );  
  MatrixXd xq2 = Vq2 * xq_sbp;
  MatrixXd yq2 = Vq2 * yq_sbp;
  MatrixXd btmq = Vq2 * btm;

  MatrixXd hex, uex, vex;

  SWE_solution( xq2, yq2, FinalTime, g, hex, uex, vex, btmq );
  MatrixXd huex = hex.array() * uex.array();
  MatrixXd hvex = hex.array() * vex.array();

  MatrixXd werr = wJq.array()*(
                  ( hex  - Vq2 * h ).array().square() +
                  ( huex - Vq2 * hu ).array().square() +
                  ( hvex - Vq2 * hv ).array().square() );
  L2_err = sqrt( werr.sum() );
  printf( "L2 error for SWE = %g\n", L2_err );
  delete( mesh_err );
}