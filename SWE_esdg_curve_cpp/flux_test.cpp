#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <unistd.h>

#include "fem.h"
#include "SWE_test_case.hpp"
#include "SWE_fileio.hpp"
// SWESOLUTION: name of the initial condition to be tested on
#define tau 1.0
#define DEBUG 0
#define TOL_ZERO 1.e-14

int main( int argc, char **argv ){
  int N = 6;
  int K = 8192;
  int V1 = 1; int V2 = 1;
  // if more arguments are provided, set N and K1D to user defined values.
  if ( argc > 2 ){
    N = atoi( argv[1] );
    K = atoi( argv[2] );
    // print N and K1D for user confirmation
    printf( "setting N = %d, K = %d", N, K ) ;
  }
  if ( argc > 4 ){
    V1 = atoi( argv[3] );
    V2 = atoi( argv[4] );
    // print N and K1D for user confirmation
    printf( "setting V1 = %d, V2 = %d", V1, V2 ) ;
  }

  int Nfields = 3;
  int Nsteps = 100;

  MatrixXd Q = MatrixXd::Zero( N, N ).array() + 1.0 / N;

  // MatrixXd F( N * N, k );
  // MatrixXd f( N, Nfields );
  MatrixXd U( N * Nfields, K );
  MatrixXd hu = MatrixXd::Zero( N, K );
  MatrixXd hv = MatrixXd::Zero( N, K );
  MatrixXd h = hu.array() + 1;
  U << h, hu, hv;


  // cout << "Q: row: " << Q.rows() << " . col: "<< Q.cols() << ". " << endl << Q << endl << endl;
  // cout << "U: row: " << U.rows() << " . col: "<< U.cols() << ". " << endl << U << endl << endl;

  double g = 2.0;

  // create an App object to setup OCCA
  App *app = new App;
  
  // app->device.setup("mode: 'Serial'");
  // select the device we have
  app->device.setup( "mode: 'CUDA', device_id: 0" );
  cout << app->device.mode();

  // define macros for OCCA kernels
  // p_Np: number of dims
  app->props[ "defines/p_Nq" ] = N; 
  // p_Nq: number of quadrature points
  app->props[ "defines/p_K" ] = K;
  app->props[ "defines/p_Nfields" ] = Nfields;
  
  app->props[ "defines/p_KblkV1" ] = V1; 
  app->props[ "defines/p_KblkV2" ] = V2; 
  app->props[ "defines/p_g" ]      = ( double ) g;
  app->props[ "defines/p_g_half" ] = ( double ) 0.5 * g;
  // build occa kernels  
  string path = "okl/flux_test.okl";
  
  //testing
  occa::kernel flux_test1, flux_test2;
  flux_test1 = app->device.buildKernel( path.c_str(), "flux_test1", app->props );
  flux_test2 = app->device.buildKernel( path.c_str(), "flux_test2", app->props );

  
  occa::memory o_Q, o_U, o_rhs; 
  setOccaArray( app, Q, o_Q );
  setOccaArray( app, U, o_U );
  setOccaArray( app, MatrixXd::Zero( N * Nfields, K ), o_rhs ); 
  MatrixXd rhs = MatrixXd::Zero( N * Nfields, K );
  
  clock_t time_flux1 = clock();
  for ( int i = 0; i < Nsteps ; ++i ){
    
    flux_test1( K, o_Q, o_U, o_rhs );
 
  }
  app->device.finish();
  time_flux1 = clock() - time_flux1;

  getOccaArray( app, o_rhs, rhs );
  // cout << "rhs: row: " << rhs.rows() << " . col: "<< rhs.cols() << ". " << endl << rhs << endl << endl;

  clock_t time_flux2 = clock();
  for ( int i = 0; i < Nsteps ; ++i ){
    
    flux_test2( K, o_Q, o_U, o_rhs );

  }
  app->device.finish();
  time_flux2 = clock() - time_flux2;

  getOccaArray( app, o_rhs, rhs );
  // cout << "rhs: row: " << rhs.rows() << " . col: "<< rhs.cols() << ". " << endl << rhs << endl << endl;

  flux_test_fileio( N, K, time_flux1, time_flux2 );

  cout << "Time for the flux test 1 is " << (float) time_flux1 / CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Time for the flux test 2 is " << (float) time_flux2 / CLOCKS_PER_SEC << " seconds" << endl;
  return 0;
}
