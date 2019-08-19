#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <unistd.h>

#include "fem.h"
#include "SWE_test_case.hpp"
#include "SWE_fileio.hpp"
// SWESOLUTION: name of the initial condition to be tested on
#define SWESOLUTION VortexSWESolution2d
#define tau 1.0
#define TOL_ZERO 1.e-14
#define CAL_RES 0
#define PRINT_ERR 0

int main( int argc, char **argv ){
  clock_t time_setup = clock();
  // N: default the degree of the polynomial
  int N = 3;
  // K1D: default number of elements in one direction.
  //      in 2D, there will be 16 by 16 quadrilateral or 16 by 16*2 triangles
  int K1D = 8;

  //occa::printModeInfo();   return 0;
  // if more arguments are provided, set N and K1D to user defined values.
  if ( argc > 2 ){
    N = atoi( argv[1] );
    K1D = atoi( argv[2] );
    // print N and K1D for user confirmation
    printf( "setting N = %d, K1D = %d\n", N, K1D) ;
  }
  // FinalTime: final time of the simulation
  double FinalTime = 0.5;
  // CFL: Courant–Friedrichs–Lewy (CFL) condition for time step
  double CFL = .125;
  // a: curved warping, a = 0 means not curved
  double curve = 0; // curved warping
  
  // if more arguments are provided, set CFL, FinalTime and a to user defined values.
  if ( argc > 4 ){
    CFL = atof( argv[3] );
    FinalTime = atof( argv[4] );
    curve = atof( argv[5] );
    // print CFL, FinalTime and a for user confirmation
    printf( "setting CFL = %f, T = %f, curved warping a = %f\n", CFL, FinalTime, curve );
  }
  
  // create a new mesh pointer
  Mesh *mesh = new Mesh;
 
  // create a triangular mesh
  TriMesh2d( mesh, 2 * K1D, K1D ); // make Cartesian mesh
  // [0,20] x [-5,5] for vortex
  // Lx: length in x direction
  // Ly: length in y direction
  double Lx = 10;
  double Ly = 5;
  // VX: x coordination of elements vertices
  // VY: y coordination of elements vertices
  mesh->VX = ( mesh->VX ) * Lx;
  mesh->VY = ( mesh->VY ) * Ly;
  /*
  cout << "VX = " << mesh->VX << endl;
  cout << "VY = " << mesh->VY << endl;
  */

  // ============ physics independent stuff ===========

  printf( "N = %d, K = %d\n", N, mesh->K );
  
  // initialize reference triangle element
  InitRefTri( mesh, N );

  // make physical nodes + geofacs
  MapTriNodes( mesh ); // low order mapping
  // initialize geometric factors
  GeometricFactors2d( mesh ); 
  // initialize normals in 2D
  Normals2d( mesh );

  // makes EToE, EToF
  ConnectElems( mesh );

  // data
  // Nfields: for 2D SWE, we have for variables to solve
  int Nfields = 3;
  //  
  mesh->Nfields = Nfields;
  // K: number of elements in the entire domain.
  int K = mesh->K;
  // Np: number of volume points in one elements
  int Np = mesh->Np;
  // Nq: number of quadrature points in one elements
  int Nq = mesh->Nq;    
  // NfpNfaces: number of face points multiplied by number of faces in one elements
  int NfqNfaces = mesh->Nfq * mesh->Nfaces;
  
  // build node maps
  // get the x and y coordinate of the face nodes
  // xf: x coordinate of the face nodes
  // yf: y coordinate of the face nodes
  // zf: filled with 0 so the BuildFaceNodeMaps know it's 2D
  MatrixXd xf = ( mesh->Vf ) * ( mesh->x );
  MatrixXd yf = ( mesh->Vf ) * ( mesh->y );
  MatrixXd zf( xf.rows(), xf.cols() ); 
  zf.fill(0.0);
  MatrixXi mapPq;
  // Build mapPq to store the connection between face nodes
  BuildFaceNodeMaps( mesh, xf, yf, zf, mapPq );

  // builds periodicity into the node maps
  double DX = mesh->VX.maxCoeff() - mesh->VX.minCoeff();
  double DY = mesh->VY.maxCoeff() - mesh->VY.minCoeff();
  // create the periodic boundary condition for 2D map
  // and store the face nodes connections for the boundary nodes
  MakeNodeMapsPeriodic2d( mesh, xf, yf, DX, DY, mapPq );

  // correct mapPq for Nfields > 1  
  MatrixXi mapPqNfields( NfqNfaces, K );
  // loop over all face nodes in all elements
  
  for( int i = 0; i < NfqNfaces * K; ++i ){
    int idP = mapPq( i );
    // e: element of the neighbor node
    int e = idP / NfqNfaces;
    // fidP: index of the neighbor node in its element
    int fidP = idP % NfqNfaces;
    // store the first index of the variable 
    mapPqNfields( i ) = fidP + e * NfqNfaces * ( Nfields + 1 );
  }  
  mesh->mapPq = mapPqNfields;

  // =================== set initial condition

  // extract the mesh points
  // x: x coordinates of the mesh points
  // y: y coordinates of the mesh points
  MatrixXd x = mesh->x;
  MatrixXd y = mesh->y;  
 
  // smooth isentropic vortex
  MatrixXd h, u, v, btm;  
  double time = 0.0;
  dfloat g = -1.0;
  // compute the initial condition  
  SWESOLUTION( x, y, time, g, h, u, v, btm );
  // hu: h * u
  MatrixXd hu = h.array() * u.array();
  // hv: h * v
  MatrixXd hv = h.array() * v.array();
  // stack all variable in to one matrix Q
  MatrixXd Q( ( Nfields + 1 ) * Np, K );
  Q << h, hu, hv, btm;

  // ========================== set up OCCA application

  // create an App object to setup OCCA
  App *app = new App;
  //  app->device.setup("mode: 'Serial'");
  // select the device we have
  app->device.setup( "mode: 'CUDA', device_id: 0" );

  // define macros for OCCA kernels
  app->props[ "defines/p_Np" ] = mesh->Np; // number of dims
  app->props[ "defines/p_Nq" ] = Nq;
  app->props[ "defines/p_NfqNfaces" ] = NfqNfaces;  
  app->props[ "defines/p_NqT" ] = Nq + NfqNfaces; // total quadrature point 
  app->props[ "defines/p_T" ] = max( mesh->Nfq * mesh->Nfaces, mesh->Np );
  
  // number of geometric terms
  // Nvgeo: rxJ, ryJ, sxJ, syJ, J
  int Nvgeo = 5; 
  // Nsgeo: nxJ,nyJ, sJ
  int Nfgeo = 3; 
  app->props[ "defines/p_Nvgeo" ] = Nvgeo;
  app->props[ "defines/p_Nfgeo" ] = Nfgeo;

  int KblkP = 8;
  int KblkV = 8;
  int KblkS = 10;
  int KblkU = 10;
  // number of elements in one group for project kernel
  app->props[ "defines/p_KblkP" ] = KblkP;
  // number of elements in one group for volume kernel
  app->props[ "defines/p_KblkV" ] = KblkV;
  // number of elements in one group for surface kernel
  app->props[ "defines/p_KblkS" ] = KblkS;
  // number of elements in one group for update kernel
  app->props[ "defines/p_KblkU" ] = KblkU;


  // switch dfloat type (double/float) in types.h
  // switch between single precision and double precision
  if ( sizeof( dfloat ) == 4 ){
    app->props[ "defines/USE_DOUBLE" ] = 0;
  }else{
    app->props[ "defines/USE_DOUBLE" ] = 1;
  }
  // 1.f is single precision floating point
  // 1.0 is double precision floating point
  app->props[ "defines/p_Nfields" ] = Nfields;
  if (sizeof( dfloat ) == 4 ){
    app->props[ "defines/p_tau" ]    = ( float ) tau;
    app->props[ "defines/p_g" ]      = ( float ) g;
    app->props[ "defines/p_g_half" ] = ( float ) 0.5 * g;
  }else{
    app->props[ "defines/p_tau" ]    = ( double ) tau;
    app->props[ "defines/p_g" ]      = ( double ) g;
    app->props[ "defines/p_g_half" ] = ( double ) 0.5 * g;
  }

  // interpolate to quad pts and store  
  //  // remove eventually...  
  //  setupOccaMesh2d(mesh,app); // store and pack mesh geofacs in OCCA mem

    // build occa kernels  
  string path = "okl/SWE2DTri_group.okl";


  //testing
  occa::kernel volume, surface, update, project;
  volume = app->device.buildKernel( path.c_str(), "volume", app->props );
  surface = app->device.buildKernel( path.c_str(), "surface", app->props );
  update = app->device.buildKernel( path.c_str(), "update", app->props );
  project = app->device.buildKernel( path.c_str(), "project", app->props );
  
  // ================= set node maps + boundary condition
  
  
  //bcFlag: Boundary condition flag. 0 for periodic boundary condition
  MatrixXi bcFlag( xf.rows(), xf.cols() );  
  bcFlag.fill( 0 ); // no BCs

  occa::memory o_bcFlag, o_mapPq;
  setOccaIntArray( app, bcFlag, o_bcFlag );
  setOccaIntArray( app, mesh->mapPq, o_mapPq );  

  // ================= construct discretization matrices
  // Vq
  MatrixXd Vq = mesh->Vq;
  MatrixXd Vf = mesh->Vf;
  MatrixXd M = Vq.transpose() * mesh->wq.asDiagonal() * Vq;
  MatrixXd VqTW = Vq.transpose() * mesh->wq.asDiagonal();
  MatrixXd VfTWf = Vf.transpose() * mesh->wf.asDiagonal();
  MatrixXd Pq = mldivide( M, VqTW );
  MatrixXd Lf = mldivide( M, VfTWf );

  // hybridized SBP ops
  int NqT = Nq + NfqNfaces;
  MatrixXd Ef = Vf * Pq;
  MatrixXd Wf = mesh->wf.asDiagonal();
  MatrixXd Br = Wf * mesh->nrJ.asDiagonal();
  MatrixXd Bs = Wf * mesh->nsJ.asDiagonal();
  MatrixXd Qr = Pq.transpose() * ( M * mesh->Dr ) * Pq;
  MatrixXd Qs = Pq.transpose() * ( M * mesh->Ds ) * Pq;  
  MatrixXd Zf = MatrixXd::Zero( NfqNfaces, NfqNfaces );

  // skew form of the operators
  MatrixXd QNr( NqT, NqT );
  //QNr << Qr - Qr.transpose(), Ef.transpose() * Br, -Br * Ef, Zf;
  QNr << Qr-.5 * Ef.transpose() * Br * Ef, .5 * Ef.transpose() * Br, -.5 * Br * Ef, .5 * Br;
  MatrixXd QNs( NqT, NqT );
  //QNs << Qs - Qs.transpose(), Ef.transpose() * Bs, -Bs * Ef, Zf;
  QNs << Qs-.5 * Ef.transpose() * Bs * Ef, .5 * Ef.transpose() * Bs, -.5 * Bs * Ef, .5 * Bs;  

  MatrixXd VN( NqT, Np );
  VN << Vq,Vf;
  MatrixXd VNT = VN.transpose();
  MatrixXd VNP = VN * Pq;
  MatrixXd PN = mldivide( M, VNT ); 

  // BTMq:  Bottom gemetry interpolated to quadrature points
  MatrixXd BTMq = VN * btm;
  // BTMx:  Stores interoplated bottom geometry derivatives in x direction
  //        To get the BTMx for the e th element:
  //        BTMx.middleCols( e, 1 )
  MatrixXd gBTMx( NqT, K  );
  // BTMy:  Stores interoplated bottom geometry derivatives in y direction
  //        To get the BTMx for the e th element:
  //        BTMy.middleCols( e, 1 )
  MatrixXd gBTMy( NqT, K  );

  // Building elementwise Mass matrices inverses and BTMx and BTMy
  for( int e = 0; e < K; ++e ){
      // get information for element e
      // MatrixXd JJi = JJ.col( e ).asDiagonal();
      // MatrixXd WqJJi = mesh->wq.asDiagonal() * JJi;
      double rxJe = mesh->rxJ( 0 , e );
      double ryJe = mesh->ryJ( 0 , e );
      double sxJe = mesh->sxJ( 0 , e );
      double syJe = mesh->syJ( 0 , e );

      // build mass matrix, QNx and QNy for element e
      MatrixXd QNx_e = rxJe * QNr + sxJe * QNs;
      MatrixXd QNy_e = ryJe * QNr + syJe * QNs;
      
      gBTMx.middleCols( e, 1 )  = g * QNx_e * BTMq.col( e );
      gBTMy.middleCols( e, 1 )  = g * QNy_e * BTMq.col( e );
  }

  QNr = .5 * ( QNr - QNr.transpose() );
  QNs = .5 * ( QNs - QNs.transpose() );


  occa::memory o_Q, o_Qv, o_Qf; // solution   
  occa::memory o_vgeo, o_fgeo; // geofacs

  MatrixXd vgeo( Nvgeo*Nq, K );  
  vgeo << Vq * ( mesh->rxJ ),Vq * ( mesh->ryJ ),
    Vq * ( mesh->sxJ ), Vq * ( mesh->syJ ), Vq * mesh->J;

  setOccaArray( app,vgeo, o_vgeo );

  MatrixXd fgeo( Nfgeo * NfqNfaces, K );  
  fgeo << mesh->nxJ, mesh->nyJ, mesh->sJ; // already at quad pts
  setOccaArray( app, fgeo, o_fgeo );

  occa::memory o_rhs, o_res;  // timestep stuff
  setOccaArray( app, MatrixXd::Zero( Np * Nfields, K ), o_rhs );
  setOccaArray( app, MatrixXd::Zero( Np * Nfields, K ), o_res );  

  occa::memory o_VNP, o_QNr, o_QNs, o_PN, o_Lf, o_Vq; // operators
  occa::memory o_gBTMx, o_gBTMy; // precomputed data

  // set solution arrays
  setOccaArray( app, Q, o_Q );

  MatrixXd Qv( ( Nfields + 1 ) * Nq, K );
  Qv << Vq * h, Vq * hu, Vq * hv, Vq * btm;
  setOccaArray( app, Qv, o_Qv );
  // cout << "Qv block: row: " <<Qv.rows() << " . col: "<<Qv.cols() << ". " << endl << Qv << endl << endl;
  MatrixXd Qf( ( Nfields + 1 ) * NfqNfaces, K );
  Qf << Vf * h, Vf * hu, Vf * hv, Vf * btm;
  setOccaArray( app, Qf, o_Qf );
  // cout << "Qf block: row: " <<Qf.rows() << " . col: "<<Qf.cols() << ". " << endl << Qf << endl << endl;

  //setOccaArray( app, MatrixXd::Zero( Nfields * NfqNfaces, K ), o_Qf );    

  // set operators
//  setOccaArray( app, VN, o_VN );
  setOccaArray( app, Vq, o_Vq );
  setOccaArray( app, Lf, o_Lf );
  setOccaArray( app, QNr, o_QNr );
  setOccaArray( app, QNs, o_QNs );
  setOccaArray( app, PN, o_PN );
  setOccaArray( app, VNP, o_VNP );
  setOccaArray( app, gBTMx, o_gBTMx );
  setOccaArray( app, gBTMy, o_gBTMy );

  printf( "All matrices are set up\n" );
  
  // ============== run RK solver ==================

  //double h = 2.0 / (double) K1D; //min(DX,DY) / (double) K1D;
  double h_mesh = mesh->J.maxCoeff() / mesh->sJ.maxCoeff(); // J = O(h^d), Jf = O(h^{d-1}) in d dims
  h_mesh = 2.0 / (double) K1D;
  double CN = ( double ) ( ( N + 1 ) * ( N + 2 ) ) / 2.0; // trace constant for GQ hexes
  // dt: time step size
  double dt = CFL * h_mesh / CN;
  
  // Nsteps: number of time steps need to reach the FinalTime
  int Nsteps = ( int ) ceil( FinalTime / dt );
  dt = FinalTime / ( double ) Nsteps;

  printf( "dt = %f, FinalTime = %f, Nsteps = %d\n", dt, FinalTime, Nsteps );  
  
  // interval: number of steps to take to print out a new solution
  int interval = max( ( int ) ceil( Nsteps / 10 ), 1 );  
  printf( "Interval = %d\n", interval );
  
  // RK4 for time integration step
  int NINT = mesh->rk4a.size(); // num RK steps
  time_setup = clock() - time_setup;
  cout << "Time for the set up is " << (float) time_setup / CLOCKS_PER_SEC << " seconds" << endl;
  clock_t time_iteration = clock();
  // outside loop to loop over all time steps
  for ( int i = 0; i < Nsteps; ++i ){
    // inside loop to loop over all RK4 steps
    for ( int INTRK = 0; INTRK < NINT; ++INTRK ){

      const dfloat fdt = ( dfloat ) dt;
      const dfloat fa  = ( dfloat ) mesh->rk4a[ INTRK ];
      const dfloat fb  = ( dfloat ) mesh->rk4b[ INTRK ];

      // entropy projection
      project( K, o_VNP, o_Qv, o_Qf );
      
      // getOccaArray( app, o_Qv, Qv );
      // cout << "Qv block: row: " <<Qv.rows() << " . col: "<<Qv.cols() << ". " << endl << Qv << endl << endl;

      // getOccaArray( app, o_Qf, Qf );
      // cout << "Qf block: row: " <<Qf.rows() << " . col: "<<Qf.cols() << ". " << endl << Qf << endl << endl;

      // compute the volume term with computing device
      volume( K, o_vgeo, o_gBTMx, o_gBTMy, o_QNr, o_QNs, o_PN, o_Qv, o_Qf, o_rhs );
      // MatrixXd rhs( Np * Nfields, K );
      // getOccaArray( app, o_rhs, rhs );
      // cout << "rhs block: row: " <<rhs.rows() << " . col: "<<rhs.cols() << ". " << endl << rhs << endl << endl;

      // compute the suface term with computing device
      surface( K, o_fgeo, o_mapPq, o_bcFlag, o_Lf, o_Qf, o_rhs ); 

      // getOccaArray( app, o_rhs, rhs );
      // cout << "rhs block: row: " <<rhs.rows() << " . col: "<<rhs.cols() << ". " << endl << rhs << endl << endl;

      // combine the surface and volume to update the right hand side of the equation
      update( K, fa, fb, fdt, o_vgeo,o_Vq, o_Q, o_Qv, o_rhs, o_res );
      // getOccaArray( app, o_Qv, Qv );
      // cout << "Qv block: row: " <<Qv.rows() << " . col: "<<Qv.cols() << ". " << endl << Qv << endl << endl;

      // getOccaArray( app, o_Q, Q );
      // cout << "Q block: row: " <<Q.rows() << " . col: "<<Q.cols() << ". " << endl << Q << endl << endl;
    }

    if ( i % interval == 0 ){
      printf( "on timestep %d out of %d\n", i, Nsteps ) ;
    }
  }
  // calculate the time for iterations
  time_iteration = clock() - time_iteration;
  cout << "Time for the loop is " << (float) time_iteration / CLOCKS_PER_SEC << " seconds" << endl;
  // get the solution from computing device back to CPU
  getOccaArray( app, o_Q, Q );
  
  // calculate the time for Linf and L2 error
  clock_t time_error = clock();
  double Linf_err, L2_err;
  SWE_sol_err( N, K1D, g, FinalTime, mesh, Q, btm, &SWESOLUTION, Linf_err, L2_err );
  time_error = clock() - time_error;
  clock_t time_total = time_setup + time_iteration + time_error;
  cout << "Time for the Linf and L2 error calculation is " << (float) time_error / CLOCKS_PER_SEC;
  cout << " seconds" << endl;
  // // calculate the time for total
  // cout << "total time is " << (float) time_total / CLOCKS_PER_SEC << " seconds" << endl;
#if CAL_RES
  SWE_fileio( Lx, Ly, N, K1D, K, Nsteps, 
              curve, CFL, FinalTime, dt, L2_err, Linf_err, 
              time_setup, time_iteration, time_error, time_total );
#endif

#if PRINT_ERR
  SWE_error_output( mesh, Q, Vq, btm, g, FinalTime, CFL, N, K1D, &SWESOLUTION );
#endif

  return 0;
  
}