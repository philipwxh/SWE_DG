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
#define DEBUG 0
#define TOL_ZERO 1.e-14
#define PRINT_ERR 0
#define CAL_RES 0

int main( int argc, char **argv ){
  clock_t time_setup = clock();
  // N: default the degree of the polynomial
  int N = 3;
  // K1D: default number of elements in one direction.
  //      in 2D, there will be 16 by 16 quadrilateral or 16 by 16*2 triangles
  int K1D = 8;

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
  //    when a is large, there might be distorted curved mesh
  //    always verifiy the mesh!
  double curve = .2; // curved warping
  
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
  TriMesh2d( mesh, 2 * K1D, K1D );
  // [-2,2] x [-1,1] for vortex
  // Lx: length in x direction
  // Ly: length in y direction
  double Lx = 10;
  double Ly = 5;
  // VX: x coordination of elements vertices
  // VY: y coordination of elements vertices
  mesh->VX = ( mesh->VX ) * Lx;
  mesh->VY = ( mesh->VY ) * Ly;

  // ============ physics independent stuff ===========

  printf( "N = %d, K = %d\n", N, mesh->K );
  
  // initialize reference triangle element
  InitRefTri_sbp( mesh, N );

  // make physical nodes
  // low order mapping
  MapTriNodes( mesh );
  // curved the node by sin and cos function
  // mesh->x = mesh->x.array() + curve * ( 1.0 / 4.0 * PI * mesh->x.array() ).cos() 
  //                                 * ( PI * mesh->y.array() ).sin();
  // mesh->y = mesh->y.array() + curve * ( 1.0 / 2.0 * PI * mesh->y.array() ).cos() 
  //                                 * ( PI * mesh->x.array() ).sin();

  MatrixXd dx = ( PI / Lx * mesh->x.array() ).cos() * 
                ( 1.5 * PI / Ly * mesh->y.array() ).cos();
  mesh->x = mesh->x + Lx * curve * dx;
  MatrixXd dy = ( 2 * PI / Lx * mesh->x.array() ).sin() * 
                ( PI / Ly * mesh->y.array() ).cos();
  mesh->y = mesh->y + Ly * curve * dy;
 
  // initialize geometric factors
  GeometricFactors2d( mesh );

  // initialize normals in 2D
  Normals2d( mesh );

  // makes EToE, EToF
  ConnectElems( mesh );

  // data
  // Nfields: for 2D SWE, we have for variables to solve
  int Nfields = 3;
  mesh->Nfields = Nfields;
  // K: number of elements in the entire domain.
  int K = mesh->K;
  // Np: number of volume points in one element
  int Np = mesh->Np;
  // Nq: number of quadrature points in one element
  int Nq = mesh->Nq;    
  // NfpNfaces: number of face points multiplied by number of faces in one elements
  int NfqNfaces = mesh->Nfq * mesh->Nfaces;
  // NqT: total number of quadrature points in one element
  int NqT = Nq + NfqNfaces;

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
    // cout << "mapPq block: row: " << mapPq.rows() << " . col: "<< mapPq.cols() << ". " << endl << mapPq << endl << endl;

  // builds periodicity into the node maps
  double DX = mesh->VX.maxCoeff() - mesh->VX.minCoeff();
  double DY = mesh->VY.maxCoeff() - mesh->VY.minCoeff();
  // create the periodic boundary condition for 2D map
  // and store the face nodes connections for the boundary nodes
  MakeNodeMapsPeriodic2d( mesh, xf, yf, DX, DY, mapPq );
  // cout << "mapPq block: row: " << mapPq.rows() << " . col: "<< mapPq.cols() << ". " << endl << mapPq << endl << endl;

  // =================== set initial condition

  // extract the mesh points
  // x: x coordinates of the mesh points
  // y: y coordinates of the mesh points
  MatrixXd x = mesh->x;
  MatrixXd y = mesh->y;  
 
  // h:   SWE variable water height
  // u:   SWE variable velocity in x direction
  // v:   SWE wariable velocity in y direction
  // btm: SWE bottom geomtry
  MatrixXd h, u, v, btm;  

  // Vq: interpolate polynomial cofficient to volume quadrature points
  MatrixXd Vq = mesh->Vq;

  double time = 0.0;
  // g: gravitational constant for SWE
  //    might not be 9.8
  //    when g is negative, SWESOLUTION should initialize it problem specific value
  dfloat g = -1.0;
  // compute the initial condition  
  SWESOLUTION( Vq * x, Vq * y, time, g, h, u, v, btm );
  // cout << "h block: row: " << h.rows() << " . col: "<< h.cols() << ". " << endl << h << endl << endl;
  // hu: h * u
  MatrixXd hu = h.array() * u.array();
  // hv: h * v
  MatrixXd hv = h.array() * v.array();
  // stack all variable in to one matrix Q
  MatrixXd Q( ( Nfields + 1 ) * Nq, K );
  Q << h, hu, hv, btm;


  // ========================== set up OCCA application

  // create an App object to setup OCCA
  App *app = new App;
  
  // app->device.setup("mode: 'Serial'");
  // select the device we have
  app->device.setup( "mode: 'CUDA', device_id: 0" );
  cout << app->device.mode();
  // app->props["compiler_flags"] = "-arch=sm_30";


  // define macros for OCCA kernels
  // p_Np: number of dims
  app->props[ "defines/p_Np" ] = mesh->Np; 
  // p_Nq: number of quadrature points
  app->props[ "defines/p_Nq" ] = Nq;
  // p_NfqNfaces: number of face quadrature points in one face in one element
  app->props[ "defines/p_Nfq" ] = mesh->Nfq; 
  // p_NqT: number of total face and volume quadrature points in one element 
  app->props[ "defines/p_NfqNfaces" ] = NfqNfaces; 
  // p_NqT: number of total face and volume quadrature points in one element 
  app->props[ "defines/p_NqT" ] = NqT; 
  // p_T: max of total face quadrature points and dimension
  app->props[ "defines/p_T" ] = max( NfqNfaces, mesh->Np );
  
  cout << "Np = " << mesh->Np << "; Nq = " << Nq << ";NfqNfaces = " << NfqNfaces << "; NqT = " << NqT << endl;
  // number of geometric terms
  // Nvgeo: rxJ, ryJ, sxJ, syJ, J
  int Nvgeo = 5; 
  // Nsgeo: nxJ,nyJ, sJ
  int Nfgeo = 3; 
  app->props[ "defines/p_Nvgeo" ] = Nvgeo;
  app->props[ "defines/p_Nfgeo" ] = Nfgeo;

  int KblkP = 1; int KblkV = 1; int KblkS = 1; int KblkU = 1;
  // number of elements in one group for project kernel
  app->props[ "defines/p_KblkP" ] = KblkP;
  // number of elements in one group for volume kernel
  app->props[ "defines/p_KblkV" ] = KblkV;
  // number of elements in one group for surface kernel
  app->props[ "defines/p_KblkS" ] = KblkS;
  // number of elements in one group for update kernel
  app->props[ "defines/p_KblkU" ] = KblkU;

  // 1.f is single precision floating point
  // 1.0 is double precision floating point
  app->props[ "defines/p_Nfields" ] = Nfields;
  if (sizeof( dfloat ) == 4 ){
    app->props[ "defines/USE_DOUBLE" ] = 0;
    app->props[ "defines/p_tau" ]      = ( float )  tau;
    app->props[ "defines/p_g" ]        = ( float )  g;
    app->props[ "defines/p_g_half" ]   = ( float ) 0.5 * g;
  }else{
    app->props[ "defines/USE_DOUBLE" ] = 1;
    app->props[ "defines/p_tau" ]      = ( double ) tau;
    app->props[ "defines/p_g" ]        = ( double ) g;
    app->props[ "defines/p_g_half" ]   = ( double ) 0.5 * g;
  }

#if DEBUG
  occa::memory o_rt;
  // for right hand side test only
  MatrixXd rt( K, 1 );
  rt.fill( 0.0 );
  occa::memory o_rhstest;
  setOccaArray( app, rt, o_rt );
  app->props[ "defines/DEBUG" ] = 1;
#endif

  // build occa kernels  
  string path = "okl/SWE2DCurve_sbp.okl";
  // volume  : Calculate the volume contribution for SWE
  // surface : Calculate the surface flux contribution for SWE
  // update  : Multiple right hand side by Mass matrix inverse and perform RK4 over time
  // project : Calculate the entropy projection for SWE
  //testing
  occa::kernel volume, surface, update, rhstest;
  volume = app->device.buildKernel( path.c_str(), "volume", app->props );
  surface = app->device.buildKernel( path.c_str(), "surface", app->props );
  update  = app->device.buildKernel( path.c_str(), "update" , app->props );
#if DEBUG
  rhstest = app->device.buildKernel( path.c_str(), "rhstest", app->props );
#endif
  

  // ================= construct discretization matrices
  // Vq: interpolate polynomial cofficient to volume quadrature points
  // MatrixXd Vq = mesh->Vq;
  // cout << "Vq block: row: " << Vq.rows() << " . col: "<< Vq.cols() << ". " << endl << Vq << endl << endl;
  // Vf: interpolate polynomial cofficient to face quadrature points
  MatrixXd Vf = mesh->Vf;
  // M: mass matrix for reference element
  MatrixXd M = Vq.transpose() * mesh->wq.asDiagonal() * Vq;
  // VqTW: Vq^T * diag( wq )
  // wq:   volume quadrature weights
  MatrixXd VqTW = Vq.transpose() * mesh->wq.asDiagonal();
  // VfTWf: Vf^T * diag( wf )
  // wf:    face quadrature weights
  MatrixXd VfTWf = Vf.transpose() * mesh->wf.asDiagonal();
  // Pq:    inv( M ) * VqTW
  //        quadrature based L2 projector to degree N polynomial
  MatrixXd Pq = mldivide( M, VqTW );
  // Lf:    Lifting matrtix, 
  //        not used in curved case because we need different one for each element
  // MatrixXd Lf = mldivide( M, VfTWf );

  // print statement left here for furture debugging purpose
  // can be used to print the above matrices for correctness checking
  // cout << "VfTWf block: row: " << VfTWf.rows() << " . col: "<< VfTWf.cols() << ". " << endl << VfTWf << endl << endl;

  // hybridized SBP ops
  MatrixXd Ef = Vf * Pq;
  MatrixXd Wf = mesh->wf.asDiagonal();
  // nrJ: [0*e; e; -e]; e = ones(length(number of face quadrature points in one face))
  MatrixXd Br = Wf * mesh->nrJ.asDiagonal();
  // nsJ: [-e ; e; 0*e];
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
  // VN: interpolate coefficient to both volume and face quadrature points
  MatrixXd VN( NqT, Np );
  VN << Vq,Vf;
  // VNT: VN^T
  // MatrixXd VNT = VN.transpose();
  // VNP: not used in curved meshes since Pq is different for each element
  // MatrixXd VNP = VN * Pq;
  // PN: not used in curved meshes since mass matrix is different for each element
  // MatrixXd PN = mldivide( M, VNT ); 

  cout << "to compute sbp operators" <<endl;
  MatrixXd E_sbp = MatrixXd::Zero( mesh->rf.size(), mesh->rq.size() );
  MatrixXi E_sbpi = MatrixXi::Zero( mesh->rf.size(), NfqNfaces - 3 );

  MatrixXi E_sbp_idx( mesh->rf.size(), 1 );
  // cout << "mesh->rq block: row: " << mesh->rq.rows() << " . col: "<< mesh->rq.cols() << ". " << endl << mesh->rq << endl << endl;
  //   cout << "mesh->rf block: row: " << mesh->rf.rows() << " . col: "<< mesh->rf.cols() << ". " << endl << mesh->rf << endl << endl;

  // cout << "mesh->sq block: row: " << mesh->sq.rows() << " . col: "<< mesh->sq.cols() << ". " << endl << mesh->sq << endl << endl;
  //   cout << "mesh->sf block: row: " << mesh->sf.rows() << " . col: "<< mesh->sf.cols() << ". " << endl << mesh->sf << endl << endl;


  for( int i = 0; i < mesh->rf.size(); ++i ){
    for( int j = 0; j < mesh->rq.size(); ++j){
      if( abs( mesh->rq( j ) - mesh->rf( i ) ) + 
          abs( mesh->sq( j ) - mesh->sf( i ) ) < 1e-10 ){
        E_sbp( i, j ) = 1;
        E_sbpi( i, j ) = 1;
        E_sbp_idx( i, 0 ) = j;
      }
    }
  }
  // cout << "E_sbp block: row: " << E_sbp.rows() << " . col: "<< E_sbp.cols() << ". " << endl << E_sbp << endl << endl;
  // cout << "E_sbp_idx block: row: " << E_sbp_idx.rows() << " . col: "<< E_sbp_idx.cols() << ". " << endl << E_sbp_idx << endl << endl;
  MatrixXd VN_sbp( mesh->rf.size() + mesh->rq.size(), mesh->rq.size() );
  VN_sbp.topLeftCorner( mesh->rq.size(), mesh->rq.size() )  = MatrixXd::Identity( mesh->rq.size(), mesh->rq.size() );
  VN_sbp.bottomLeftCorner( mesh->rf.size(), mesh->rq.size() )  = E_sbp;

  MatrixXd QNr_sbp = VN_sbp.transpose() * QNr * VN_sbp;
  MatrixXd QNs_sbp = VN_sbp.transpose() * QNs * VN_sbp;
  MatrixXd Pf = E_sbp.transpose() * mesh->wf.asDiagonal();


 // ================= set node maps + boundary condition
  
  
  //bcFlag: Boundary condition flag. 0 for periodic boundary condition
  MatrixXi bcFlag( xf.rows(), xf.cols() );  
  bcFlag.fill( 0 ); // no BCs

  occa::memory o_bcFlag, o_mapPq, o_E_sbp_idx;
  setOccaIntArray( app, bcFlag, o_bcFlag );

  MatrixXi fnodes_idx = MatrixXi::Zero( NfqNfaces - 3, K );;

  for( int i = 0; i < (NfqNfaces - 3) * K; ++i ){
    fnodes_idx( i ) = i;
  }
  // cout << "fnodes_idx block: row: " << fnodes_idx.rows() << " . col: "<< fnodes_idx.cols() << ". " << endl << fnodes_idx << endl << endl;

  MatrixXi fnodes_idx_sbp = E_sbpi * fnodes_idx;
  MatrixXi mapPq_sbp = fnodes_idx_sbp;

  for( int i = 0; i < NfqNfaces * K; ++i ){
    mapPq_sbp( i ) = fnodes_idx_sbp( mapPq( i ) ) ;
  }
 
  mapPq = mapPq_sbp;
  // cout << "mapPq block: row: " << mapPq.rows() << " . col: "<< mapPq.cols() << ". " << endl << mapPq << endl << endl;

  // correct mapPq for Nfields > 1  
  MatrixXi mapPqNfields( NfqNfaces, K );
  // loop over all face nodes in all elements
  
  for( int i = 0; i < NfqNfaces * K; ++i ){
    int idP = mapPq( i );
    // e: element of the neighbor node
    int e = idP / (NfqNfaces-3);
    // fidP: index of the neighbor node in its element
    int fidP = idP % (NfqNfaces-3);
    // store the first index of the variable 
    mapPqNfields( i ) = fidP + e * Nq * ( Nfields + 1 );
  }  
  mesh->mapPq = mapPqNfields;
  // cout << "mapPqNfields block: row: " << mapPqNfields.rows() << " . col: "<< mapPqNfields.cols() << ". " << endl << mapPqNfields << endl << endl;
  setOccaIntArray( app, mesh->mapPq, o_mapPq );  

  cout << "mapPq built"<< endl;
  
  // print statement left here for furture debugging purpose
  // can be used to print the above matrices for correctness checking
  // cout << "Pf block: row: " << Pf.rows() << " . col: "<< Pf.cols() << ". " << endl << Pf << endl << endl;

  // calculate interpolated geometric cofficients
  MatrixXd JJ   = Vq * mesh->J;
  MatrixXd rxJJ = Vq * ( mesh->rxJ );
  MatrixXd ryJJ = Vq * ( mesh->ryJ );
  MatrixXd sxJJ = Vq * ( mesh->sxJ );
  MatrixXd syJJ = Vq * ( mesh->syJ );
  MatrixXd WqJJ = mesh->wq.asDiagonal() * JJ;

  // M_inv: stores the mass matrix inverse for all elements.
  //        To get the mass matrix for the e th element:
  //        M_inv.middleRows( e * Np, Np )
  MatrixXd M_inv( K * Nq, Nq );
  // QNx:   stores the QNx for all elements
  //        To get the QNx for the e th element:
  //        QNx.middleRows( e * NqT, NqT )
  // MatrixXd QNx( K * NqT, NqT  );
  // Qny:   stores the QNy for all elements
  //        To get the QNx for the e th element:
  //        QNy.middleRows( e * NqT, NqT )
  // MatrixXd QNy( K * NqT, NqT  );
  // BTMq:  Bottom gemetry interpolated to quadrature points
  MatrixXd BTMq = btm;
  // BTMx:  Stores interoplated bottom geometry derivatives in x direction
  //        To get the BTMx for the e th element:
  //        BTMx.middleCols( e, 1 )
  MatrixXd gBTMx( Nq, K  );
  // BTMy:  Stores interoplated bottom geometry derivatives in y direction
  //        To get the BTMx for the e th element:
  //        BTMy.middleCols( e, 1 )
  MatrixXd gBTMy( Nq, K  );

  // Building elementwise Mass matrices inverses and BTMx and BTMy
  for( int e = 0; e < K; ++e ){
      // get information for element e
      // MatrixXd JJi = JJ.col( e ).asDiagonal();
      // MatrixXd WqJJi = mesh->wq.asDiagonal() * JJi;
      MatrixXd M_e = WqJJ.col( e ).asDiagonal();
      MatrixXd rxJJi = rxJJ.col( e ).asDiagonal();
      MatrixXd ryJJi = ryJJ.col( e ).asDiagonal();
      MatrixXd sxJJi = sxJJ.col( e ).asDiagonal();
      MatrixXd syJJi = syJJ.col( e ).asDiagonal();

      // build mass matrix, QNx and QNy for element e
      // MatrixXd M_e = WqJJi.array().;
      MatrixXd QNx_e = 0.5 * ( rxJJi * QNr_sbp + QNr_sbp * rxJJi + sxJJi * QNs_sbp + QNs_sbp * sxJJi );
      MatrixXd QNy_e = 0.5 * ( ryJJi * QNr_sbp + QNr_sbp * ryJJi + syJJi * QNs_sbp + QNs_sbp * syJJi );
      
      // stack information to larger matrices
      M_inv.middleRows( e * Nq, Nq ) = M_e.inverse();
      gBTMx.middleCols( e, 1 )  = g * QNx_e * BTMq.col( e );
      gBTMy.middleCols( e, 1 )  = g * QNy_e * BTMq.col( e );
      // QNx.middleRows( e * NqT, NqT  ) = QNx_e;
      // QNy.middleRows( e * NqT, NqT  ) = QNy_e;
  }
  // print statement left here for furture debugging purpose
  // can be used to print the above matrices for correctness checking
  // cout << "M_inv block: row: " << M_inv.rows() << " . col: "<< M_inv.cols() << ". " << endl << M_inv << endl << endl;
  // cout << "BTMx block: row: " << BTMx.rows() << " . col: "<< BTMx.cols() << ". " << endl << BTMx << endl << endl;

  MatrixXd QNr_SBP = .5 * ( QNr_sbp - QNr_sbp.transpose() );
  MatrixXd QNs_SBP = .5 * ( QNs_sbp - QNs_sbp.transpose() );
  // cout << "QNr_sbp block: row: " << QNr_sbp.rows() << " . col: "<< QNr_sbp.cols() << ". " << endl << QNr_sbp << endl << endl;
  // cout << "QNs block: row: " << QNs.rows() << " . col: "<< QNs.cols() << ". " << endl << QNs << endl << endl;
  
  // solution:
  // o_Q:  contains h, hu, hv at mesh points
  // o_Qv: contains h, hu, hv interpolated to volume quadrature points
  // o_Qf: contains h, hu, hv interpolated to face quadrature points
  occa::memory o_Q, o_Qv, o_Qf; 
  occa::memory o_vgeo, o_fgeo; // geofacs

  // MatrixXd vgeo( Nvgeo * NqT, K );  
  // vgeo << rxJJ, ryJJ, sxJJ, syJJ, JJ;
  // move volume geometric factors to occa device
  occa::memory o_rxJJ, o_ryJJ, o_sxJJ, o_syJJ, o_WqJJ; //o_JJ;
  setOccaArray( app, rxJJ, o_rxJJ );
  setOccaArray( app, ryJJ, o_ryJJ );
  setOccaArray( app, sxJJ, o_sxJJ );
  setOccaArray( app, syJJ, o_syJJ );
  setOccaArray( app, WqJJ, o_WqJJ );

  // setOccaArray( app, JJ  , o_JJ   );

  // move face geometric factors to occa device
  MatrixXd fgeo( Nfgeo * NfqNfaces, K );  
  fgeo << mesh->nxJ, mesh->nyJ, mesh->sJ; // already at quad pts
  setOccaArray( app, fgeo, o_fgeo );

  // calculated right hand side timestep stuff
  occa::memory o_rhs, o_res, o_rhsv;
  setOccaArray( app, MatrixXd::Zero( Nq * Nfields, K ), o_rhs );
  setOccaArray( app, MatrixXd::Zero( Nq * Nfields, K ), o_res );  
  setOccaArray( app, MatrixXd::Zero( Nq * Nfields, K ), o_rhsv ); 

  // occa operators
  occa::memory o_M_inv, o_wq, o_QNr, o_QNs, o_Pf, o_Vq, o_VN; 
  // operators not used to curved mesh:
  // occa::memory o_Lf, o_PN, o_VNP
  occa::memory o_gBTx, o_gBTy; // precomputed data

  // set solution arrays
  setOccaArray( app, Q, o_Q );
  // MatrixXd Qv( ( Nfields + 1 ) * Nq, K );
  // Qv << Vq * h, Vq * hu, Vq * hv, Vq * btm;
  // setOccaArray( app, Qv, o_Qv );
  // // cout << "Qv block: row: " << Qv.rows() << " . col: "<< Qv.cols() << ". " << endl << Qv << endl << endl;
  // MatrixXd Qf( ( Nfields + 1 ) * NfqNfaces, K );
  // Qf << Vf * h, Vf * hu, Vf * hv, Vf * btm;
  // setOccaArray( app, Qf, o_Qf );
  // cout << "Qf block: row: " << Qf.rows() << " . col: "<< Qf.cols() << ". " << endl << Qf << endl << endl;

  // set operators
  setOccaArray( app, M_inv, o_M_inv );
  setOccaArray( app, mesh->wq, o_wq );
  setOccaIntArray( app, E_sbp_idx, o_E_sbp_idx );
  // setOccaArray( app, Vq, o_Vq );
  setOccaArray( app, QNr_SBP, o_QNr );
  setOccaArray( app, QNs_SBP, o_QNs );
  setOccaArray( app, Pf, o_Pf );
  setOccaArray( app, VN, o_VN );
  setOccaArray( app, gBTMx, o_gBTx );
  setOccaArray( app, gBTMy, o_gBTy );
  // not used for curved mesh
  // setOccaArray( app, PN, o_PN );
  // setOccaArray( app, VNP, o_VNP );
  // setOccaArray( app, Lf, o_Lf );

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
  // outside loop to loop over all time steps
  time_setup = clock() - time_setup;
  cout << "Time for the set up is " << (float) time_setup / CLOCKS_PER_SEC << " seconds" << endl;
  clock_t time_iteration = clock();

  for ( int i = 0; i < Nsteps ; ++i ){
    // inside loop to loop over all RK4 steps
    for ( int INTRK = 0; INTRK < NINT; ++INTRK ){

      const dfloat fdt = ( dfloat ) dt;
      const dfloat fa  = ( dfloat ) mesh->rk4a[ INTRK ];
      const dfloat fb  = ( dfloat ) mesh->rk4b[ INTRK ];
      // cout << "Qv block: row: " << Qv.rows() << " . col: "<< Qv.cols() << ". " << endl << Qv << endl << endl;
      // entropy projection
      // project( K, o_M_inv, o_VN, o_WqJJ, o_Qv, o_Qf, o_rhsv_DEBUG );

      // getOccaArray( app, o_rhsv_DEBUG, rhsv_DEBUG );
      // cout << "rhsv_DEBUG block: row: " << rhsv_DEBUG.rows() << " . col: "<< rhsv_DEBUG.cols() << ". " << endl << rhsv_DEBUG << endl << endl;

      // getOccaArray( app, o_Qv, Qv );
      // cout << "Qv block: row: " << Qv.rows() << " . col: "<< Qv.cols() << ". " << endl << Qv << endl << endl;
      // getOccaArray( app, o_Qf, Qf );
      // cout << "Qf block: row: " << Qf.rows() << " . col: "<< Qf.cols() << ". " << endl << Qf << endl << endl;
      // return 0;

      // compute the volume term with computing device
      volume( K, o_rxJJ, o_ryJJ, o_sxJJ, o_syJJ,  o_gBTx, o_gBTy, o_QNr, o_QNs, o_Q, o_rhsv );

      // MatrixXd rhsv( Nq * Nfields, K );
      // getOccaArray( app, o_rhsv, rhsv );
      // cout << "rhsv block: row: " << rhsv.rows() << " . col: "<< rhsv.cols() << ". " << endl << rhsv << endl << endl;
      // return 0;
      // cout << "volume done"<<endl;
      // compute the suface term with computing device
      // getOccaArray( app, o_Q, Q );
      // cout << "Q block: row: " << Q.rows() << " . col: "<< Q.cols() << ". " << endl << Q << endl << endl;
      surface( K, o_E_sbp_idx, o_fgeo, o_mapPq, o_bcFlag, o_Pf, o_Q, o_rhsv ); 
      // cout << "surface done"<<endl;

      // MatrixXd rhsv( Nq * Nfields, K );
      // getOccaArray( app, o_rhsv, rhsv );
      // cout << "rhsv block: row: " << rhsv.rows() << " . col: "<< rhsv.cols() << ". " << endl << rhsv << endl << endl;
      // return 0;
#if DEBUG
    if( INTRK == 0 ){
      // right hand side test
      rhstest( K, o_M_inv, o_WqJJ, o_Q, o_rhsv, o_rt );
      getOccaArray( app, o_rt, rt );
      cout << "right hand side test is " << rt.sum() << endl;
    }
#endif

      // combine the surface and volume to update the right hand side of the equation
      update( K, fa, fb, fdt, o_M_inv, o_Q, o_rhsv, o_res );
      // MatrixXd rhs( Np * Nfields, K );
      // getOccaArray( app, o_rhsv, rhsv );
      // cout << "rhsv block: row: " << rhsv.rows() << " . col: "<< rhsv.cols() << ". " << endl << rhsv << endl << endl;

      // getOccaArray( app, o_Qv, Qv );
      // cout << "Qv block: row: " << Qv.rows() << " . col: " << Qv.cols() << ". " << endl << Qv << endl << endl;
      // getOccaArray( app, o_Q, Q );
      // cout << "Q block: row: " << Q.rows() << " . col: " << Q.cols() << ". " << endl << Q << endl << endl;
      // return 0;
      
    }
// #if DEBUG
//     // right hand side test
//     rhstest( K, o_M_inv, o_WqJJ, o_Q, o_rhsv, o_rt );
//     getOccaArray( app, o_rt, rt );
//     cout << "right hand side test is " << rt.sum() << endl;
// #endif
    // cout << "right hand side test is " << rt2.sum() << endl;
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
  SWE_sol_err_sbp( N, K1D, g, FinalTime, mesh, Q, btm, &SWESOLUTION, Linf_err, L2_err );
  time_error = clock() - time_error;
  clock_t time_total = time_setup + time_iteration + time_error;
  cout << "Time for the Linf and L2 error calculation is " << (float) time_error / CLOCKS_PER_SEC;
  cout << " seconds" << endl;
  // // calculate the time for total
  // cout << "total time is " << (float) time_total / CLOCKS_PER_SEC << " seconds" << endl;
#if CAL_RES
  // wrting the result to a file
  SWE_fileio( Lx, Ly, N, K1D, K, Nsteps, 
              curve, CFL, FinalTime, dt, L2_err, Linf_err, 
              time_setup, time_iteration, time_error, time_total );
#endif

#if PRINT_ERR
  // writing the error to a file
  SWE_error_output( mesh, Q, Vq, btm, g, FinalTime, CFL, N, K1D, &SWESOLUTION );
#endif

  return 0;
}
