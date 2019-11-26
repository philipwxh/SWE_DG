#include "SWE_fileio.hpp"

void SWE_fileio( double Lx, double Ly, int N, int K1D, int K, int Nsteps, 
     double curve, double CFL, double FinalTime, double dt, double L2_err, double Linf_err, 
     clock_t time_setup, clock_t time_iteration, clock_t time_error, clock_t time_total ){  
 // write results to .csv file

  string filename;
  if( curve == 0 ){
    filename= "SWE_affine_result_Lx_";
  }else{
    filename= "SWE_curve_result_Lx_";
  }

  filename += to_string( Lx ) + "_Ly" + to_string( Ly ) + ".csv";
  ifstream file_test( filename.c_str() );
  ofstream outfile;
  // check if the file already exist
  bool file_exist = file_test.good();
  file_test.close();
  // if the file already exist, open to append
  if( file_exist ){
    outfile.open( filename.c_str(), ios_base::app );
  }else{
    // if the file does not exist, create and open to write
    outfile.open( filename.c_str() );
  }
  // check if the file opened correctly
  if( outfile ){
    // if the file does not exist, add a row of column description
    if( !file_exist ){
      outfile << "N,K1D,K,CFL,Final Time,Curve cofficient,dt,Nsteps,L2 Error, Linf_err,";
      outfile << "Setup Time,Iteration Time,L2 Error time,Total Time,Finished Time\n";
    }
    // write out all the result and parameter
    outfile << N << "," << K1D << ","  << K << "," << CFL << "," << FinalTime << ",";
    outfile << curve << "," << dt << "," << Nsteps << "," << L2_err << "," << Linf_err << ",";
    outfile << (float) time_setup / CLOCKS_PER_SEC << ",";
    outfile << (float) time_iteration / CLOCKS_PER_SEC << ",";
    outfile << (float) time_error / CLOCKS_PER_SEC << ",";
    outfile << (float) time_total / CLOCKS_PER_SEC << ",";
    time_t cur_time = clock();
    outfile << ctime( &cur_time );
  }else{
    cerr << "Cannot open file results.csv to write" << endl;
  }
  outfile.close();
}

void flux_test_fileio( int N, int K, clock_t flux_time_t1, clock_t flux_time_t2 ){  
 // write results to .csv file

  string filename = "flux_test_serial.csv";
  ifstream file_test( filename.c_str() );
  ofstream outfile;
  // check if the file already exist
  bool file_exist = file_test.good();
  file_test.close();
  // if the file already exist, open to append
  if( file_exist ){
    outfile.open( filename.c_str(), ios_base::app );
  }else{
    // if the file does not exist, create and open to write
    outfile.open( filename.c_str() );
  }
  // check if the file opened correctly
  if( outfile ){
    // if the file does not exist, add a row of column description
    if( !file_exist ){
      outfile << "N, K, flux_time_t1, flux_time_t2\n";
    }
    // write out all the result and parameter
    outfile << N << "," << K << ",";
    outfile << (float) flux_time_t1 / CLOCKS_PER_SEC << ",";
    outfile << (float) flux_time_t2 / CLOCKS_PER_SEC << endl;
  }else{
    cerr << "Cannot open file results.csv to write" << endl;
  }
  outfile.close();
}

void SWE_error_output( Mesh *mesh, MatrixXd &Q, MatrixXd &Vq, MatrixXd &btm,
                       double g, double FinalTime, double CFL, int N, int K1D,
     void ( *SWE_solution )( MatrixXd, MatrixXd, double, double &, MatrixXd &, MatrixXd &, MatrixXd &, MatrixXd &)){  
 // write results to .csv file
  MatrixXd h  = Q.middleRows( 0, mesh->Np );
  MatrixXd hu = Q.middleRows( mesh->Np, mesh->Np );
  MatrixXd hv = Q.middleRows( 2 * mesh->Np, mesh->Np );

  MatrixXd hex, uex, vex;
  MatrixXd xq   = Vq * mesh->x;
  MatrixXd yq   = Vq * mesh->y;
  MatrixXd btmq = Vq * btm;
  
  // compute the exact solution
  SWE_solution( xq, yq, FinalTime, g, hex, uex, vex, btmq );
  
  MatrixXd huex = hex.array() * uex.array();
  MatrixXd hvex = hex.array() * vex.array();
  // evaluate the L2 error
  MatrixXd SWE_hError  = hex  - Vq * h;
  MatrixXd SWE_huError = huex - Vq * hu;
  MatrixXd SWE_hvError = hvex - Vq * hv;
   
  string error_filename = "SWE_affine_h_" + to_string( N ) + "_" + to_string( K1D );
  error_filename += "_" + to_string( CFL ) + "_" + to_string( FinalTime );
  ofstream SWE_error_file( error_filename.c_str() );
  for( int i = 0; i < SWE_hError.rows(); ++i ){
    for( int j = 0; j < SWE_hError.cols(); ++j ){
      SWE_error_file << SWE_hError( i, j ) << ",";
    }
    SWE_error_file << "\n";
  }
  SWE_error_file.close();

  error_filename = "SWE_affine_hu_" + to_string( N ) + "_" + to_string( K1D );
  error_filename += "_" + to_string( CFL ) + "_" + to_string( FinalTime );
  SWE_error_file.open( error_filename.c_str() );
  for( int i = 0; i < SWE_huError.rows(); ++i ){
    for( int j = 0; j < SWE_huError.cols(); ++j ){
      SWE_error_file << SWE_huError( i, j ) << ",";
    }
    SWE_error_file << "\n";
  }
  SWE_error_file.close();

  error_filename = "SWE_affine_hv_" + to_string( N ) + "_" + to_string( K1D );
  error_filename += "_" + to_string( CFL ) + "_" + to_string( FinalTime );
  SWE_error_file.open( error_filename.c_str() );
  for( int i = 0; i < SWE_hvError.rows(); ++i ){
    for( int j = 0; j < SWE_hvError.cols(); ++j ){
      SWE_error_file << SWE_hvError( i, j ) << ",";
    }
    SWE_error_file << "\n";
  }
  SWE_error_file.close();
}