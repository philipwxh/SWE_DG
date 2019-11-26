#include <cstdlib>
#include <string>
#include <iostream>

using namespace std;
int main( int argc, char **argv ) {
  string program = "./flux_test ";
  if ( argc == 2 ){
    program =  string( argv[1] );
    cout << "calling program "<< program << endl;
  }
  for( int N = 210; N <= 300; N += 10 ){
	string flux_test_Call = program + " " + to_string( N )+ " 1024 1 1";
	system( flux_test_Call.c_str() );	
  }
  return 0;
}