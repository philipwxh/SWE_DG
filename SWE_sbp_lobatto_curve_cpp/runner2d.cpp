#include <cstdlib>
#include <string>
#include <iostream>

using namespace std;
int main( int argc, char **argv ) {
  string program = "./main2d ";
  if ( argc == 2 ){
    program =  string( argv[1] );
    cout << "calling program "<< program << endl;
  }
	for( int N = 1; N <= 5; N += 1 ){
		for( int K1D = 8; K1D <= 128; K1D *= 2  ){
			for( double CFL = .125; CFL >= 0.125; CFL /= 2 ){
				for( double FT = 0.5; FT <= 0.5; FT += 0.5 ){
					for( double curve = 0.1; curve <= 0.1; curve += 0.2 ){
						string SWE_Call = program + " " + to_string( N )+ " ";
						SWE_Call += to_string( K1D ) + " " + to_string( CFL ) + " ";
						SWE_Call += to_string( FT ) + " " + to_string( curve );
						system( SWE_Call.c_str() );
					}
				}
			}
		}
	}
   return 0;
}