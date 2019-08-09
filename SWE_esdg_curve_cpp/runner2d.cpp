#include <cstdlib>
#include <string>
using namespace std;
int main() {
	for( int N = 1; N <= 6; N += 1 ){
		for( int K1D = 4; K1D <= 128; K1D *= 2  ){
			for( double CFL = .125; CFL >= 0.125; CFL /= 2 ){
				for( double FT = 0.5; FT <= 0.5; FT += 0.5 ){
					for( double curve = 0.1; curve <= 0.1; curve += 0.2 ){
						string SWE_Call = "./main2d " + to_string( N )+ " ";
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