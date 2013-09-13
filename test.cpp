#include "euler_curve.h"

int main(int argc, char *argv[]){
	
	
	double k,l,g=1000;
	EulerCurve e;
	e.SolveEuler( 0, 0, PI/2, -1*2, 8, PI/2, 200, k, l, g);
	std::cout << k <<" "<< l <<" "<< g <<std::endl;
	
return 0;
}
