#ifndef _EULER_CURVE_H_
#define _EULER_CURVE_H_



#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <cmath>

#define PI 3.14159265359


/* @Brief : This class is used to find the parameters of a euler curve 
 * 			between two points given the tangents at the points
 * 			Since there are infinite solutions so the total curvature is 
 * 			minimised, this is achieved by keeping the change <= 2 * pi 			
 * 			Also since no straight away numerical solution exists for the
 * 			euler spiral equations gradient descent is used get the best fit.
 * 			For starting the estimation the initial conditions are to be specified.
 * 			The following code estimate the initial parameters by fitting a Biarc
 * 			on the given data.
 * 			The Curvature in Euler curve is represented as K(s) = gamma*s + k
 * 			So the purpose is to find k, gamma and the length of the curve
 */

class EulerCurve{
	
private:
	
	int integral_steps;
	
	// TODO: Replace the following by some calculus algebra libraries for more efficiency
	/* @Brief : Calculates the fresnel integral of x using trapezoidal
	 * 			method of numerical integration
	 * 			Cos term = integration of cos(s^2) with 0 < s < x
	 * 			Sin term = integration of sin(s^2) with 0 < s < x
	 * 
	 * @params  [in] x       : The variable whose fresnel integral is to be found
	 * @params [out] costerm : The value of the cos term of the integration
	 * @params [out] sinterm : The value of the sin term of the integration
	 */ 
	int fresnel(double x, double &costerm, double &sinterm);
	
	/* @Brief : Calculates the minimum of the 4 given numbers
	 * 			Returns the minimum of the 4 numbers.
	 * @params  [in] a  : The first number
	 * @params  [in] b  : The Second number
	 * @params  [in] c  : The Third number
	 * @params  [in] d  : The Fourth number
	 */ 
	double min4(double a, double b, double c, double d);
	
	/* @Brief : Calculates the minimum of the 3 given numbers
	 * 			Returns the minimum of the 3 numbers.
	 * @params  [in] a  : The first number
	 * @params  [in] b  : The Second number
	 * @params  [in] c  : The Third number
	 */ 
	double min3(double a, double b, double c);
	
	// ---------------------------------------------------------------------------
	// 							Biarc Estimator
	// ---------------------------------------------------------------------------
	
	/* @Brief : This returns the angle b/w two vectors 
	 * 			( Sample configuration in complex space - y + ix)
	 * @params  [in] x  : The imaginary part of the vector
	 * @params  [in] y  : The real part of the vector
	 */ 
	double AngleOfVector(double x, double y);
	
	/* @Brief : This return the angle at the meeting point of two Biarcs
	 * 			
	 * @params  [in] x0      : Starting X co-ordinate for biarc formation
	 * @params  [in] y0      : Starting Y co-ordinate for biarc formation
	 * @params  [in] theta0  : The angle of the tangent at the starting point (Starting orientation)
	 * @params  [in] x2      : Ending X co-ordinate for biarc formation
	 * @params  [in] y2      : Ending Y co-ordinate for biarc formation
	 * @params  [in] theta2  : The angle of the tangent at the Ending point (End point orientation)
	 * @params  [in] k1      : The curvature at tha starting point
	 * @params  [in] k2      : The curvature at tha ending point
	 */ 
	double ComputeJoinTheta( double x0, double y0, double theta0, double x2, double y2, double theta2, double k1, double k2);
	
	/* @Brief : This returns the arc length for a curve with starting tagential angle theta0 to the 
	 * 			point where the angle of the tangent vector reaches theta1 with a curve using constant
	 * 			curvature k
	 * @params  [in] theta0  : The angle of the tangent at the starting point (Starting orientation)
	 * @params  [in] theta1  : The angle of the tangent at the Ending point (End point orientation)
	 * @params  [in] k       : The curvature of the curve ( Remains constant)
	 */ 
	double ComputeArcLength(double theta0, double theta1, double k);

	/* @Brief : Fits a Biarc b/w two points given the tangent vector at both the points.
	 * 			
	 * @params  [in] x0        : Starting X co-ordinate for biarc formation
	 * @params  [in] y0        : Starting Y co-ordinate for biarc formation
	 * @params  [in] x2        : Ending X co-ordinate for biarc formation
	 * @params  [in] y2        : Ending Y co-ordinate for biarc formation
	 * @params  [in] theta0    : The angle of the tangent at the starting point (Starting orientation)
	 * @params  [in] theta2    : The angle of the tangent at the Ending point (End point orientation)
	 * @params [out] estimateK : The estimate of the curvature from the Biarc fitting
	 * @params [out] estimatel : The estimate of the length of the Biarc curve
	 */ 
	int ComputeBiarcSolution(double x0, double y0, // Start point coordinates
							 double x2, double y2, // End point coordinates
							 double theta0, double theta2, // The initial and final angles respectively
							 double &estimateK, double &estimatel); // To be calculated



	// ---------------------------------------------------------------------------
	// 						Euler Curve Formation
	// ---------------------------------------------------------------------------

	/* @Brief : Function to apply gradient descent for calculation of the euler curve parameters
	 * 			This takes the initial values from the Biarc estimation
	 * @params  [in] x0        : Starting X co-ordinate
	 * @params  [in] y0        : Starting Y co-ordinate
	 * @params  [in] theta0    : The angle of the tangent at the starting point (Starting orientation)
	 * @params  [in] x2        : The required ending X co-ordinate
	 * @params  [in] y2        : The required ending Y co-ordinate
	 * @params  [in] theta2    : The required angle of the tangent at the Ending point (End point orientation)
	 * @params  [in] estimateK : The estimate of the curvature from the Biarc fitting
	 * @params  [in] estimatel : The estimate of the length of the Biarc curve
	 * @params  [in] iternum   : The number of iterations to be done in gradient descent
	 * @params [out] Kfin      : The curvature parameter of euler curve
	 * @params [out] Lfin      : The length of the euler curve
	 */ 
	int SolveIteratively(double x0, double y0, double theta0, 
						double x1, double y1, double theta2, 
						double estimateK, double estimateL,
						int iternum,
						double &Kfin, double &Lfin);

	/* @Brief : This function calculates the end points X,Y fromt the curve paramets k & l 
	 * 			and than return the error as the distance as the straight line distance b/w 
	 * 			(X,Y) and (x2,y2)
	 * @params  [in] x0        : Starting X co-ordinate
	 * @params  [in] y0        : Starting Y co-ordinate
	 * @params  [in] theta0    : The angle of the tangent at the starting point (Starting orientation)
	 * @params  [in] x2        : The required ending X co-ordinate
	 * @params  [in] y2        : The required ending Y co-ordinate
	 * @params  [in] theta2    : The required angle of the tangent at the Ending point (End point orientation)
	 * @params  [in] k         : The curvature parameter of euler curve
	 * @params  [in] l         : The length of the euler curve
	 */ 
	double ComputeError(double x0, double y0, double theta0, double x2, double y2, double theta2, double k, double l);
	
	/* @Brief : This function calculates the end points X,Y fromt the curve paramets k & l 
	 * 			and than return points X,Y
	 * @params  [in] a         : Starting X co-ordinate
	 * @params  [in] b         : Starting Y co-ordinate
	 * @params  [in] theta     : The angle of the tangent at the starting point (Starting orientation)
	 * @params  [in] k         : The curvature parameter of euler curve
	 * @params  [in] gamma     : the gamma as in K(s) = gamma*s + k
	 * @params  [in] l         : The length of the euler curve
	 * @params [out] x         : The X co-ordinate of the end point calculated from the given curve parameters
	 * @params [out] y         : The Y co-ordinate of the end point calculated from the given curve parameters
	 */ 
	int EulerSpiralEndPoint(double a, double b, double theta, double k, double gamma, double s, double &x, double &y);					
	
public:
	EulerCurve(){
		integral_steps=1000; // Setting the no. of divisions for numerical integraion of fresnel integration
	}
	
	/* @Brief : Function to get the euler cuve parameters from the given starting and end points
	 * 			along with the orientations at both points.
	 * @params  [in] x0        : Starting X co-ordinate
	 * @params  [in] y0        : Starting Y co-ordinate
	 * @params  [in] theta0    : The angle of the tangent at the starting point (Starting orientation)
	 * @params  [in] x2        : The required ending X co-ordinate
	 * @params  [in] y2        : The required ending Y co-ordinate
	 * @params  [in] theta2    : The required angle of the tangent at the Ending point (End point orientation)
	 * @params  [in] iternum   : The number of iterations to be done in gradient descent
	 * @params [out] Kfin      : The curvature parameter of euler curve
	 * @params [out] Lfin      : The length of the euler curve
	 */
	int SolveEuler( double x0, double y0, double theta0, 
					double x1, double y1, double theta2, 
					int iternum,
					double &Kfin, double &Lfin, double &Gfin);
	
	/* @Brief : Function to set the no. of divisions for numerical integraion of fresnel integration
	 * 			This is required as for large difference b/w (x2,y2) and (x0,y0) the default value
	 * 			1000 is not enough
	 * @params  [in] no_of_integral_steps    : No. of divisions for numerical integraion of fresnel integration
	 */
	void setFresenlIntegrationSteps(int no_of_integral_steps){
		integral_steps=no_of_integral_steps;
	}
};











// ---------------------------------------------------------------------------------------
// 								Function Definitions
// ---------------------------------------------------------------------------------------








int EulerCurve::fresnel(double x, double &costerm, double &sinterm){
	costerm=0;
	sinterm=0;
	costerm += ( cos(0) + cos( (x*x) ) ) /2;
	sinterm += ( sin(0) + sin( (x*x) ) ) /2;
	
	for(int i=0;i<integral_steps;i++){
		costerm += cos( ((x)*i/integral_steps)*((x)*i/integral_steps)  ) ;
		sinterm += sin( ((x)*i/integral_steps)*((x)*i/integral_steps)  ) ;
	}

	costerm = costerm * (x)/integral_steps;
	sinterm = sinterm * (x)/integral_steps;
return 0;
}

double EulerCurve::min4(double a, double b, double c, double d){

	if( a <= b && a <= c && a<=d )
		return a;
		
	if( b <= a && b <= c && b<=d )
		return b;

	if( c <= a && c <= b && c<=d )
		return c;

return d;
}

double EulerCurve::min3(double a, double b, double c){

	if( a <= b && a <= c)
		return a;
		
	if( b <= a && b <= c)
		return b;

return c;
}
	 
double EulerCurve::AngleOfVector(double x, double y){
	double angle = atan2(y,x);
	if( angle < 0 ){
		angle = angle + 2*PI;
	}
return angle;
}

double EulerCurve::ComputeJoinTheta( double x0, double y0, double theta0, double x2, double y2, double theta2, double k1, double k2){
	double sin_thetaJoin = ( k1*k2*(x2-x0) + k2*sin(theta0) -k1*sin(theta2) )/(k2-k1);
	double cos_thetaJoin = ( k1*k2*(y2-y0) + k2*cos(theta0) -k1*cos(theta2) )/(k2-k1);

return AngleOfVector(sin_thetaJoin,cos_thetaJoin);
}

double EulerCurve::ComputeArcLength(double theta0, double theta1, double k){
	double numerator = theta1-theta0;
	if( k < 0 &&  numerator > 0){
		numerator = numerator - 2*PI;
	}
	else if( k > 0 &&  numerator < 0){
		numerator = numerator + 2*PI;
	}
return numerator/k;
}



int EulerCurve::ComputeBiarcSolution(double x0, double y0, // Start point coordinates
						 double x2, double y2, // End point coordinates
						 double theta0, double theta2, // The initial and final angles respectively
						 double &estimateK, double &estimatel) // To be calculated
						 
{						 
	double L = sqrt( (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) );
	double psi = AngleOfVector( x2-x0 , y2-y0 );
	
	double smallAngle = 0.001;
	double smallCurvature = 0.001;
	int biarc_flag = 0;
	double L1 = -10;
	double L2 = -10;
	
	double k1, k2;
	
	double L3,L4,E1,E2;
	
	double check = (psi - (theta2-theta0)/2 );
	while (check > PI || check < -1*PI){
		if(check > PI){
			check = check - PI;
		}
		else if (check < (-1*PI)){
			check = check + PI;
		}
	}
	
	
	if( check > (-1*smallAngle) && check < (smallAngle) ){
		if( (psi-theta0) > (-1*smallAngle) && (psi-theta0) < (smallAngle) ){
			k1 = 0 ;
			k2 = k1 ;
			L1 = L/2;
			L2 = L1;
		}
		else{
			k1 = -2*sin((theta2-theta0)/2) / L ;
			k2 = k1 ;
			L1 = std::abs(theta0 * L * sin((PI/2 - theta0)/sin(2*theta0)));
			L2 = L1;
		}
		estimateK = k1;
		estimatel = L1+L2;
		
	}// if( check > (-1*smallAngle) && check < (smallAngle) )
	else{
		biarc_flag=1;
		k1 = -4*sin( (3*theta0 + theta2)/4 - psi )*cos( (theta2-theta0)/4 )/L;
		k2 =  4*sin( (theta0 + 3*theta2)/4 - psi )*cos( (theta2-theta0)/4 )/L;
		double theta_join = ComputeJoinTheta(x0,y0,theta0,x2,y2,theta2,k1,k2);
		
		if( k1 == 0){
			L1 = L* sin( (theta2+theta0)/2 - psi) / sin( (theta0-theta2)/2 );
		}
		else{
			L1 = ComputeArcLength(theta0,theta_join, k1);
		}
		
		if( k2 == 0){
			L2 = L* sin( (theta2+theta0)/2 - psi) / sin( (theta0-theta2)/2 );
		}
		else{
			L2 = ComputeArcLength(theta_join, theta2,k2);
		}
		
		double k3 = (4/L)*cos((3*theta0 + theta2)/4 - psi)*sin( (theta2-theta0)/4 );
		double k4 = (4/L)*cos((theta0 + 3*theta2)/4 - psi)*sin( (theta2-theta0)/4 );
		
		
		if( k3 == 0 && k4 == 0 ){
			estimatel = L;
			estimateK = 0;
		}
		else{
			if(k3 != k4){
				theta_join = ComputeJoinTheta(x0,y0,theta0,x2,y2,theta2,k3,k4);
				
				if( k3 == 0){
					L3 = L* sin( (theta2+theta0)/2 - psi) / sin( (theta2-theta0)/2 );
				}
				else{
					L3 = ComputeArcLength(theta0,theta_join , k3);
				}
				
				if( k4 == 0){
					L4 = L* sin( (theta2+theta0)/2 - psi) / sin( (theta0-theta2)/2 );
				}
				else{
					L4 = ComputeArcLength(theta_join, theta2 , k4);
				}
			}
		}
		
		E1 = (k2-k1)*(k2-k1);
		E2 = (k3-k4)*(k3-k4);
		
		
		if ( (L1+L2) < (L3+L4) ){
			estimateK = k1;
			estimatel = L1 + L2;
		}
		else{
			estimateK = k3;
			estimatel = L3 + L4;
		}
		
	}
return 0;
}

double EulerCurve::ComputeError(double x0, double y0, double theta0, double x2, double y2, double theta2, double k, double L){
	double ex, ey;
	double gamma = 2*(theta2-theta0 - k*L) /(L*L);
	EulerSpiralEndPoint(x0, y0, theta0, k, gamma, L, ex, ey);
return sqrt( (ex-x2)*(ex-x2) + (ey-y2)*(ey-y2) );
}


int EulerCurve::EulerSpiralEndPoint(double a, double b, double theta, double k, double gamma, double s, double &x, double &y){
	double epsilon = 0.00001;
	double constTerm = 0;
	double fc, fs, sc, ss;
	double cosTerm, sinTerm;
	double A,B;
	if( (gamma > 0 && gamma <epsilon ) || (gamma < 0 && gamma > -1*epsilon ) ){
		gamma = 0;
	}
	else{
		if(gamma > 0){
			fresnel((k+gamma*s)/(sqrt(PI*gamma)),fc,fs);
			fresnel((k)/(sqrt(PI*gamma)),sc,ss);
			A = fc - sc;
			B = fs - ss;
			cosTerm = cos(theta - k*k*0.5/gamma);
			sinTerm = sin(theta - k*k*0.5/gamma);
			constTerm = sqrt(PI/gamma);
			x = a + (constTerm) * (cosTerm * A - sinTerm * B);
			y = b + (constTerm) * (sinTerm * A + cosTerm * B);
		}
		if(gamma < 0){
			fresnel((k*-1+gamma*s*-1)/(sqrt(PI*gamma*-1)),fc,fs);
			fresnel((k*-1)/(sqrt(PI*gamma*-1)),sc,ss);
			A = fc - sc;
			B = (fs - ss)*-1;
			cosTerm = cos(theta - k*k*0.5/gamma);
			sinTerm = sin(theta - k*k*0.5/gamma);
			constTerm = sqrt(PI/gamma*-1);
			x = a + (constTerm) * (cosTerm * A - sinTerm * B);
			y = b + (constTerm) * (sinTerm * A + cosTerm * B);
		}
	}
	if(gamma == 0){
		if( k ==0 ){
			x= a + s*cos(theta);
			y= b + s*sin(theta);
		}
		else{
			constTerm = 0.1 / k ;
			x= a + (constTerm)    * ( sin(k*s+theta) - sin(theta));
			y= b + (-1*constTerm) * ( cos(k*s+theta) - cos(theta));
		}
	}
return 0;
}


int EulerCurve::SolveIteratively(double x0, double y0, double theta0, 
				     double x1, double y1, double theta2, 
				     double estimateK, double estimateL,
				     int iternum,
				     double &Kfin, double &Lfin)
{
	
	double error = ComputeError(x0,y0,theta0,x1,y1,theta2,estimateK,estimateL);
	double prevError = 1000;
	double errorStep = 0.1;
	double epsilon = 0.001;
	double epsilonError = 0.1;
	Kfin = estimateK;
	Lfin = estimateL;
	int i=0;
	
	double error0, error1, error2, error3;
	error0 = error1 = error2 = error3 = 999;
	do{
		if(i == (iternum-1) || error<epsilon){
			return 1;
		}
		
		error0 = ComputeError(x0,y0,theta0,x1,y1,theta2, Kfin + errorStep, Lfin );
		error1 = ComputeError(x0,y0,theta0,x1,y1,theta2, Kfin - errorStep, Lfin );
		error2 = ComputeError(x0,y0,theta0,x1,y1,theta2, Kfin            , Lfin + errorStep );
		int f=0;
		if( Lfin > errorStep ){
				f=1;
				error3 = ComputeError(x0,y0,theta0,x1,y1,theta2, Kfin, Lfin - errorStep );
		}
		if(f==1){
			error = min4 ( error0,error1,error2,error3 );
		}
		else{
			error = min3 ( error0,error1,error2 );
		}
		
		if(error >= prevError ){
			errorStep = errorStep/2;
		}
		else if(error == error0){
			Kfin = Kfin + errorStep;
		}
		else if(error == error1){
			Kfin = Kfin - errorStep;
		}
		else if(error == error2){
			Lfin = Lfin + errorStep;
		}
		else if(Lfin > errorStep){
			Lfin = Lfin - errorStep;
		}
		prevError = error;
		i++;
		printf("DEBUG (Line %d): Error : %lf %lf %lf %lf\n",__LINE__,error0,error1,error2,error3);
	}while ( i<=iternum );
	
	return 0;
}


int EulerCurve::SolveEuler( double x0, double y0, double theta0, 
				double x2, double y2, double theta2, 
				int iternum,
				double &Kfin, double &Lfin, double &Gfin){
	double ek, el;
	ComputeBiarcSolution(x0, y0,
						 x2, y2,
						 theta0, theta2,
						 ek, el);
						 
	SolveIteratively(x0, y0, theta0, 
				     x2, y2, theta2, 
				     ek, el,
				     iternum,
				     Kfin, Lfin);
	Gfin = 2 * (theta2 - theta0 - Kfin*Lfin) / (Lfin*Lfin); 
}



#endif
