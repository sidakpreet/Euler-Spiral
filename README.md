Euler-Spiral
============

Source:
Euler Spiral for Shape Completion
BENJAMIN B. KIMIA, ILANA FRANKEL AND ANA-MARIA POPESCU
Division of Engineering, Brown University, Providence, R1 02912, USA

Takes the initial and final postion+orientataion. 
Returns the Euler Spiral parameters k, gamma and curve length with the equations:
    Curvature K(s) = gamma*s + k (Here s is the arc length)

/* @Brief : This is used to find the parameters of a euler spiral 
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
