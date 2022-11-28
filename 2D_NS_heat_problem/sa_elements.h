#pragma once
#include "utils.h"
#include "geom3.h"
#include <iostream>
#include <tuple>
#include <vector>
#include <cmath>


// First Row-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inline
double A00_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m,  
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
	return f0*((1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*f0
	        +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*f0x+dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*f0y);
}	
inline
double A00_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m,  
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
	return f0*((1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*f1
	        +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*f1x+dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*f1y);
}	
inline
double A00_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m,  
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
	return f0*((1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*f2
	        +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*f2x+dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*f2y);
}	
inline
double A00_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m,  
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
	return f1*((1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*f0
	        +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*f0x+dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*f0y);
}	
inline
double A00_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m,  
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
	return f1*((1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*f1
	        +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*f1x+dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*f1y);
}	
inline
double A00_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m,  
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
	return f1*((1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*f2
	        +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*f2x+dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*f2y);
}	
inline
double A00_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m,  
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
	return f2*((1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*f0
	        +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*f0x+dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*f0y);
}	
inline
double A00_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m,  
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
	return f2*((1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*f1
	        +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*f1x+dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*f1y);
}	
inline
double A00_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m,  
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
	return f2*((1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*f2
	        +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*f2x+dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*f2y);
}	

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A01_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f0*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g0+(r0_m*f0+r1_m*f1+r2_m*f2)*(g0x));
}	
inline
double A01_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f0*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g1+(r0_m*f0+r1_m*f1+r2_m*f2)*(g1x));
}	
inline
double A01_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f0*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g2+(r0_m*f0+r1_m*f1+r2_m*f2)*(g2x));
}	
inline
double A01_03(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f0*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g3+(r0_m*f0+r1_m*f1+r2_m*f2)*(g3x));
}	
inline
double A01_04(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f0*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g4+(r0_m*f0+r1_m*f1+r2_m*f2)*(g4x));
}	
inline
double A01_05(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f0*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g5+(r0_m*f0+r1_m*f1+r2_m*f2)*(g5x));
}	
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
inline
double A01_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f1*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g0+(r0_m*f0+r1_m*f1+r2_m*f2)*(g0x));
}	
inline
double A01_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f1*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g1+(r0_m*f0+r1_m*f1+r2_m*f2)*(g1x));
}	
inline
double A01_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f1*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g2+(r0_m*f0+r1_m*f1+r2_m*f2)*(g2x));
}	
inline
double A01_13(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f1*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g3+(r0_m*f0+r1_m*f1+r2_m*f2)*(g3x));
}	
inline
double A01_14(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f1*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g4+(r0_m*f0+r1_m*f1+r2_m*f2)*(g4x));
}	
inline
double A01_15(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f1*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g5+(r0_m*f0+r1_m*f1+r2_m*f2)*(g5x));
}	
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
inline
double A01_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f2*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g0+(r0_m*f0+r1_m*f1+r2_m*f2)*(g0x));
}	
inline
double A01_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f2*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g1+(r0_m*f0+r1_m*f1+r2_m*f2)*(g1x));
}	
inline
double A01_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f2*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g2+(r0_m*f0+r1_m*f1+r2_m*f2)*(g2x));
}	
inline
double A01_23(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f2*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g3+(r0_m*f0+r1_m*f1+r2_m*f2)*(g3x));
}	
inline
double A01_24(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f2*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g4+(r0_m*f0+r1_m*f1+r2_m*f2)*(g4x));
}	
inline
double A01_25(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f2*dt*((r0_m*f0x+r1_m*f1x+r2_m*f2x)*g5+(r0_m*f0+r1_m*f1+r2_m*f2)*(g5x));
}	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A02_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f0*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g0+(r0_m*f0+r1_m*f1+r2_m*f2)*(g0y));
}	
inline
double A02_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f0*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g1+(r0_m*f0+r1_m*f1+r2_m*f2)*(g1y));
}	
inline
double A02_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f0*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g2+(r0_m*f0+r1_m*f1+r2_m*f2)*(g2y));
}	
inline
double A02_03(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f0*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g3+(r0_m*f0+r1_m*f1+r2_m*f2)*(g3y));
}	
inline
double A02_04(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f0*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g4+(r0_m*f0+r1_m*f1+r2_m*f2)*(g4y));
}	
inline
double A02_05(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f0*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g5+(r0_m*f0+r1_m*f1+r2_m*f2)*(g5y));
}	
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
inline
double A02_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f1*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g0+(r0_m*f0+r1_m*f1+r2_m*f2)*(g0y));
}	
inline
double A02_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f1*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g1+(r0_m*f0+r1_m*f1+r2_m*f2)*(g1y));
}	
inline
double A02_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f1*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g2+(r0_m*f0+r1_m*f1+r2_m*f2)*(g2y));
}	
inline
double A02_13(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f1*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g3+(r0_m*f0+r1_m*f1+r2_m*f2)*(g3y));
}	
inline
double A02_14(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f1*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g4+(r0_m*f0+r1_m*f1+r2_m*f2)*(g4y));
}	
inline
double A02_15(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f1*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g5+(r0_m*f0+r1_m*f1+r2_m*f2)*(g5y));
}	
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
inline
double A02_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f2*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g0+(r0_m*f0+r1_m*f1+r2_m*f2)*(g0y));
}	
inline
double A02_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f2*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g1+(r0_m*f0+r1_m*f1+r2_m*f2)*(g1y));
}	
inline
double A02_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f2*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g2+(r0_m*f0+r1_m*f1+r2_m*f2)*(g2y));
}	
inline
double A02_23(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f2*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g3+(r0_m*f0+r1_m*f1+r2_m*f2)*(g3y));
}	
inline
double A02_24(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f2*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g4+(r0_m*f0+r1_m*f1+r2_m*f2)*(g4y));
}	
inline
double A02_25(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return f2*dt*((r0_m*f0y+r1_m*f1y+r2_m*f2y)*g5+(r0_m*f0+r1_m*f1+r2_m*f2)*(g5y));
}	
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A10_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m, double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f0
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g0*f0x;
}	 
inline
double A10_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f1
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g0*f1x;
}	 
inline
double A10_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f2
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g0*f2x;
}	 
inline
double A10_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f0
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g1*f0x;
}	 
inline
double A10_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f1
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g1*f1x;
}
inline
double A10_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f2
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g1*f2x;
}
inline
double A10_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f0
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g2*f0x;
}
inline
double A10_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f1
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g2*f1x;
}
inline
double A10_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f2
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g2*f2x;
}
inline
double A10_30(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f0
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g3*f0x;
}
inline
double A10_31(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f1
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g3*f1x;
}
inline
double A10_32(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f2
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g3*f2x;
}
inline
double A10_40(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f0
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g4*f0x;
}
inline
double A10_41(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f1
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g4*f1x;
}
inline
double A10_42(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f2
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g4*f2x;
}
inline
double A10_50(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f0
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g5*f0x;
}
inline
double A10_51(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f1
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g5*f1x;
}
inline
double A10_52(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0x+r1_m*f1x+r2_m*f2x))*f2
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g5*f2x;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A11_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+mu*dt*(2*g0x*g0x+g0y*g0y)+2*dt*(mu+lambda)*(g0x*g0x);
}
inline
double A11_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+mu*dt*(2*g1x*g0x+g1y*g0y)+2*dt*(mu+lambda)*(g1x*g0x);
}	
inline
double A11_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+mu*dt*(2*g2x*g0x+g2y*g0y)+2*dt*(mu+lambda)*(g2x*g0x);
}	
inline
double A11_03(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+mu*dt*(2*g3x*g0x+g3y*g0y)+2*dt*(mu+lambda)*(g3x*g0x);
}	
inline
double A11_04(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+mu*dt*(2*g4x*g0x+g4y*g0y)+2*dt*(mu+lambda)*(g4x*g0x);
}	
inline
double A11_05(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+mu*dt*(2*g5x*g0x+g5y*g0y)+2*dt*(mu+lambda)*(g5x*g0x);
}
//#############################
inline
double A11_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+mu*dt*(2*g0x*g1x+g0y*g1y)+2*dt*(mu+lambda)*(g0x*g1x);
}	
inline
double A11_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+mu*dt*(2*g1x*g1x+g1y*g1y)+2*dt*(mu+lambda)*(g1x*g1x);
}	
inline
double A11_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+mu*dt*(2*g2x*g1x+g2y*g1y)+2*dt*(mu+lambda)*(g2x*g1x);
}	
inline
double A11_13(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+mu*dt*(2*g3x*g1x+g3y*g1y)+2*dt*(mu+lambda)*(g3x*g1x);
}	
inline
double A11_14(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+mu*dt*(2*g4x*g1x+g4y*g1y)+2*dt*(mu+lambda)*(g4x*g1x);
}	
inline
double A11_15(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+mu*dt*(2*g5x*g1x+g5y*g1y)+2*dt*(mu+lambda)*(g5x*g1x);
}
//#######################
inline
double A11_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+mu*dt*(2*g0x*g2x+g0y*g2y)+2*dt*(mu+lambda)*(g0x*g2x);
}	
inline
double A11_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+mu*dt*(2*g1x*g2x+g1y*g2y)+2*dt*(mu+lambda)*(g1x*g2x);
}	
inline
double A11_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+mu*dt*(2*g2x*g2x+g2y*g2y)+2*dt*(mu+lambda)*(g2x*g2x);
}	
inline
double A11_23(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+mu*dt*(2*g3x*g2x+g3y*g2y)+2*dt*(mu+lambda)*(g3x*g2x);
}	
inline
double A11_24(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+mu*dt*(2*g4x*g2x+g4y*g2y)+2*dt*(mu+lambda)*(g4x*g2x);
}	
inline
double A11_25(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+mu*dt*(2*g5x*g2x+g5y*g2y)+2*dt*(mu+lambda)*(g5x*g2x);
}

//#######################
inline
double A11_30(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+mu*dt*(2*g0x*g3x+g0y*g3y)+2*dt*(mu+lambda)*(g0x*g3x);
}	
inline
double A11_31(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+mu*dt*(2*g1x*g3x+g1y*g3y)+2*dt*(mu+lambda)*(g1x*g3x);
}	
inline
double A11_32(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+mu*dt*(2*g2x*g3x+g2y*g3y)+2*dt*(mu+lambda)*(g2x*g3x);
}	
inline
double A11_33(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+mu*dt*(2*g3x*g3x+g3y*g3y)+2*dt*(mu+lambda)*(g3x*g3x);
}	
inline
double A11_34(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+mu*dt*(2*g4x*g3x+g4y*g3y)+2*dt*(mu+lambda)*(g4x*g3x);
}	
inline
double A11_35(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+mu*dt*(2*g5x*g3x+g5y*g3y)+2*dt*(mu+lambda)*(g5x*g3x);
}

//#######################
inline
double A11_40(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+mu*dt*(2*g0x*g4x+g0y*g4y)+2*dt*(mu+lambda)*(g0x*g4x);
}	
inline
double A11_41(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+mu*dt*(2*g1x*g4x+g1y*g4y)+2*dt*(mu+lambda)*(g1x*g4x);
}	
inline
double A11_42(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+mu*dt*(2*g2x*g4x+g2y*g4y)+2*dt*(mu+lambda)*(g2x*g4x);
}	
inline
double A11_43(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+mu*dt*(2*g3x*g4x+g3y*g4y)+2*dt*(mu+lambda)*(g3x*g4x);
}	
inline
double A11_44(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+mu*dt*(2*g4x*g4x+g4y*g4y)+2*dt*(mu+lambda)*(g4x*g4x);
}	
inline
double A11_45(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+mu*dt*(2*g5x*g4x+g5y*g4y)+2*dt*(mu+lambda)*(g5x*g4x);
}

//#######################
inline
double A11_50(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+mu*dt*(2*g0x*g5x+g0y*g5y)+2*dt*(mu+lambda)*(g0x*g5x);
}	
inline
double A11_51(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+mu*dt*(2*g1x*g5x+g1y*g5y)+2*dt*(mu+lambda)*(g1x*g5x);
}	
inline
double A11_52(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+mu*dt*(2*g2x*g5x+g2y*g5y)+2*dt*(mu+lambda)*(g2x*g5x);
}	
inline
double A11_53(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+mu*dt*(2*g3x*g5x+g3y*g5y)+2*dt*(mu+lambda)*(g3x*g5x);
}	
inline
double A11_54(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+mu*dt*(2*g4x*g5x+g4y*g5y)+2*dt*(mu+lambda)*(g4x*g5x);
}
inline
double A11_55(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+mu*dt*(2*g5x*g5x+g5y*g5y)+2*dt*(mu+lambda)*(g5x*g5x);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A12_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g0+mu*dt*g0y*g0x+2*dt*(mu+lambda)*(g0y*g0x);
}	
inline
double A12_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g1+mu*dt*g1y*g0x+2*dt*(mu+lambda)*(g1y*g0x);
}	
inline
double A12_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g2+mu*dt*g2y*g0x+2*dt*(mu+lambda)*(g2y*g0x);
}	
inline
double A12_03(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g3+mu*dt*g3y*g0x+2*dt*(mu+lambda)*(g3y*g0x);
}	
inline
double A12_04(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g4+mu*dt*g4y*g0x+2*dt*(mu+lambda)*(g4y*g0x);
}	
inline
double A12_05(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g5+mu*dt*g5y*g0x+2*dt*(mu+lambda)*(g5y*g0x);
}
//#####################
inline
double A12_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g0+mu*dt*g0y*g1x+2*dt*(mu+lambda)*(g0y*g1x);
}	
inline
double A12_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g1+mu*dt*g1y*g1x+2*dt*(mu+lambda)*(g1y*g1x);
}	
inline
double A12_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g2+mu*dt*g2y*g1x+2*dt*(mu+lambda)*(g2y*g1x);
}	
inline
double A12_13(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g3+mu*dt*g3y*g1x+2*dt*(mu+lambda)*(g3y*g1x);
}	
inline
double A12_14(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g4+mu*dt*g4y*g1x+2*dt*(mu+lambda)*(g4y*g1x);
}	
inline
double A12_15(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g5+mu*dt*g5y*g1x+2*dt*(mu+lambda)*(g5y*g1x);
}
//#####################
inline
double A12_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g0+mu*dt*g0y*g2x+2*dt*(mu+lambda)*(g0y*g2x);
}	
inline
double A12_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g1+mu*dt*g1y*g2x+2*dt*(mu+lambda)*(g1y*g2x);
}	
inline
double A12_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g2+mu*dt*g2y*g2x+2*dt*(mu+lambda)*(g2y*g2x);
}	
inline
double A12_23(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g3+mu*dt*g3y*g2x+2*dt*(mu+lambda)*(g3y*g2x);
}	
inline
double A12_24(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g4+mu*dt*g4y*g2x+2*dt*(mu+lambda)*(g4y*g2x);
}	
inline
double A12_25(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g5+mu*dt*g5y*g2x+2*dt*(mu+lambda)*(g5y*g2x);
}
//#######################
inline
double A12_30(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g0+mu*dt*g0y*g3x+2*dt*(mu+lambda)*(g0y*g3x);
}	
inline
double A12_31(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g1+mu*dt*g1y*g3x+2*dt*(mu+lambda)*(g1y*g3x);
}	
inline
double A12_32(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g2+mu*dt*g2y*g3x+2*dt*(mu+lambda)*(g2y*g3x);
}	
inline
double A12_33(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g3+mu*dt*g3y*g3x+2*dt*(mu+lambda)*(g3y*g3x);
}	
inline
double A12_34(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g4+mu*dt*g4y*g3x+2*dt*(mu+lambda)*(g4y*g3x);
}	
inline
double A12_35(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g5+mu*dt*g5y*g3x+2*dt*(mu+lambda)*(g5y*g3x);
}
//###################################
inline
double A12_40(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g0+mu*dt*g0y*g4x+2*dt*(mu+lambda)*(g0y*g4x);
}	
inline
double A12_41(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g1+mu*dt*g1y*g4x+2*dt*(mu+lambda)*(g1y*g4x);
}	
inline
double A12_42(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g2+mu*dt*g2y*g4x+2*dt*(mu+lambda)*(g2y*g4x);
}	
inline
double A12_43(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g3+mu*dt*g3y*g4x+2*dt*(mu+lambda)*(g3y*g4x);
}	
inline
double A12_44(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g4+mu*dt*g4y*g4x+2*dt*(mu+lambda)*(g4y*g4x);
}	
inline
double A12_45(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g5+mu*dt*g5y*g4x+2*dt*(mu+lambda)*(g5y*g4x);
}
//#########################
inline
double A12_50(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g0+mu*dt*g0y*g5x+2*dt*(mu+lambda)*(g0y*g5x);
}	
inline
double A12_51(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g1+mu*dt*g1y*g5x+2*dt*(mu+lambda)*(g1y*g5x);
}	
inline
double A12_52(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g2+mu*dt*g2y*g5x+2*dt*(mu+lambda)*(g2y*g5x);
}	
inline
double A12_53(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g3+mu*dt*g3y*g5x+2*dt*(mu+lambda)*(g3y*g5x);
}	
inline
double A12_54(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g4+mu*dt*g4y*g5x+2*dt*(mu+lambda)*(g4y*g5x);
}	
inline
double A12_55(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g5+mu*dt*g5y*g5x+2*dt*(mu+lambda)*(g5y*g5x);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A20_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f0
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g0*f0y;
}	 
inline
double A20_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f1
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g0*f1x;
}	
inline
double A20_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f2
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g0*f2x;
}	
inline
double A20_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f0
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g1*f0x;
}	
inline
double A20_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f1
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g1*f1x;
}	
inline
double A20_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f2
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g1*f2x;
}	
inline
double A20_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f0
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g2*f0x;
}	
inline
double A20_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f1
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g2*f1x;
}	
inline
double A20_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f2
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g2*f2x;
}	
inline
double A20_30(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f0
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g3*f0x;
}	
inline
double A20_31(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f1
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g3*f1x;
}	
inline
double A20_32(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f2
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g3*f2x;
}	
inline
double A20_40(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f0
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g4*f0x;
}	
inline
double A20_41(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f1
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g4*f1x;
}	
inline
double A20_42(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f2
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g4*f2x;
}	
inline
double A20_50(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f0
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g5*f0x;
}	
inline
double A20_51(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f1
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g5*f1x;
}	
inline
double A20_52(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m,  double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*((v0_m-v0_n)*g0+(v1_m-v1_n)*g1+(v2_m-v2_n)*g2+(v3_m-v3_n)*g3+(v4_m-v4_n)*g4+(v5_m-v5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)
	          +dt*a*gamma*(gamma-1)*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-2)*(r0_m*f0y+r1_m*f1y+r2_m*f2y))*f2
	          +dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*g5*f2x;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A21_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g0+mu*dt*g0y*g0x+2*dt*(mu+lambda)*(g0y*g0x);
}	
inline
double A21_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g1+mu*dt*g1y*g0x+2*dt*(mu+lambda)*(g1y*g0x);
}	
inline
double A21_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g2+mu*dt*g2y*g0x+2*dt*(mu+lambda)*(g2y*g0x);
}	
inline
double A21_03(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g3+mu*dt*g3y*g0x+2*dt*(mu+lambda)*(g3y*g0x);
}	
inline
double A21_04(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g4+mu*dt*g4y*g0x+2*dt*(mu+lambda)*(g4y*g0x);
}	
inline
double A21_05(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g5+mu*dt*g5y*g0x+2*dt*(mu+lambda)*(g5y*g0x);
}
//#####################
inline
double A21_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g0+mu*dt*g0y*g1x+2*dt*(mu+lambda)*(g0y*g1x);
}	
inline
double A21_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g1+mu*dt*g1y*g1x+2*dt*(mu+lambda)*(g1y*g1x);
}	
inline
double A21_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g2+mu*dt*g2y*g1x+2*dt*(mu+lambda)*(g2y*g1x);
}	
inline
double A21_13(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g3+mu*dt*g3y*g1x+2*dt*(mu+lambda)*(g3y*g1x);
}	
inline
double A21_14(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g4+mu*dt*g4y*g1x+2*dt*(mu+lambda)*(g4y*g1x);
}	
inline
double A21_15(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g5+mu*dt*g5y*g1x+2*dt*(mu+lambda)*(g5y*g1x);
}
//#####################
inline
double A21_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g0+mu*dt*g0y*g2x+2*dt*(mu+lambda)*(g0y*g2x);
}	
inline
double A21_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g1+mu*dt*g1y*g2x+2*dt*(mu+lambda)*(g1y*g2x);
}	
inline
double A21_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g2+mu*dt*g2y*g2x+2*dt*(mu+lambda)*(g2y*g2x);
}	
inline
double A21_23(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g3+mu*dt*g3y*g2x+2*dt*(mu+lambda)*(g3y*g2x);
}	
inline
double A21_24(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g4+mu*dt*g4y*g2x+2*dt*(mu+lambda)*(g4y*g2x);
}	
inline
double A21_25(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g5+mu*dt*g5y*g2x+2*dt*(mu+lambda)*(g5y*g2x);
}
//#######################
inline
double A21_30(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g0+mu*dt*g0y*g3x+2*dt*(mu+lambda)*(g0y*g3x);
}	
inline
double A21_31(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g1+mu*dt*g1y*g3x+2*dt*(mu+lambda)*(g1y*g3x);
}	
inline
double A21_32(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g2+mu*dt*g2y*g3x+2*dt*(mu+lambda)*(g2y*g3x);
}	
inline
double A21_33(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g3+mu*dt*g3y*g3x+2*dt*(mu+lambda)*(g3y*g3x);
}	
inline
double A21_34(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g4+mu*dt*g4y*g3x+2*dt*(mu+lambda)*(g4y*g3x);
}	
inline
double A21_35(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g5+mu*dt*g5y*g3x+2*dt*(mu+lambda)*(g5y*g3x);
}
//###################################
inline
double A21_40(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g0+mu*dt*g0y*g4x+2*dt*(mu+lambda)*(g0y*g4x);
}	
inline
double A21_41(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g1+mu*dt*g1y*g4x+2*dt*(mu+lambda)*(g1y*g4x);
}	
inline
double A21_42(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g2+mu*dt*g2y*g4x+2*dt*(mu+lambda)*(g2y*g4x);
}	
inline
double A21_43(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g3+mu*dt*g3y*g4x+2*dt*(mu+lambda)*(g3y*g4x);
}	
inline
double A21_44(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g4+mu*dt*g4y*g4x+2*dt*(mu+lambda)*(g4y*g4x);
}	
inline
double A21_45(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g5+mu*dt*g5y*g4x+2*dt*(mu+lambda)*(g5y*g4x);
}
//#########################
inline
double A21_50(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g0+mu*dt*g0y*g5x+2*dt*(mu+lambda)*(g0y*g5x);
}	
inline
double A21_51(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g1+mu*dt*g1y*g5x+2*dt*(mu+lambda)*(g1y*g5x);
}	
inline
double A21_52(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g2+mu*dt*g2y*g5x+2*dt*(mu+lambda)*(g2y*g5x);
}	
inline
double A21_53(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g3+mu*dt*g3y*g5x+2*dt*(mu+lambda)*(g3y*g5x);
}	
inline
double A21_54(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g4+mu*dt*g4y*g5x+2*dt*(mu+lambda)*(g4y*g5x);
}	
inline
double A21_55(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g5+mu*dt*g5y*g5x+2*dt*(mu+lambda)*(g5y*g5x);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A22_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+mu*dt*(g0x*g0x+2*g0y*g0y)+2*dt*(mu+lambda)*(g0y*g0y);
}	
inline
double A22_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+mu*dt*(g1x*g0x+2*g1y*g0y)+2*dt*(mu+lambda)*(g1y*g0y);
}	
inline
double A22_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+mu*dt*(g2x*g0x+2*g2y*g0y)+2*dt*(mu+lambda)*(g2y*g0y);
}	
inline
double A22_03(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+mu*dt*(g3x*g0x+2*g3y*g0y)+2*dt*(mu+lambda)*(g3y*g0y);
}	
inline
double A22_04(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+mu*dt*(g4x*g0x+2*g4y*g0y)+2*dt*(mu+lambda)*(g4y*g0y);
}	
inline
double A22_05(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+mu*dt*(g5x*g0x+2*g5y*g0y)+2*dt*(mu+lambda)*(g5y*g0y);
}
//#############################
inline
double A22_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+mu*dt*(g0x*g1x+2*g0y*g1y)+2*dt*(mu+lambda)*(g0y*g1y);
}	
inline
double A22_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+mu*dt*(g1x*g1x+2*g1y*g1y)+2*dt*(mu+lambda)*(g1y*g1y);
}	
inline
double A22_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+mu*dt*(g2x*g1x+2*g2y*g1y)+2*dt*(mu+lambda)*(g2y*g1y);
}	
inline
double A22_13(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+mu*dt*(g3x*g1x+2*g3y*g1y)+2*dt*(mu+lambda)*(g3y*g1y);
}	
inline
double A22_14(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+mu*dt*(g4x*g1x+2*g4y*g1y)+2*dt*(mu+lambda)*(g4y*g1y);
}	
inline
double A22_15(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+mu*dt*(g5x*g1x+2*g5y*g1y)+2*dt*(mu+lambda)*(g5y*g1y);
}
//#######################
inline
double A22_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+mu*dt*(g0x*g2x+2*g0y*g2y)+2*dt*(mu+lambda)*(g0y*g2y);
}	
inline
double A22_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+mu*dt*(g1x*g2x+2*g1y*g2y)+2*dt*(mu+lambda)*(g1y*g2y);
}	
inline
double A22_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+mu*dt*(g2x*g2x+2*g2y*g2y)+2*dt*(mu+lambda)*(g2y*g2y);
}	
inline
double A22_23(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+mu*dt*(g3x*g2x+2*g3y*g2y)+2*dt*(mu+lambda)*(g3y*g2y);
}	
inline
double A22_24(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+mu*dt*(g4x*g2x+2*g4y*g2y)+2*dt*(mu+lambda)*(g4y*g2y);
}	
inline
double A22_25(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+mu*dt*(g5x*g2x+2*g5y*g2y)+2*dt*(mu+lambda)*(g5y*g2y);
}

//#######################
inline
double A22_30(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+mu*dt*(g0x*g3x+2*g0y*g3y)+2*dt*(mu+lambda)*(g0y*g3y);
}	
inline
double A22_31(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+mu*dt*(g1x*g3x+2*g1y*g3y)+2*dt*(mu+lambda)*(g1y*g3y);
}	
inline
double A22_32(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+mu*dt*(g2x*g3x+2*g2y*g3y)+2*dt*(mu+lambda)*(g2y*g3y);
}	
inline
double A22_33(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+mu*dt*(g3x*g3x+2*g3y*g3y)+2*dt*(mu+lambda)*(g3y*g3y);
}	
inline
double A22_34(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+mu*dt*(g4x*g3x+2*g4y*g3y)+2*dt*(mu+lambda)*(g4y*g3y);
}	
inline
double A22_35(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+mu*dt*(g5x*g3x+2*g5y*g3y)+2*dt*(mu+lambda)*(g5y*g3y);
}

//#######################
inline
double A22_40(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+mu*dt*(g0x*g4x+2*g0y*g4y)+2*dt*(mu+lambda)*(g0y*g4y);
}	
inline
double A22_41(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+mu*dt*(g1x*g4x+2*g1y*g4y)+2*dt*(mu+lambda)*(g1y*g4y);
}	
inline
double A22_42(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+mu*dt*(g2x*g4x+2*g2y*g4y)+2*dt*(mu+lambda)*(g2y*g4y);
}	
inline
double A22_43(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+mu*dt*(g3x*g4x+2*g3y*g4y)+2*dt*(mu+lambda)*(g3y*g4y);
}	
inline
double A22_44(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+mu*dt*(g4x*g4x+2*g4y*g4y)+2*dt*(mu+lambda)*(g4y*g4y);
}	
inline
double A22_45(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+mu*dt*(g5x*g4x+2*g5y*g4y)+2*dt*(mu+lambda)*(g5y*g4y);
}

//#######################
inline
double A22_50(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+mu*dt*(g0x*g5x+2*g0y*g5y)+2*dt*(mu+lambda)*(g0y*g5y);
}	
inline
double A22_51(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+mu*dt*(g1x*g5x+2*g1y*g5y)+2*dt*(mu+lambda)*(g1y*g5y);
}	
inline
double A22_52(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+mu*dt*(g2x*g5x+2*g2y*g5y)+2*dt*(mu+lambda)*(g2y*g5y);
}	
inline
double A22_53(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+mu*dt*(g3x*g5x+2*g3y*g5y)+2*dt*(mu+lambda)*(g3y*g5y);
}	
inline
double A22_54(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+mu*dt*(g4x*g5x+2*g4y*g5y)+2*dt*(mu+lambda)*(g4y*g5y);
}	
inline
double A22_55(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+mu*dt*(g5x*g5x+2*g5y*g5y)+2*dt*(mu+lambda)*(g5y*g5y);
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A30_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f0;
}	 
inline
double A30_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f1;
}	
inline
double A30_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f2;
}	
inline
double A30_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f0;
}	
inline
double A30_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f1;
}	
inline
double A30_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f2;
}	
inline
double A30_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f0;
}	
inline
double A30_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
return g2*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f1;
}	
inline
double A30_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f2;
}	
inline
double A30_30(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f0;
}	
inline
double A30_31(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f1;
}	
inline
double A30_32(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f2;
}	
inline
double A30_40(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f0;
}	
inline
double A30_41(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f1;
}	
inline
double A30_42(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f2;
}	
inline
double A30_50(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f0;
}	
inline
double A30_51(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f1;
}	
inline
double A30_52(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*((w0_m-w0_n)*g0+(w1_m-w1_n)*g1+(w2_m-w2_n)*g2+(w3_m-w3_n)*g3+(w4_m-w4_n)*g4+(w5_m-w5_n)*g5
	          +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	          +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*f2;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A31_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g0;
}	
inline
double A31_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g1;
}	
inline
double A31_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g2;
}	
inline
double A31_03(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g3;
}	
inline
double A31_04(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g4;
}	
inline
double A31_05(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g5;
}
//#####################
inline
double A31_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g0;
}	
inline
double A31_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g1;
}	
inline
double A31_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g2;
}	
inline
double A31_13(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g3;
}	
inline
double A31_14(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g4;
}	
inline
double A31_15(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g5;
}
//#####################
inline
double A31_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g0;
}	
inline
double A31_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g1;
}	
inline
double A31_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g2;
}	
inline
double A31_23(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g3;
}	
inline
double A31_24(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g4;
}	
inline
double A31_25(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g5;
}
//#######################
inline
double A31_30(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g0;
}	
inline
double A31_31(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g1;
}	
inline
double A31_32(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g2;
}	
inline
double A31_33(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g3;
}	
inline
double A31_34(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g4;
}	
inline
double A31_35(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g5;
}
//###################################
inline
double A31_40(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g0;
}	
inline
double A31_41(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g1;
}	
inline
double A31_42(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g2;
}	
inline
double A31_43(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g3;
}	
inline
double A31_44(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g4;
}	
inline
double A31_45(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g5;
}
//#########################
inline
double A31_50(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g0;
}	
inline
double A31_51(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g1;
}	
inline
double A31_52(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g2;
}	
inline
double A31_53(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g3;
}	
inline
double A31_54(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g4;
}	
inline
double A31_55(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g5;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A32_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g0;
}	
inline
double A32_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g1;
}	
inline
double A32_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g2;
}	
inline
double A32_03(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g3;
}	
inline
double A32_04(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g4;
}	
inline
double A32_05(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g5;
}
//#############################
inline
double A32_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g0;
}	
inline
double A32_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g1;
}	
inline
double A32_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g2;
}	
inline
double A32_13(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g3;
}	
inline
double A32_14(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g4;
}	
inline
double A32_15(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g5;
}
//#######################
inline
double A32_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g0;
}	
inline
double A32_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g1;
}	
inline
double A32_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g2;
}	
inline
double A32_23(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g3;
}	
inline
double A32_24(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g4;
}	
inline
double A32_25(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g5;
}

//#######################
inline
double A32_30(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g0;
}	
inline
double A32_31(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g1;
}	
inline
double A32_32(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g2;
}	
inline
double A32_33(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g3;
}	
inline
double A32_34(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g4;
}	
inline
double A32_35(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g5;
}

//#######################
inline
double A32_40(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g0;
}	
inline
double A32_41(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g1;
}	
inline
double A32_42(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g2;
}	
inline
double A32_43(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g3;
}	
inline
double A32_44(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g4;
}	
inline
double A32_45(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g5;
}

//#######################
inline
double A32_50(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g0;
}	
inline
double A32_51(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g1;
}	
inline
double A32_52(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g2;
}	
inline
double A32_53(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g3;
}	
inline
double A32_54(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g4;
}	
inline
double A32_55(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g5;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A33_00(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x+dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+(kt/cp)*dt*(g0x*g0x+g0y*g0y);
}	
inline
double A33_01(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+(kt/cp)*dt*(g1x*g0x+g1y*g0y);
}	
inline
double A33_02(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+(kt/cp)*dt*(g2x*g0x+g2y*g0y);
}	
inline
double A33_03(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+(kt/cp)*dt*(g3x*g0x+g3y*g0y);
}	
inline
double A33_04(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+(kt/cp)*dt*(g4x*g0x+g4y*g0y);
}	
inline
double A33_05(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g0*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+(kt/cp)*dt*(g5x*g0x+g5y*g0y);
}
//#############################
inline
double A33_10(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+(kt/cp)*dt*(g0x*g1x+g0y*g1y);
}	
inline
double A33_11(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+(kt/cp)*dt*(g1x*g1x+g1y*g1y);
}	
inline
double A33_12(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+(kt/cp)*dt*(g2x*g1x+g2y*g1y);
}	
inline
double A33_13(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+(kt/cp)*dt*(g3x*g1x+g3y*g1y);
}	
inline
double A33_14(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+(kt/cp)*dt*(g4x*g1x+g4y*g1y);
}	
inline
double A33_15(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g1*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+(kt/cp)*dt*(g5x*g1x+g5y*g1y);
}
//#######################
inline
double A33_20(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+(kt/cp)*dt*(g0x*g2x+g0y*g2y);
}	
inline
double A33_21(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+(kt/cp)*dt*(g1x*g2x+g1y*g2y);
}	
inline
double A33_22(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+(kt/cp)*dt*(g2x*g2x+g2y*g2y);
}	
inline
double A33_23(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+(kt/cp)*dt*(g3x*g2x+g3y*g2y);
}	
inline
double A33_24(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+(kt/cp)*dt*(g4x*g2x+g4y*g2y);
}	
inline
double A33_25(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g2*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+(kt/cp)*dt*(g5x*g2x+g5y*g2y);
}

//#######################
inline
double A33_30(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+(kt/cp)*dt*(g0x*g3x+g0y*g3y);
}	
inline
double A33_31(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+(kt/cp)*dt*(g1x*g3x+g1y*g3y);
}	
inline
double A33_32(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+(kt/cp)*dt*(g2x*g3x+g2y*g3y);
}	
inline
double A33_33(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+(kt/cp)*dt*(g3x*g3x+g3y*g3y);
}	
inline
double A33_34(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+(kt/cp)*dt*(g4x*g3x+g4y*g3y);
}	
inline
double A33_35(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g3*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+(kt/cp)*dt*(g5x*g3x+g5y*g3y);
}

//#######################
inline
double A33_40(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+(kt/cp)*dt*(g0x*g4x+g0y*g4y);
}	
inline
double A33_41(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+(kt/cp)*dt*(g1x*g4x+g1y*g4y);
}	
inline
double A33_42(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+(kt/cp)*dt*(g2x*g4x+g2y*g4y);
}	
inline
double A33_43(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+(kt/cp)*dt*(g3x*g4x+g3y*g4y);
}	
inline
double A33_44(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+(kt/cp)*dt*(g4x*g4x+g4y*g4y);
}	
inline
double A33_45(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g4*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+(kt/cp)*dt*(g5x*g4x+g5y*g4y);
}

//#######################
inline
double A33_50(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g0+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g0x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g0y)+(kt/cp)*dt*(g0x*g5x+g0y*g5y);
}	
inline
double A33_51(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g1+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g1x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g1y)+(kt/cp)*dt*(g1x*g5x+g1y*g5y);
}	
inline
double A33_52(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g2+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g2x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g2y)+(kt/cp)*dt*(g2x*g5x+g2y*g5y);
}	
inline
double A33_53(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g3+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g3x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g3y)+(kt/cp)*dt*(g3x*g5x+g3y*g5y);
}	
inline
double A33_54(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g4+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g4x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g4y)+(kt/cp)*dt*(g4x*g5x+g4y*g5y);
}	
inline
double A33_55(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return g5*(r0_m*f0+r1_m*f1+r2_m*f2)*(g5+dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*g5x
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*g5y)+(kt/cp)*dt*(g5x*g5x+g5y*g5y);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double B0_0(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_n, double r1_n, double r2_n, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return ((r0_n-r0_m)*f0+(r1_n-r1_m)*f1+(r2_n-r2_m)*f2-dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(r0_m*f0x+r1_m*f1x+r2_m*f2x)-dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(r0_m*f0y+r1_m*f1y+r2_m*f2y)
	       -dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)-dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*f0;
}	
inline
double B0_1(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_n, double r1_n, double r2_n, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return ((r0_n-r0_m)*f0+(r1_n-r1_m)*f1+(r2_n-r2_m)*f2-dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(r0_m*f0x+r1_m*f1x+r2_m*f2x)-dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(r0_m*f0y+r1_m*f1y+r2_m*f2y)
	       -dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)-dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*f1;
}	
inline
double B0_2(double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_n, double r1_n, double r2_n, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return ((r0_n-r0_m)*f0+(r1_n-r1_m)*f1+(r2_n-r2_m)*f2-dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(r0_m*f0x+r1_m*f1x+r2_m*f2x)-dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(r0_m*f0y+r1_m*f1y+r2_m*f2y)
	       -dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)-dt*(r0_m*f0+r1_m*f1+r2_m*f2)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*f2;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double B1_0(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5-(u0_n*g0+u1_n*g1+u2_n*g2+u3_n*g3+u4_n*g4+u5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y))*g0
	       +mu*dt*(2*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g0x
	       +(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g0y
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g0x)
	       +2*dt*(mu+lambda)*((u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g0x
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g0x)
	       -dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*(r0_m*f0x+r1_m*f1x+r2_m*f2x)*g0);//-(M_PI*cos(M_PI*x)*sin(M_PI*x)*exp(-4*M_PI*M_PI*t_ni))*g0;
}	 
inline
double B1_1(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5-(u0_n*g0+u1_n*g1+u2_n*g2+u3_n*g3+u4_n*g4+u5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y))*g1
	       +mu*dt*(2*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g1x
	       +(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g1y
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g1x)
	       +2*dt*(mu+lambda)*((u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g1x
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g1x)
	       -dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*(r0_m*f0x+r1_m*f1x+r2_m*f2x)*g1);//-(M_PI*cos(M_PI*x)*sin(M_PI*x)*exp(-4*M_PI*M_PI*t_ni))*g1;
}	
inline
double B1_2(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5-(u0_n*g0+u1_n*g1+u2_n*g2+u3_n*g3+u4_n*g4+u5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y))*g2
	       +mu*dt*(2*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g2x
	       +(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g2y
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g2x)
	       +2*dt*(mu+lambda)*((u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g2x
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g2x)
	       -dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*(r0_m*f0x+r1_m*f1x+r2_m*f2x)*g2);//-(M_PI*cos(M_PI*x)*sin(M_PI*x)*exp(-4*M_PI*M_PI*t_ni))*g2;
}	  
inline
double B1_3(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5-(u0_n*g0+u1_n*g1+u2_n*g2+u3_n*g3+u4_n*g4+u5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y))*g3
	       +mu*dt*(2*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g3x
	       +(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g3y
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g3x)
	       +2*dt*(mu+lambda)*((u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g3x
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g3x)
	       -dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*(r0_m*f0x+r1_m*f1x+r2_m*f2x)*g3);//-(M_PI*cos(M_PI*x)*sin(M_PI*x)*exp(-4*M_PI*M_PI*t_ni))*g3;
}	 
inline
double B1_4(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5-(u0_n*g0+u1_n*g1+u2_n*g2+u3_n*g3+u4_n*g4+u5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y))*g4
	       +mu*dt*(2*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g4x
	       +(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g4y
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g4x)
	       +2*dt*(mu+lambda)*((u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g4x
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g4x)
	       -dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*(r0_m*f0x+r1_m*f1x+r2_m*f2x)*g4);//-(M_PI*cos(M_PI*x)*sin(M_PI*x)*exp(-4*M_PI*M_PI*t_ni))*g4;
}	 
inline
double B1_5(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5-(u0_n*g0+u1_n*g1+u2_n*g2+u3_n*g3+u4_n*g4+u5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y))*g5
	       +mu*dt*(2*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g5x
	       +(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g5y
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g5x)
	       +2*dt*(mu+lambda)*((u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x)*g5x
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g5x)
	       -dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*(r0_m*f0x+r1_m*f1x+r2_m*f2x)*g5);//-(M_PI*cos(M_PI*x)*sin(M_PI*x)*exp(-4*M_PI*M_PI*t_ni))*g5;
}	 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double B2_0(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5-(v0_n*g0+v1_n*g1+v2_n*g2+v3_n*g3+v4_n*g4+v5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*g0
	       +mu*dt*((v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g0x
	       +2*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g0y
	       +(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g0x)
	       +2*dt*(mu+lambda)*((u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g0x
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g0y)
	       -dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*(r0_m*f0y+r1_m*f1y+r2_m*f2y)*g0);//-(M_PI*cos(M_PI*y)*sin(M_PI*y)*exp(-4*M_PI*M_PI*t_ni))*g0;
}	 
inline
double B2_1(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5-(v0_n*g0+v1_n*g1+v2_n*g2+v3_n*g3+v4_n*g4+v5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*g1
	       +mu*dt*((v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g1x
	       +2*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g1y
	       +(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g1x)
	       +2*dt*(mu+lambda)*((u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g1x
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g1y)
	       -dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*(r0_m*f0y+r1_m*f1y+r2_m*f2y)*g1);//-(M_PI*cos(M_PI*y)*sin(M_PI*y)*exp(-4*M_PI*M_PI*t_ni))*g1;
}	 
inline
double B2_2(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5-(v0_n*g0+v1_n*g1+v2_n*g2+v3_n*g3+v4_n*g4+v5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*g2
	       +mu*dt*((v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g2x
	       +2*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g2y
	       +(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g2x)
	       +2*dt*(mu+lambda)*((u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g2x
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g2y)
	       -dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*(r0_m*f0y+r1_m*f1y+r2_m*f2y)*g2);//-(M_PI*cos(M_PI*y)*sin(M_PI*y)*exp(-4*M_PI*M_PI*t_ni))*g2;
}	 
inline
double B2_3(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5-(v0_n*g0+v1_n*g1+v2_n*g2+v3_n*g3+v4_n*g4+v5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*g3
	       +mu*dt*((v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g3x
	       +2*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g3y
	       +(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g3x)
	       +2*dt*(mu+lambda)*((u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g3x
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g3y)
	       -dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*(r0_m*f0y+r1_m*f1y+r2_m*f2y)*g3);//-(M_PI*cos(M_PI*y)*sin(M_PI*y)*exp(-4*M_PI*M_PI*t_ni))*g3;
}	 
inline
double B2_4(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5-(v0_n*g0+v1_n*g1+v2_n*g2+v3_n*g3+v4_n*g4+v5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*g4
	       +mu*dt*((v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g4x
	       +2*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g4y
	       +(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g4x)
	       +2*dt*(mu+lambda)*((u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g4x
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g4y)
	       -dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*(r0_m*f0y+r1_m*f1y+r2_m*f2y)*g4);//-(M_PI*cos(M_PI*y)*sin(M_PI*y)*exp(-4*M_PI*M_PI*t_ni))*g4;
}	 
inline
double B2_5(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5-(v0_n*g0+v1_n*g1+v2_n*g2+v3_n*g3+v4_n*g4+v5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y))*g5
	       +mu*dt*((v0_m*g0x+v1_m*g1x+v2_m*g2x+v3_m*g3x+v4_m*g4x+v5_m*g5x)*g5x
	       +2*(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g5y
	       +(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g5x)
	       +2*dt*(mu+lambda)*((u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y)*g5x
	       +(v0_m*g0y+v1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y)*g5y)
	       -dt*a*gamma*pow(r0_m*f0+r1_m*f1+r2_m*f2,gamma-1)*(r0_m*f0y+r1_m*f1y+r2_m*f2y)*g5);//-(M_PI*cos(M_PI*y)*sin(M_PI*y)*exp(-4*M_PI*M_PI*t_ni))*g5;
}	 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double B3_0(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5-(w0_n*g0+w1_n*g1+w2_n*g2+w3_n*g3+w4_n*g4+w5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*g0
	       +(kt/cp)*dt*((w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g0x
	              +(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g0y));
}	 
inline
double B3_1(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5-(w0_n*g0+w1_n*g1+w2_n*g2+w3_n*g3+w4_n*g4+w5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*g1
	       +(kt/cp)*dt*((w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g1x
	              +(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g1y));
}	 
inline
double B3_2(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5-(w0_n*g0+w1_n*g1+w2_n*g2+w3_n*g3+w4_n*g4+w5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*g2
	       +(kt/cp)*dt*((w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g2x
	              +(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g2y));
}	
inline
double B3_3(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5-(w0_n*g0+w1_n*g1+w2_n*g2+w3_n*g3+w4_n*g4+w5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*g3
	       +(kt/cp)*dt*((w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g3x
	              +(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g3y));
}	
inline
double B3_4(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
   return -((r0_m*f0+r1_m*f1+r2_m*f2)*((w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5-(w0_n*g0+w1_n*g1+w2_n*g2+w3_n*g3+w4_n*g4+w5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*g4
	       +(kt/cp)*dt*((w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g4x
	              +(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g4y));
}	
inline
double B3_5(double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m)
{ 
     double V = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
     
     double f0 = (1.0/V)*(a2*b3 - a3*b2 - a2*y + b2*x + a3*y - b3*x);
     double f1 = (1.0/V)*(a3*b1 - a1*b3 + a1*y - b1*x - a3*y + b3*x); 
     double f2 = (1.0/V)*(a1*b2 - a2*b1 - a1*y + b1*x + a2*y - b2*x); 
     
     
     double f0x = (1.0/V)*(b2 - b3); double f0y = (1.0/V)*(a3 - a2); 
     double f1x = (1.0/V)*(b3 - b1); double f1y = (1.0/V)*(a1 - a3); 
     double f2x = (1.0/V)*(b1 - b2); double f2y = (1.0/V)*(a2 - a1); 
     
     
     double g0 = f0*(2*f0-1); double g1 = f1*(2*f1-1); double g2 = f2*(2*f2-1); double g3 = 4*f0*f1; double g4 = 4*f1*f2; double g5 = 4*f0*f2; 
     
     double g0x = f0x*(2*f0-1)+f0*(2*f0x); double g1x = f1x*(2*f1-1)+f1*(2*f1x); double g2x = f2x*(2*f2-1)+f2*(2*f2x); double g3x = 4*f0x*f1+4*f0*f1x; double g4x = 4*f1x*f1+4*f1*f2x; double g5x = 4*f0x*f2+4*f0*f2x; 
     double g0y = f0y*(2*f0-1)+f0*(2*f0y); double g1y = f1y*(2*f1-1)+f1*(2*f1y); double g2y = f2y*(2*f2-1)+f2*(2*f2y); double g3y = 4*f0y*f1+4*f0*f1y; double g4y = 4*f1y*f1+4*f1*f2y; double g5y = 4*f0y*f2+4*f0*f2y; 
     
     
    	       
    return -((r0_m*f0+r1_m*f1+r2_m*f2)*((w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5-(w0_n*g0+w1_n*g1+w2_n*g2+w3_n*g3+w4_n*g4+w5_n*g5))
	       +dt*(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5)*(w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)
	       +dt*(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5)*(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y))*g5
	       +(kt/cp)*dt*((w0_m*g0x+w1_m*g1x+w2_m*g2x+w3_m*g3x+w4_m*g4x+w5_m*g5x)*g5x
	              +(w0_m*g0y+w1_m*g1y+w2_m*g2y+w3_m*g3y+w4_m*g4y+w5_m*g5y)*g5y));
}	



