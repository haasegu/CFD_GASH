#ifndef GETMATRIX_FILE
#define GETMATRIX_FILE

#include "geom3.h"
#include <cassert>
#include <string>
#include <iostream>
#include <vector>

/**
 * Calculates the element stiffness matrix @p ske and the element load vector @p fe
 * of one triangular element with linear shape functions.
 * @param[in]	ial	node indices of the three element vertices
 * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coodinates of node k
 * @param[out] ske	element stiffness matrix
 * @param[out] fe	element load vector
 */
 
 
void CalcElem(int const ial[4], double const xc[], double ske[4][4], double fe[4]);

//first row
void CalcElem_Navier_Stokes_A00(int const ial[4], double const xc[], double ske[4][4], double fe[4], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

void CalcElem_Navier_Stokes_A01(int const ial[10], double const xc[], double ske[4][10], double fe[4], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

void CalcElem_Navier_Stokes_A02(int const ial[10], double const xc[], double ske[4][10], double fe[4], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

void CalcElem_Navier_Stokes_A03(int const ial[10], double const xc[], double ske[4][10], double fe[4], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

//second row
void CalcElem_Navier_Stokes_A10(int const ial[10], double const xc[], double ske[10][4], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

void CalcElem_Navier_Stokes_A11(int const ial[10], double const xc[], double ske[10][4], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

void CalcElem_Navier_Stokes_A12(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

void CalcElem_Navier_Stokes_A13(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

//third row
void CalcElem_Navier_Stokes_A20(int const ial[10], double const xc[], double ske[10][4], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

void CalcElem_Navier_Stokes_A21(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

void CalcElem_Navier_Stokes_A22(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

void CalcElem_Navier_Stokes_A23(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

//fourth row
void CalcElem_Navier_Stokes_A30(int const ial[10], double const xc[], double ske[10][4], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

void CalcElem_Navier_Stokes_A31(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

void CalcElem_Navier_Stokes_A31(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);

void CalcElem_Navier_Stokes_A31(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp);



// First Row-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inline
double A01_00(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g0+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g0x));
}	 

inline
double A01_01(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g1+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g1x));
}	 

inline
double A01_02(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g2+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g2x));
}	 

inline
double A01_03(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g3+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g3x));
}	 

inline
double A01_04(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g4+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g4x));
}	 

inline
double A01_05(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g5+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g5x));
}	

inline
double A01_06(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g6+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g6x));
}	  

inline
double A01_07(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g7+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g7x));
}	

inline
double A01_08(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g8+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g8x));
}	  

inline
double A01_09(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g9+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g9x));
}	 

inline
double A01_10(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g0+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g0x));
}	 

inline
double A01_11(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g1+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g1x));
}	 

inline
double A01_12(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g2+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g2x));
}	

inline
double A01_13(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g3+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g3x));
}	 

inline
double A01_14(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g4+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g4x));
}

inline
double A01_15(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g5+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g5x));
}	 

inline
double A01_16(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g6+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g6x));
}	   

inline
double A01_17(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g7+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g7x));
}	

inline
double A01_18(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g8+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g8x));
}	  

inline
double A01_19(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g9+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g9x));
}	

inline
double A01_20(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g0+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g0x));
}	  

inline
double A01_21(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g1+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g1x));
}	  

inline
double A01_22(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g2+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g2x));
}	

inline
double A01_23(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g3+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g3x));
}	  

inline
double A01_24(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g4+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g4x));
}	

inline
double A01_25(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g5+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g5x));
}	  

inline
double A01_26(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g6+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g6x));
}	 

inline
double A01_27(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g7+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g7x));
}

inline
double A01_28(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g8+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g8x));
}	

inline
double A01_29(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g9+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g9x));
}	 

inline
double A01_30(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g0+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g0x));
}	  

inline
double A01_31(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g1+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g1x));
}	 

inline
double A01_32(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g2+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g2x));
}	  

inline
double A01_33(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g3+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g3x));
}		

inline
double A01_34(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g4+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g4x));
}	   

inline
double A01_35(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g5+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g5x));
}	  

inline
double A01_36(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g6+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g6x));
}	

inline
double A01_37(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g7+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g7x));
}	

inline
double A01_38(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g8+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g8x));
}	  

inline
double A01_39(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)*g9+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g9x));
}	
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
inline
double A02_00(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g0+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g0y));
}	

inline
double A02_01(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g1+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g1y));
}	

inline
double A02_02(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g2+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g2y));
}	

inline
double A02_03(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g4+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g3y));
}	

inline
double A02_04(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g4+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g4y));
}	

inline
double A02_05(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g5+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g5y));
}	

inline
double A02_06(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g6+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g6y));
}	

inline
double A02_07(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g7+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g7y));
}	

inline
double A02_08(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g8+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g8y));
}	

inline
double A02_09(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g9+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g9y));
}	
//#################################################################################################################################################################################
inline
double A02_10(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g0+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g0y));
}	

inline
double A02_11(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g1+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g1y));
}	

inline
double A02_12(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g2+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g2y));
}	

inline
double A02_13(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g3+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g3y));
}	

inline
double A02_14(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g4+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g4y));
}	

inline
double A02_15(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g5+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g5y));
}	

inline
double A02_16(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g6+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g6y));
}	

inline
double A02_17(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g7+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g7y));
}	

inline
double A02_18(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g8+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g8y));
}	

inline
double A02_19(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g9+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g9y));
}	
//#################################################################################################################################################################################
inline
double A02_20(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g0+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g0y));
}	

inline
double A02_21(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g1+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g1y));
}	

inline
double A02_22(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g2+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g2y));
}	

inline
double A02_23(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g3+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g3y));
}	

inline
double A02_24(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g4+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g4y));
}	

inline
double A02_25(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g5+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g5y));
}	

inline
double A02_26(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g6+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g6y));
}	

inline
double A02_27(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g7+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g7y));
}	

inline
double A02_28(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g8+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g8y));
}	

inline
double A02_29(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g9+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g9y));
}	
//#################################################################################################################################################################################
inline
double A02_30(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g0+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g0y));
}	

inline
double A02_31(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g1+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g1y));
}	

inline
double A02_32(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g2+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g2y));
}	

inline
double A02_33(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g3+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g3y));
}	

inline
double A02_34(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g4+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g4y));
}	

inline
double A02_35(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g5+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g5y));
}	

inline
double A02_36(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g6+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g6y));
}	

inline
double A02_37(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g7+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g7y));
}	

inline
double A02_38(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g8+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g8y));
}	

inline
double A02_39(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0y+r1_m*f1y+r2_m*f2y+r3_m*f3y)*g9+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g9y));
}	
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
inline
double A03_00(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g0)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g0z));
}	

inline
double A03_01(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g1)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g1z));
}	

inline
double A03_02(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g2)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g2z));
}	

inline
double A03_03(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g3)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g3z));
}	

inline
double A03_04(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g4)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g4z));
}	

inline
double A03_05(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g5)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g5z));
}	

inline
double A03_06(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g6)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g6z));
}	

inline
double A03_07(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g7)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g7z));
}	

inline
double A03_08(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g8)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g8z));
}	

inline
double A03_09(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f0*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g9)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g9z));
}	
//#################################################################################################################################################################################
inline
double A03_10(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g0)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g0z));
}	

inline
double A03_11(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g1)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g1z));
}	

inline
double A03_12(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g2)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g2z));
}	

inline
double A03_13(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g3)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g3z));
}	

inline
double A03_14(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g4)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g4z));
}	

inline
double A03_15(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g5)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g5z));
}	

inline
double A03_16(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g6)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g6z));
}	

inline
double A03_17(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g7)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g7z));
}	

inline
double A03_18(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g8)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g8z));
}	

inline
double A03_19(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f1*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g9)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g9z));
}	
//#################################################################################################################################################################################
inline
double A03_20(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g0)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g0z));
}	

inline
double A03_21(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g1)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g1z));
}	

inline
double A03_22(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g2)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g2z));
}	

inline
double A03_23(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g3)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g3z));
}	

inline
double A03_24(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g4)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g4z));
}	

inline
double A03_25(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g5)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g5z));
}	

inline
double A03_26(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g6)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g6z));
}	

inline
double A03_27(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g7)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g7z));
}	

inline
double A03_28(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g8)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g8z));
}	

inline
double A03_29(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f2*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g9)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g9z));
}	
//#################################################################################################################################################################################
inline
double A03_30(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g0)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g0z));
}	

inline
double A03_31(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g1)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g1z));
}	

inline
double A03_32(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g2)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g2z));
}	

inline
double A03_33(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g3)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g3z));
}	

inline
double A03_34(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g4)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g4z));
}	

inline
double A03_35(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g5)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g5z));
}	

inline
double A03_36(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g6)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g6z));
}	

inline
double A03_37(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g7)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g7z));
}	

inline
double A03_38(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g8)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g8z));
}	

inline
double A03_39(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return f3*((r0_m*f0z+r1_m*f1z+r2_m*f2z+r3_m*f3z)*(g9)+(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(g9z));
}	
//#################################################################################################################################################################################
inline
double B0_0(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
            u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
            u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
            u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
            u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
            u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
            u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
            u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
            u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
            u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
            
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0 = f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z, // g0x = (\partial g0)/(\partial x) = a1(\partial g0)/(\partial \xi)+a2(\partial g0)/(\partial \eta) +a3(\partial g0)/(\partial \zeta)
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z, // g0y = (\partial g0)/(\partial y) = b1(\partial g0)/(\partial \xi)+b2(\partial g0)/(\partial \eta) +b3(\partial g0)/(\partial \zeta)
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z; // g0z = (\partial g0)/(\partial z) = c1(\partial g0)/(\partial \xi)+c2(\partial g0)/(\partial \eta) +c3(\partial g0)/(\partial \zeta)
                 
	return ((1/dt)*((r0_n-r0_m)*f0+(r1_n-r1_m)*f1+(r2_n-r2_m)*f2+(r3_n-r3_m)*f3)
	       -(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)
	       -(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(r0_m*f0y+r1_m*f2y+r2_m*f2y+r3_m*f3y)
	       -(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(r0_m*f0z+r1_m*f2z+r2_m*f2z+r3_m*f3z)
	       -(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	       -(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(v0_m*g0x+u1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y+v6_m*g6y+v7_m*g7y+v8_m*g8y+v9_m*g9y)
	       -(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(w0_m*g0x+u1_m*g1z+w2_m*g2z+w3_m*g3z+w4_m*g4z+w5_m*g5z+v6_m*g6z+w7_m*g7z+w8_m*g8z+w9_m*g9z))*f0;
}	 

inline
double B0_1(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z, // g0x = (\partial g0)/(\partial x) = a1(\partial g0)/(\partial \xi)+a2(\partial g0)/(\partial \eta) +a3(\partial g0)/(\partial \zeta)
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z, // g0y = (\partial g0)/(\partial y) = b1(\partial g0)/(\partial \xi)+b2(\partial g0)/(\partial \eta) +b3(\partial g0)/(\partial \zeta)
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z; // g0z = (\partial g0)/(\partial z) = c1(\partial g0)/(\partial \xi)+c2(\partial g0)/(\partial \eta) +c3(\partial g0)/(\partial \zeta)
                 
               
	return ((1/dt)*((r0_n-r0_m)*f0+(r1_n-r1_m)*f1+(r2_n-r2_m)*f2+(r3_n-r3_m)*f3)
	       -(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)
	       -(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(r0_m*f0y+r1_m*f2y+r2_m*f2y+r3_m*f3y)
	       -(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(r0_m*f0z+r1_m*f2z+r2_m*f2z+r3_m*f3z)
	       -(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	       -(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(v0_m*g0x+u1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y+v6_m*g6y+v7_m*g7y+v8_m*g8y+v9_m*g9y)
	       -(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(w0_m*g0x+u1_m*g1z+w2_m*g2z+w3_m*g3z+w4_m*g4z+w5_m*g5z+v6_m*g6z+w7_m*g7z+w8_m*g8z+w9_m*g9z))*f1;
}	 

inline
double B0_2(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
                
	return ((1/dt)*((r0_n-r0_m)*f0+(r1_n-r1_m)*f1+(r2_n-r2_m)*f2+(r3_n-r3_m)*f3)
	       -(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)
	       -(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(r0_m*f0y+r1_m*f2y+r2_m*f2y+r3_m*f3y)
	       -(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(r0_m*f0z+r1_m*f2z+r2_m*f2z+r3_m*f3z)
	       -(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	       -(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(v0_m*g0x+u1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y+v6_m*g6y+v7_m*g7y+v8_m*g8y+v9_m*g9y)
	       -(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(w0_m*g0x+u1_m*g1z+w2_m*g2z+w3_m*g3z+w4_m*g4z+w5_m*g5z+v6_m*g6z+w7_m*g7z+w8_m*g8z+w9_m*g9z))*f2;
}	 

inline
double B0_3(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return ((1/dt)*((r0_n-r0_m)*f0+(r1_n-r1_m)*f1+(r2_n-r2_m)*f2+(r3_n-r3_m)*f3)
	       -(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(r0_m*f0x+r1_m*f1x+r2_m*f2x+r3_m*f3x)
	       -(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(r0_m*f0y+r1_m*f2y+r2_m*f2y+r3_m*f3y)
	       -(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(r0_m*f0z+r1_m*f2z+r2_m*f2z+r3_m*f3z)
	       -(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	       -(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(v0_m*g0x+u1_m*g1y+v2_m*g2y+v3_m*g3y+v4_m*g4y+v5_m*g5y+v6_m*g6y+v7_m*g7y+v8_m*g8y+v9_m*g9y)
	       -(r0_m*f0+r1_m*f1+r2_m*f2+r3_m*f3)*(w0_m*g0x+u1_m*g1z+w2_m*g2z+w3_m*g3z+w4_m*g4z+w5_m*g5z+v6_m*g6z+w7_m*g7z+w8_m*g8z+w9_m*g9z))*f3;
}	 
// Second Row-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inline
double A10_00(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g0*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f0;
}	 

inline
double A10_01(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g0*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f1;
}	 

inline
double A10_02(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g0*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f2;
}	 

inline
double A10_03(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g0*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f3;
}	 
//#################################################################################################################################################################################################################################################
inline
double A10_10(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g1*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f0;
}	 

inline
double A10_11(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g1*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f1;
}	 

inline
double A10_12(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g1*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f2;
}	 

inline
double A10_13(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g1*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f3;
}	 
//#################################################################################################################################################################################################################################################
inline
double A10_20(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g2*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f0;
}	 

inline
double A10_21(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g2*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f1;
}	

inline
double A10_22(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g2*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f2;
}	 

inline
double A10_23(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g2*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f3;
}	  
//#################################################################################################################################################################################################################################################
inline
double A10_30(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g3*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f0;
}	 

inline
double A10_31(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g3*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f1;
}

inline
double A10_32(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g3*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f2;
}	 

inline
double A10_33(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g3*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f3;
}	 	 
//#################################################################################################################################################################################################################################################
inline
double A10_40(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g4*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f0;
}	 

inline
double A10_41(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g4*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f1;
}

inline
double A10_42(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g4*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f2;
}	 

inline
double A10_43(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g4*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f3;
}	 	
//#################################################################################################################################################################################################################################################
inline
double A10_50(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g5*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f0;
}	 

inline
double A10_51(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g5*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f1;
}

inline
double A10_52(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g5*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f2;
}	 

inline
double A10_53(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g5*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f3;
}	 	
//#################################################################################################################################################################################################################################################
inline
double A10_60(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g6*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f0;
}	 

inline
double A10_61(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g6*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f1;
}

inline
double A10_62(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g6*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f2;
}	 

inline
double A10_63(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g6*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f3;
}	 	
//#################################################################################################################################################################################################################################################
inline
double A10_70(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g7*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f0;
}	 

inline
double A10_71(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g7*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f1;
}

inline
double A10_72(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g7*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f2;
}	 

inline
double A10_73(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g7*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f3;
}	 	  
//#################################################################################################################################################################################################################################################
inline
double A10_80(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g8*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f0;
}	 

inline
double A10_81(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g8*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f1;
}

inline
double A10_82(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g8*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f2;
}	 

inline
double A10_83(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g8*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f3;
}	 	  
//#################################################################################################################################################################################################################################################
inline
double A10_90(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g9*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f0;
}	 

inline
double A10_91(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g9*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f1;
}

inline
double A10_92(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g9*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f2;
}	 

inline
double A10_93(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return g9*((1/dt)*((u0_m-u0_n)*g0+(u1_m-u1_n)*g1+(u2_m-u2_n)*g2+(u3_m-u3_n)*g3+(u4_m-u4_n)*g4+(u5_m-u5_n)*g5+(u6_m-u6_n)*g6+(u7_m-u7_n)*g7+(u8_m-u8_n)*g8+(u9_m-u9_n)*g9)
	          +(u0_m*g0+u1_m*g1+u2_m*g2+u3_m*g3+u4_m*g4+u5_m*g5+u6_m*g6+u7_m*g7+u8_m*g8+u9_m*g9)*(u0_m*g0x+u1_m*g1x+u2_m*g2x+u3_m*g3x+u4_m*g4x+u5_m*g5x+u6_m*g6x+u7_m*g7x+u8_m*g8x+u9_m*g9x)
	          +(v0_m*g0+v1_m*g1+v2_m*g2+v3_m*g3+v4_m*g4+v5_m*g5+v6_m*g6+v7_m*g7+v8_m*g8+v9_m*g9)*(u0_m*g0y+u1_m*g1y+u2_m*g2y+u3_m*g3y+u4_m*g4y+u5_m*g5y+u6_m*g6y+u7_m*g7y+u8_m*g8y+u9_m*g9y)
	          +(w0_m*g0+w1_m*g1+w2_m*g2+w3_m*g3+w4_m*g4+w5_m*g5+w6_m*g6+w7_m*g7+w8_m*g8+w9_m*g9)*(u0_m*g0z+u1_m*g1z+u2_m*g2z+u3_m*g3z+u4_m*g4z+u5_m*g5z+u6_m*g6z+u7_m*g7z+u8_m*g8z+u9_m*g9z))*f3;
}	 	 
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
inline
double B1_0(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return 0;
}	 
// Third Row-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A20_00(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return 0;
}	 

inline
double B2_0(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return 0;
}	
// Fourth Row-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline
double A30_00(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return 0;
}	 

inline
double B3_0(double x, double  y,double  z)
{ 
	 double jac, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp;
                 
     double r0_n, r0_m, r1_n , r1_m , r2_n, r2_m, r3_n, r3_m;
     
     double u0_n, u0_m, v0_n, v0_m, w0_n, w0_m,
                 u1_n, u1_m, v1_n, v1_m, w1_n, w1_m,
                 u2_n, u2_m, v2_n, v2_m, w2_n, w2_m,
                 u3_n, u3_m, v3_n, v3_m, w3_n, w3_m,
                 u4_n, u4_m, v4_n, v4_m, w4_n, w4_m,
                 u5_n, u5_m, v5_n, v5_m, w5_n, w5_m,
                 u6_n, u6_m, v6_n, v6_m, w6_n, w6_m,
                 u7_n, u7_m, v7_n, v7_m, w7_n, w7_m,
                 u8_n, u8_m, v8_n, v8_m, w8_n, w8_m,
                 u9_n, u9_m, v9_n, v9_m, w9_n, w9_m;
                 
     double f0 = 1-x-y-z, f1 = x, f2 = y, f3 = z;
     
     double f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
            f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
            f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;
            
     
     double g0=f0*(2*f0-1), g1=f1*(2*f1-1), g2=f2*(2*f2-1), g3=f3*(2*f3-1), g4=4*f0*f1, g5=4*f1*f2, g6=4*f0*f2, g7=4*f0*f3, g8=4*f1*f3, g9=4*f2*f3;
     
     
     double g0x = (a1 + a2 + a3)*(4*x + 4*y + 4*z - 3), g1x = a1*(4*x - 1), g2x = a2*(4*y - 1), g3x = a3*(4*z - 1), g4x = - 4*a2*x - 4*a3*x - a1*(8*x + 4*y + 4*z - 4), g5x = 4*a2*x + 4*a1*y, g6x = - 4*a1*y - 4*a3*y - a2*(4*x + 8*y + 4*z - 4), g7x = - 4*a1*z - 4*a2*z - a3*(4*x + 4*y + 8*z - 4), g8x = 4*a3*x + 4*a1*z, g9x = 4*a3*y + 4*a2*z,
            g0y = (b1 + b2 + b3)*(4*x + 4*y + 4*z - 3), g1y = b1*(4*x - 1), g2y = b2*(4*y - 1), g3y = b3*(4*z - 1), g4y = - 4*b2*x - 4*b3*x - b1*(8*x + 4*y + 4*z - 4), g5y = 4*b2*x + 4*b1*y, g6y = - 4*b1*y - 4*b3*y - b2*(4*x + 8*y + 4*z - 4), g7y = - 4*b1*z - 4*b2*z - b3*(4*x + 4*y + 8*z - 4), g8y = 4*b3*x + 4*b1*z, g9y = 4*b3*y + 4*b2*z,
            g0z = (c1 + c2 + c3)*(4*x + 4*y + 4*z - 3), g1z = c1*(4*x - 1), g2z = c2*(4*y - 1), g3z = c3*(4*z - 1), g4z = - 4*c2*x - 4*c3*x - c1*(8*x + 4*y + 4*z - 4), g5z = 4*c2*x + 4*c1*y, g6z = - 4*c1*y - 4*c3*y - c2*(4*x + 8*y + 4*z - 4), g7z = - 4*c1*z - 4*c2*z - c3*(4*x + 4*y + 8*z - 4), g8z = 4*c3*x + 4*c1*z, g9z = 4*c3*y + 4*c2*z;
                 
               
	return 0;
}

// End-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/**
 * Calculates the element mass matrix @p ske.
 * of one triangular element with linear shape functions.
 * @param[in]	ial	node indices of the three element vertices
 * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coordinates of node k
 * @param[out] ske	element stiffness matrix
 */
void CalcElem_Masse(int const ial[4], double const xc[], double ske[4][4]);

/**
 * Calculates the matrix @p ske of our system which consists of element stiffness matrix and mass matrix and the element load vector @p fe
 * of one triangular element with linear shape functions.
 * @param[in] ial node indices of the three element vertices
 * @param[in] xc  vector of node coordinates with x(2*k,2*k+1) as coordinates of node k
 * @param[in] u_old  solution of previous time step
 * @param[in] dt  size of time step
 * @param[in] t  current time
 * @param[in] c  velocity in heat equation
 * @param[out] ske  element stiffness matrix
 * @param[out] fe element load vector
 */
void CalcElem_heat_equation_crank_nichelson(int const ial[4], double const xc[], double ske[4][4], double fe[4], const std::vector<double> &u_old, 
const double dt, const double t, const double c);

/**
 * Adds the element stiffness matrix @p ske and the element load vector @p fe
 * of one triangular element with linear shape functions to the appropriate positions in
 * the symmetric stiffness matrix, stored as CSR matrix K(@p sk,@p id, @p ik)
 *
 * @param[in] ial   node indices of the three element vertices
 * @param[in] ske	element stiffness matrix
 * @param[in] fe	element load vector
 * @param[out] sk	vector non-zero entries of CSR matrix
 * @param[in] id	index vector containing the first entry in a CSR row
 * @param[in] ik	column index vector of CSR matrix
 * @param[out] f	distributed local vector storing the right hand side
 *
 * @warning Algorithm requires indices in connectivity @p ial in ascending order.
 *          Currently deprecated.
*/
void AddElem(int const ial[4], double const ske[4][4], double const fe[4],
             int const id[], int const ik[], double sk[], double f[]);


// #####################################################################
/**
 * Square matrix in CRS format (compressed row storage; also named CSR),
 * see an <a href="https://en.wikipedia.org/wiki/Sparse_matrix">introduction</a>.
 */



class CRS_Matrix
{
    public:
       /**
        * Intializes the CRS matrix structure from the given discetization in @p mesh.
        *
        * The sparse matrix pattern is generated but the values are 0.
        *
        * @param[in] mesh given discretization
        *
        * @warning A reference to the discretization @p mesh is stored inside this class.
        *          Therefore, changing @p mesh outside requires also
        *          to call method @p Derive_Matrix_Pattern explicitely.
        *
        * @see Derive_Matrix_Pattern
        */
       explicit CRS_Matrix(Mesh const & mesh,int ndof_v=2);
       //explicit FEM_Matrix(Mesh const & mesh, int ndof_v=1);
      

      /**
        * Destructor.
        */
       ~CRS_Matrix()
       {}

       /**
        * Generates the sparse matrix pattern and overwrites the existing pattern.
        *
        * The sparse matrix pattern is generated but the values are 0.
       */
       void Derive_Matrix_Pattern();

        /**
        * Calculates the entries of f.e. stiffness matrix and load/rhs vector @p f for the Laplace operator in 2D.
        * No memory is allocated.
        *
        * @param[in,out] f (preallocated) rhs/load vector
        */
       void CalculateLaplace(std::vector<double> &f);
       
       void Skalar2VectorMatrix(int ndof_v);
       
       /**
        * Calculates the entries of f.e. matrix (stiffnes matrix + mass matrix using Crank-Nichelson implicit scheme) load/rhs vector @p f.
        *
        * @param[in,out] f (preallocated) rhs/load vector
        */
       void CalculateLaplace_heat_equation(std::vector<double> &f, const std::vector<double> &u_old, const double dt, const double t, const double c);
       //void CalculateLaplace_heat_equation(vector<double> &f, const vector<double> &u_old, const double dt, const double t, const double c);

       /**
        * Applies Dirichlet boundary conditions to stiffness matrix and to load vector @p f.
        * The <a href="https://www.jstor.org/stable/2005611?seq=1#metadata_info_tab_contents">penalty method</a>
        * is used for incorporating the given values @p u.
        *
        * @param[in]     u (global) vector with Dirichlet data
        * @param[in,out] f load vector
        */
       void ApplyDirichletBC(std::vector<double> const &u, std::vector<double> &f);

       /**
        * Extracts the diagonal elemenst of the sparse matrix.
        *
        * @param[in,out]  d  (prellocated) vector of diagonal elements
        */
       void GetDiag(std::vector<double> &d) const;

       /**
        * Performs the matrix-vector product  w := K*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     u vector
        */
       void Mult(std::vector<double> &w, std::vector<double> const &u) const;

        /**
        * Calculates the defect/residuum w := f - K*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     f load vector
        * @param[in]     u vector
        */
       void Defect(std::vector<double> &w,
                   std::vector<double> const &f, std::vector<double> const &u) const;

       /**
		 * Number rows in matrix.
		 * @return number of rows.
		 */
       int Nrows() const
          {return _nrows;}

       /**
		 * Show the matrix entries.
		 */
       void Debug() const;

       /**
		 * Finds in a CRS matrix the access index for an entry at row @p row and column @p col.
		 *
		 * @param[in] row	row index
		 * @param[in] col	column index
		 * @return index for element (@p row, @p col). If no appropriate entry exists then -1 will be returned.
		 *
		 * @warning assert() stops the function in case that matrix element (@p row, @p col) doesn't exist.
		*/
       int fetch(int row, int col) const;

      /**
        * Adds the element stiffness matrix @p ske and the element load vector @p fe
        * of one triangular element with linear shape functions to the appropriate positions in
        * the stiffness matrix, stored as CSR matrix K(@p sk,@p id, @p ik).
        *
        * @param[in]     ial   node indices of the three element vertices
        * @param[in]     ske   element stiffness matrix
        * @param[in]     fe    element load vector
        * @param[in,out] f	   distributed local vector storing the right hand side
        *
        * @warning Algorithm assumes  linear triangular elements (ndof_e==3).
       */
       void AddElem_3(int const ial[4], double const ske[4][4], double const fe[4], std::vector<double> &f);

        /**
        * Compare @p this CRS matrix with an external CRS matrix stored in C-Style.
        *
        * The method prints statements on differences found.
        *
        * @param[in]     nnode  row number of external matrix
        * @param[in]     id     start indices of matrix rows of external matrix
        * @param[in]     ik     column indices of external matrix
        * @param[in]     sk     non-zero values of external matrix
        *
        * @return true iff all data are identical.
        */
       bool Compare2Old(int nnode, int const id[], int const ik[], double const sk[]) const;
       
       int NdofsVertex() const
    {
        return _ndof_v;
    }

    private:
       Mesh const & _mesh;      //!< reference to discretization
       int _nrows;              //!< number of rows in matrix
       int _ncols;              //!< number of columns in matrix
       int _nnz;                //!< number of non-zero entries
       std::vector<int> _id;    //!< start indices of matrix rows
       std::vector<int> _ik;    //!< column indices
       std::vector<double> _sk; //!< non-zero values
       int  const _ndof_v;      //!< degrees of freedom per vertex (vector valued problems)

};


// ############################################################################################################################################################################################################


class Matrix
{
    public:
       /**
		 * Constructor for abstract matrix class.
         *
         * No memory is allocated.
		 *
		 * @param[in] nrows   number of matrix rows.
		 * @param[in] ncols   number of matrix columns.
		*/
       Matrix(int nrows, int ncols);
       //Matrix();
       
       
       //explicit CRS_Matrix(Mesh const & mesh,int ndof_v=5);

       Matrix(Matrix const &) = default;
       /**
		 * Destructor.
         *
         * No memory is allocated.
		*/
       virtual ~Matrix();

       /**
		 * Checks whether the matrix is a square matrix.
		 *
		 * @return True iff square matrix.
		*/
       bool isSquare() const
       { return _nrows==_ncols;}

       /**
		 * Number of rows in matrix.
		 * @return number of rows.
		 */
       int Nrows() const
          {return _nrows;}

       /**
		 * Number of columns in matrix.
		 * @return number of columns.
		 */
       int Ncols() const
          {return _ncols;}

       /**
		 * Show the matrix entries.
		 */
       virtual void Debug() const = 0;
       
       //void Skalar2VectorMatrix(int ndof_v);

       /**
        * Extracts the diagonal elements of an inherited matrix.
        *
        * @param[in,out]  d  (prellocated) vector of diagonal elements
        */
       virtual void GetDiag(std::vector<double> &d) const = 0;

       /**
        * Extracts the diagonal elements of the matrix.
        *
        * @return  d  vector of diagonal elements
        */
       std::vector<double> const & GetDiag() const
       {
           if ( 0==_dd.size() )          // GH: better?   Nrows()>static_cast<int>(_dd.size())
           {
               _dd.resize(Nrows());
               this->GetDiag(_dd);
               std::cout << "PPPPPPPPPPPPPPPPPPPP\n";
            }
            assert( Nrows()==static_cast<int>(_dd.size()) );
            std::cout << ".";
            return _dd;
        }

       /**
        * Performs the matrix-vector product  w := K*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     u vector
        */
       virtual void Mult(std::vector<double> &w, std::vector<double> const &u) const = 0;

        /**
        * Calculates the defect/residuum w := f - K*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     f load vector
        * @param[in]     u vector
        */
       virtual void Defect(
                   std::vector<double> &w,
                   std::vector<double> const &f, std::vector<double> const &u) const = 0;

       virtual void JacobiSmoother(std::vector<double> const &f, std::vector<double> &u,
                    std::vector<double> &r, int nsmooth, double const omega, bool zero) const
       {
           std::cout << "ERROR in Matrix::JacobiSmoother" << std::endl;
           assert(false);
       }
       
      

       /**
		 * Finds in a CRS matrix the access index for an entry at row @p row and column @p col.
		 *
		 * @param[in] row	row index
		 * @param[in] col	column index
		 * @return index for element (@p row, @p col). If no appropriate entry exists then -1 will be returned.
		 *
		 * @warning assert() stops the function in case that matrix element (@p row, @p col) doesn't exist.
		*/
       virtual int fetch(int row, int col) const =0;

    protected:
       int _nrows;              //!< number of rows in matrix
       int _ncols;              //!< number of columns in matrix
       mutable std::vector<double> _dd; //!< diagonal matrix elements
};

// #####################################################################
class BisectInterpolation;  // class forward declaration
/**
 * Matrix in CRS format (compressed row storage; also named CSR),
 * see an <a href="https://en.wikipedia.org/wiki/Sparse_matrix">introduction</a>.
 */






class CRS_Matrix1: public Matrix
{
    public:
       /**
        * Constructor
        *
        */
       CRS_Matrix1();
       
//! \brief The sparse matrix in CRS format is initialized from a binary file.
//!
//!        The binary file has to store 4 Byte integers and 8 Byte doubles and contains the following data:
//!        - Number of rows
//!        - Number of non-zero elements/blocks
//!        - Number of non-zero matrix elements (= previous number * dofs per block)
//!        - [#elements per row] (counter)
//!        - [column indices]
//!        - [matrix elements]
//!
//! \param[in]   file name of binary file
//!
       CRS_Matrix1(const std::string& file);

       CRS_Matrix1(CRS_Matrix1 const &) = default;


      /**
        * Destructor.
        */
       virtual ~CRS_Matrix1() override;
       /**
        * Extracts the diagonal elements of the sparse matrix.
        *
        * @param[in,out]  d  (prellocated) vector of diagonal elements
        */
       void GetDiag(std::vector<double> &d) const override;

       /**
        * Performs the matrix-vector product  w := K*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     u vector
        */
       void Mult(std::vector<double> &w, std::vector<double> const &u) const override;

        /**
        * Calculates the defect/residuum w := f - K*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     f load vector
        * @param[in]     u vector
        */
       void Defect(std::vector<double> &w,
                   std::vector<double> const &f, std::vector<double> const &u) const override;

       void JacobiSmoother(std::vector<double> const &f, std::vector<double> &u,
                    std::vector<double> &r, int nsmooth, double const omega, bool zero) const override;

       /**
		 * Show the matrix entries.
		 */
       void Debug() const override;

       /**
		 * Finds in a CRS matrix the access index for an entry at row @p row and column @p col.
		 *
		 * @param[in] row	row index
		 * @param[in] col	column index
		 * @return index for element (@p row, @p col). If no appropriate entry exists then -1 will be returned.
		 *
		 * @warning assert() stops the function in case that matrix element (@p row, @p col) doesn't exist.
		*/
       int fetch(int row, int col) const override;

        /**
        * Compare @p this CRS matrix with an external CRS matrix stored in C-Style.
        *
        * The method prints statements on differences found.
        *
        * @param[in]     nnode  row number of external matrix
        * @param[in]     id     start indices of matrix rows of external matrix
        * @param[in]     ik     column indices of external matrix
        * @param[in]     sk     non-zero values of external matrix
        *
        * @return true iff all data are identical.
        */
       bool Compare2Old(int nnode, int const id[], int const ik[], double const sk[]) const;
              
       /**
		 * Calculates the defect and projects it to the next coarser level @f$ f_C := P^T \cdot (f_F - SK\cdot u_F) @f$.  
		 *
		 * @param[in] SK	matrix on fine mesh
		 * @param[in] P	    prolongation operator
		 * @param[in,out] fc  resulting coarse mesh vector (preallocated)
		 * @param[in] ff	r.h.s. on fine mesh
		 * @param[in] uf	status vector on fine mesh 
		 *
		*/
       friend void DefectRestrict(CRS_Matrix1 const & SK1, BisectInterpolation const& P, 
       std::vector<double> &fc, std::vector<double> &ff, std::vector<double> &uf);
       
//! \brief A sparse matrix in CRS format (counter, column index, value) is written to a binary file.
//!
//!        The binary file has to store 4 Byte integers and 8 Byte doubles and contains the following data:
//!        - Number of rows
//!        - Number of non-zero elements
//!        - Number of non-zero elements
//!        - [#elements per row]
//!        - [column indices]
//!        - [elements]
//!
//! \param[in]  file name of binary file
//!
        void writeBinary(const std::string& file);
        
    private:
//! \brief A sparse matrix in CRS format (counter, column index, value) is read from a binary file.
//!
//!        The binary file has to store 4 Byte integers and 8 Byte doubles and contains the following data:
//!        - Number of rows
//!        - Number of non-zero elements/blocks
//!        - Number of non-zero matrix elements (= previous number * dofs per block)
//!        - [#elements per row]
//!        - [column indices]
//!        - [matrix elements]
//!
//! \param[in]   file name of binary file
//!
        void readBinary(const std::string& file);       
    
    public:

    protected:
       int _nnz;                //!< number of non-zero entries
       std::vector<int> _id;    //!< start indices of matrix rows
       std::vector<int> _ik;    //!< column indices
       std::vector<double> _sk; //!< non-zero values
};







/**
 * FEM Matrix in CRS format (compressed row storage; also named CSR),
 * see an <a href="https://en.wikipedia.org/wiki/Sparse_matrix">introduction</a>.
 */
class FEM_Matrix: public CRS_Matrix1
{
    public:
       /**
        * Initializes the CRS matrix structure from the given discretization in @p mesh.
        *
        * The sparse matrix pattern is generated but the values are 0.
        *
        * @param[in] mesh given discretization
        *
        * @warning A reference to the discretization @p mesh is stored inside this class.
        *          Therefore, changing @p mesh outside requires also
        *          to call method @p Derive_Matrix_Pattern explicitly.
        *
        * @see Derive_Matrix_Pattern
        */
      
      explicit FEM_Matrix(Mesh const & mesh, int ndof_v=4);
      
      
       FEM_Matrix(FEM_Matrix const &) = default;

      /**
        * Destructor.
        */
       ~FEM_Matrix() override;

       /**
        * Generates the sparse matrix pattern and overwrites the existing pattern.
        *
        * The sparse matrix pattern is generated but the values are 0.
       */
       void Derive_Matrix_Pattern()
       {
           //Derive_Matrix_Pattern_slow();
           Derive_Matrix_Pattern_fast();
           CheckRowSum();
       }
       void Derive_Matrix_Pattern_fast();
       void Derive_Matrix_Pattern_slow();


        /**
        * Calculates the entries of f.e. stiffness matrix and load/rhs vector @p f for the Laplace operator in 2D.
        * No memory is allocated.
        *
        * @param[in,out] f (preallocated) rhs/load vector
        */
       void CalculateLaplace(std::vector<double> &f);
       
       /**
        * Calculates the entries of f.e. matrix (stiffnes matrix + mass matrix using Crank-Nichelson implicit scheme) load/rhs vector @p f.
        *
        * @param[in,out] f (preallocated) rhs/load vector
        */
       void CalculateLaplace_heat_equation(std::vector<double> &f, const std::vector<double> &u_old, const double dt, const double t, const double c);

       /**
        * Applies Dirichlet boundary conditions to stiffness matrix and to load vector @p f.
        * The <a href="https://www.jstor.org/stable/2005611?seq=1#metadata_info_tab_contents">penalty method</a>
        * is used for incorporating the given values @p u.
        *
        * @param[in]     u (global) vector with Dirichlet data
        * @param[in,out] f load vector
        */
       void ApplyDirichletBC(std::vector<double> const &u, std::vector<double> &f);
       
       void ApplyDirichletBC_Box(std::vector<double> const &u, std::vector<double> &f,
            double xl, double xh, double yl, double yh,double zl, double zh );

       /**
        * Extracts the diagonal elements of the sparse matrix.
        *
        * @param[in,out]  d  (prellocated) vector of diagonal elements
       */
       //void GetDiag(std::vector<double> &d) const;   // override in MPI parallel
       void GetDiag(std::vector<double> &d) const override   { std::cout << "GDDDDDD\n"; GetDiag_M(d); }
       
       // Solves non-M matrix problems for Jacobi iteration but not for MG
       void GetDiag_M(std::vector<double> &d) const;
       
       void Skalar2VectorMatrix(int ndof_v);
       

      /**
        * Adds the element stiffness matrix @p ske and the element load vector @p fe
        * of one triangular element with linear shape functions to the appropriate positions in
        * the stiffness matrix, stored as CSR matrix K(@p sk,@p id, @p ik).
        *
        * @param[in]     ial   node indices of the three element vertices
        * @param[in]     ske   element stiffness matrix
        * @param[in]     fe    element load vector
        * @param[in,out] f	   distributed local vector storing the right hand side
        *
        * @warning Algorithm assumes  linear triangular elements (ndof_e==3).
       */
       void AddElem_3(int const ial[4], double const ske[4][4], double const fe[4], std::vector<double> &f);
       
    /**
     * Global number of degrees of freedom (dof) for each finite element.
     * @return degrees of freedom per element.
     */
    int NdofsVertex() const
    {
        return _ndof_v;
    }
    
    //private:
       bool CheckSymmetry() const;
       bool CheckRowSum() const;
       bool CheckMproperty() const;
       bool CheckMatrix() const;
       
       bool ForceMproperty();

    private:
       Mesh const & _mesh;      //!< reference to discretization
       int  const _ndof_v;      //!< degrees of freedom per vertex (vector valued problems)

};


//**
 //* Prolongation matrix in CRS format (compressed row storage; also named CSR),
 //* see an <a href="https://en.wikipedia.org/wiki/Sparse_matrix">introduction</a>.
 //*
 //* The prolongation is applied for each node from the coarse mesh to the fine mesh and
 //* is derived only geometrically (no operator weighted prolongation).
 //*/
//class Prolongation: public CRS_Matrix
//{
    //public:
       ///**
        //* Intializes the CRS matrix structure from the given discetization in @p mesh.
        //*
        //* The sparse matrix pattern is generated but the values are 0.
        //*
        //* @param[in] cmesh coarse mesh
        //* @param[in] fmesh fine mesh
        //*
        //* @warning A reference to the discretizations @p fmesh  @p cmesh are stored inside this class.
        //*          Therefore, changing these meshes outside requires also
        //*          to call method @p Derive_Matrix_Pattern explicitely.
        //*
        //* @see Derive_Matrix_Pattern
        //*/
       //Prolongation(Mesh const & cmesh, Mesh const & fmesh);

       ///**
        //* Destructor.
        //*/
       //~Prolongation() override
       //{}

       ///**
        //* Generates the sparse matrix pattern and overwrites the existing pattern.
        //*
        //* The sparse matrix pattern is generated but the values are 0.
       //*/
       //void Derive_Matrix_Pattern() override;

    //private:
       //Mesh const & _cmesh;      //!< reference to coarse discretization
       //Mesh const & _fmesh;      //!< reference to fine discretization
//};

// *********************************************************************


/**
 * Interpolation matrix for prolongation coarse mesh (C) to a fine mesh (F)
 * generated by bisecting edges.
 *
 * All interpolation weights are 0.5 (injection points contribute twice).
 */
class BisectInterpolation: public Matrix
{
    public:
       /**
        * Generates the interpolation matrix for prolongation coarse mesh to a fine mesh
        * generated by bisecting edges.
        * The interpolation weights are all 0.5.
        *
        * @param[in] fathers vector[nnodes][2] containing
        *                    the two coarse grid fathers of a fine grid vertex
        *
        */
       explicit BisectInterpolation(std::vector<int> const & fathers);
       BisectInterpolation();

       BisectInterpolation(BisectInterpolation const &) = default;

       /**
        * Destructor.
        */
       ~BisectInterpolation() override;

       /**
        * Extracts the diagonal elements of the matrix.
        *
        * @param[in,out]  d  (prellocated) vector of diagonal elements
        */
       void GetDiag(std::vector<double> &d) const override;
       ///**
        //* Extracts the diagonal elements of the sparse matrix.
        //*
        //* @return  d  vector of diagonal elements
        //*/
       //std::vector<double> const & GetDiag() const override;

       /**
        * Performs the prolongation  @f$ w_F := P*u_C @f$.
        *
        * @param[in,out] wf resulting fine vector (preallocated)
        * @param[in]     uc coarse vector
        */
       void Mult(std::vector<double> &wf, std::vector<double> const &uc) const override;

       /**
        * Performs the restriction  @f$ u_C := P^T*w_F @f$.
        *
        * @param[in]         wf fine vector
        * @param[in,out]     uc resulting coarse vector (preallocated)
        */
       virtual void MultT(std::vector<double> const &wf, std::vector<double> &uc) const;
       
        /**
        * Performs the full restriction  @f$ u_C := F^{-1}*P^T*w_F @f$.
        * 
        * @f$ F @f$ denotes the row sum of the restriction matrix 
        * and results in restricting exactly a bilinear function from the fine grid onto 
        * the same bilinear function on the coarse grid.
        *
        * @param[in]         wf fine vector
        * @param[in,out]     uc resulting coarse vector (preallocated)
        */
       void MultT_Full(std::vector<double> const &wf, std::vector<double> &uc) const;


        /**
        * Calculates the defect/residuum w := f - P*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     f load vector
        * @param[in]     u coarse vector
        */
       void Defect(std::vector<double> &w,
                   std::vector<double> const &f, std::vector<double> const &u) const override;

       /**
		 * Show the matrix entries.
		 */
       void Debug() const override;

       /**
		 * Finds in this matrix the access index for an entry at row @p row and column @p col.
		 *
		 * @param[in] row	row index
		 * @param[in] col	column index
		 * @return index for element (@p row, @p col). If no appropriate entry exists then -1 will be returned.
		 *
		 * @warning assert() stops the function in case that matrix element (@p row, @p col) doesn't exist.
		*/
       int fetch(int row, int col) const override;

       /**
		 * Calculates the defect and projects it to the next coarser level @f$ f_C := P^T \cdot (f_F - SK\cdot u_F) @f$.  
		 *
		 * @param[in] SK	matrix on fine mesh
		 * @param[in] P	    prolongation operator
		 * @param[in,out] fc  resulting coarse mesh vector (preallocated)
		 * @param[in] ff	r.h.s. on fine mesh
		 * @param[in] uf	status vector on fine mesh 
		 *
		*/       
       friend void DefectRestrict(CRS_Matrix1 const & SK1, BisectInterpolation const& P, 
       std::vector<double> &fc, std::vector<double> &ff, std::vector<double> &uf);

    protected:
       std::vector<int> _iv;     //!< fathers[nnode][2] of fine grid nodes, double entries denote injection points
       std::vector<double> _vv;  //!< weights[nnode][2] of fathers for grid nodes
};

/**
 * Interpolation matrix for prolongation from coarse mesh (C)) to a fine mesh (F)
 * generated by bisecting edges.
 * 
 * We take into account that values at Dirichlet nodes have to be preserved, i.e.,
 * @f$ w_F = P \cdot I_D \cdot w_C @f$ and @f$ d_C = I_D  \cdot P^T \cdot  d_F@f$
 * with @f$ I_D @f$ as @f$ n_C \times n_C @f$ diagonal matrix and entries
 * @f$ I_{D(j,j)} := \left\{{\begin{array}{l@{\qquad}l} 0 & x_{j}\;\;  \textrm{is Dirichlet node} \\ 1 & \textrm{else} \end{array}}\right. @f$
 *
 * Interpolation weights are eighter 0.5 or 0.0 in case of coarse Dirichlet nodes
 * (injection points contribute twice),
 * Sets weight to zero iff (at least) one father nodes is a Dirichlet node.
 */
class BisectIntDirichlet: public BisectInterpolation
{
    public:
       /**
		 * Default constructor.
		*/
       BisectIntDirichlet()
        : BisectInterpolation(), _idxDir()
       {}

       /**
		 * Constructs interpolation from father-@p row and column @p col.
		 *
		 * @param[in] fathers	two father nodes from each fine node [nnode_f*2].
		 * @param[in] idxc_dir	vector containing the indices of coarse mesh Dirichlet nodes.
		 *
		*/
       BisectIntDirichlet(std::vector<int> const & fathers, std::vector<int> const & idxc_dir);

       BisectIntDirichlet(BisectIntDirichlet const &) = default;
       
      /**
        * Performs the restriction  @f$ u_C := P^T*w_F @f$.
        *
        * @param[in]         wf fine vector
        * @param[in,out]     uc resulting coarse vector (preallocated)
        */
       virtual void MultT(std::vector<double> const &wf, std::vector<double> &uc) const override;


       /**
        * Destructor.
        */
       ~BisectIntDirichlet() override;
       
    private:
       std::vector<int> const _idxDir;
};



// *********************************************************************

/**
 * Calculates the element stiffness matrix @p ske and the element load vector @p fe
 * of one triangular element with linear shape functions.
 * @param[in]	ial	node indices of the three element vertices
 * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coordinates of node k
 * @param[out] ske	element stiffness matrix
 * @param[out] fe	element load vector
 */
//void CalcElem(int const ial[4], double const xc[], double ske[4][4], double fe[4]);

/**
 * Calculates the element mass matrix @p ske.
 * of one triangular element with linear shape functions.
 * @param[in]	ial	node indices of the three element vertices
 * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coordinates of node k
 * @param[out] ske	element stiffness matrix
 */
//void CalcElem_Masse(int const ial[4], double const xc[], double ske[4][4]);

/**
 * Adds the element stiffness matrix @p ske and the element load vector @p fe
 * of one triangular element with linear shape functions to the appropriate positions in
 * the symmetric stiffness matrix, stored as CSR matrix K(@p sk,@p id, @p ik)
 *
 * @param[in] ial   node indices of the three element vertices
 * @param[in] ske	element stiffness matrix
 * @param[in] fe	element load vector
 * @param[out] sk	vector non-zero entries of CSR matrix
 * @param[in] id	index vector containing the first entry in a CSR row
 * @param[in] ik	column index vector of CSR matrix
 * @param[out] f	distributed local vector storing the right hand side
 *
 * @warning Algorithm requires indices in connectivity @p ial in ascending order.
 *          Currently deprecated.
*/
//void AddElem(int const ial[4], double const ske[4][4], double const fe[4],
             //int const id[], int const ik[], double sk[], double f[]);



#endif


