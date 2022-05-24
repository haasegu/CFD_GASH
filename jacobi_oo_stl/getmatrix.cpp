#include "getmatrix.h"
#include "userset.h"
#include "binaryIO.h"

#include "omp.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <list>
#include <vector>
using namespace std;

//  general routine for lin. triangular elements

void CalcElem(int const ial[4], double const xc[], double ske[4][4], double fe[4])
//void CalcElem(const int* __restrict__ ial, const double* __restrict__ xc, double* __restrict__ ske[3], double* __restrict__ fe)
{
    const int  i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];
    const double x1 = xc[i2 + 0] - xc[i1 + 0],  y1 = xc[i2 + 1] - xc[i1 + 1], z1 = xc[i2 + 2] - xc[i1 + 2],
                 x2 = xc[i3 + 0] - xc[i1 + 0],  y2 = xc[i3 + 1] - xc[i1 + 1], z2 = xc[i3 + 2] - xc[i1 + 2],
                 x3 = xc[i4 + 0] - xc[i1 + 0],  y3 = xc[i4 + 1] - xc[i1 + 1], z3 = xc[i4 + 2] - xc[i1 + 2],
                 x4 = xc[i3 + 0] - xc[i2 + 0],  y4 = xc[i3 + 1] - xc[i2 + 1], z4 = xc[i3 + 2] - xc[i2 + 2],
                 x5 = xc[i2 + 0] - xc[i4 + 0],  y5 = xc[i2 + 1] - xc[i4 + 1], z5 = xc[i2 + 2] - xc[i4 + 2],
                 x6 = xc[i4 + 0] - xc[i3 + 0],  y6 = xc[i4 + 1] - xc[i3 + 1], z6 = xc[i4 + 2] - xc[i3 + 2];
    const double jac = fabs(x3*(y1*z2-y2*z1)+y3*(x2*z1-x1*z2)+z3*(x1*y2-y1*x2));
    

    ske[0][0] = (1/(6*jac)) * ((x4*x4+y4*y4+z4*z4)*(x5*x5+y5*y5+z5*z5)-(x4*x5+y4*y5+z4*z5)*(x4*x5+y4*y5+z4*z5));
    ske[0][1] = (1/(6*jac)) * ((x2*x6+y2*y6+z2*z6)*(x4*x6+y4*y6+z4*z6)-(x2*x4+y2*y4+z2*z4)*(x6*x6+y6*y6+z6*z6));
    ske[0][2] = (1/(6*jac)) * ((x3*x4+y3*y4+z3*z4)*(x5*x5+y5*y5+z5*z5)-(x3*x5+y3*y5+z3*z5)*(x4*x5+y4*y5+z4*z5));
    ske[0][3] = (1/(6*jac)) * ((x1*x6+y1*y6+z1*z6)*(x4*x4+y4*y4+z4*z4)-(x1*x4+y1*y4+z1*z4)*(x4*x6+y4*y6+z4*z6));
    ske[1][0] = ske[0][1];
    ske[1][1] = (1/(6*jac)) * ((x2*x2+y2*y2+z2*z2)*(x3*x3+y3*y3+z3*z3)-(x2*x3+y2*y3+z2*z3)*(x2*x3+y2*y3+z2*z3));
    ske[1][2] = (1/(6*jac)) * ((x1*x6+y1*y6+z1*z6)*(x3*x3+y3*y3+z3*z3)-(x1*x3+y1*y3+z1*z3)*(x3*x6+y3*y6+z3*z6));
    ske[1][3] = (1/(6*jac)) * ((x1*x3+y1*y3+z1*z3)*(x4*x6+y4*y6+z4*z6)-(x1*x6+y1*y6+z1*z6)*(x3*x4+y3*y4+z3*z4));
    ske[2][0] = ske[0][2];
    ske[2][1] = ske[1][2];
    ske[2][2] = (1/(6*jac)) * ((x1*x1+y1*y1+z1*z1)*(x3*x3+y3*y3+z3*z3)-(x1*x3+y1*y3+z1*z3)*(x1*x3+y1*y3+z1*z3));
    ske[2][3] = (1/(6*jac)) * ((x1*x3+y1*y3+z1*z3)*(x1*x4+y1*y4+z1*z4)-(x1*x1+y1*y1+z1*z1)*(x3*x4+y3*y4+z3*z4));
    ske[3][0] = ske[0][3];
    ske[3][1] = ske[1][3];
    ske[3][2] = ske[2][3];
    ske[3][3] = (1/(6*jac)) * ((x1*x1+y1*y1+z1*z1)*(x4*x4+y4*y4+z4*z4)-(x1*x4+y1*y4+z1*z4)*(x1*x4+y1*y4+z1*z4));

    const double xm    = (xc[i1 + 0] + xc[i2 + 0] + xc[i3 + 0]+ xc[i4 + 0]) / 4.0,
                 ym    = (xc[i1 + 1] + xc[i2 + 1] + xc[i3 + 1]+ xc[i4 + 1]) / 4.0,
                 zm    = (xc[i1 + 2] + xc[i2 + 2] + xc[i3 + 2]+ xc[i4 + 2]) / 4.0;
    //fe[0] = fe[1] = fe[2] = 0.5 * jac * FunctF(xm, ym) / 3.0;
    fe[0] = fe[1] = fe[2]= fe[3] =  0*(jac/6.0) * fNice(0, xm, ym, zm) / 4.0;
}

void CalcElem_Masse(int const ial[4], double const xc[], double const cm, double ske[4][4])
{
    const int  i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];
    const double x1 = xc[i2 + 0] - xc[i1 + 0],  y1 = xc[i2 + 1] - xc[i1 + 1], z1 = xc[i2 + 2] - xc[i1 + 2],
                 x2 = xc[i3 + 0] - xc[i1 + 0],  y2 = xc[i3 + 1] - xc[i1 + 1], z2 = xc[i3 + 2] - xc[i1 + 2],
                 x3 = xc[i4 + 0] - xc[i1 + 0],  y3 = xc[i4 + 1] - xc[i1 + 1], z3 = xc[i4 + 2] - xc[i1 + 2],
                 x4 = xc[i3 + 0] - xc[i2 + 0],  y4 = xc[i3 + 1] - xc[i2 + 1], z4 = xc[i3 + 2] - xc[i2 + 2],
                 x5 = xc[i2 + 0] - xc[i4 + 0],  y5 = xc[i2 + 1] - xc[i4 + 1], z5 = xc[i2 + 2] - xc[i4 + 2],
                 x6 = xc[i4 + 0] - xc[i3 + 0],  y6 = xc[i4 + 1] - xc[i3 + 1], z6 = xc[i4 + 2] - xc[i3 + 2];
    const double jac = fabs(x3*(y1*z2-y2*z1)+y3*(x2*z1-x1*z2)+z3*(x1*y2-y1*x2));

    ske[0][0] += jac/12.0;
    ske[0][1] += jac/24.0;
    ske[0][2] += jac/24.0;
    ske[0][3] += jac/24.0;
    ske[1][0] += jac/24.0;
    ske[1][1] += jac/12.0;
    ske[1][2] += jac/24.0;
    ske[1][3] += jac/24.0;
    ske[2][0] += jac/24.0;
    ske[2][1] += jac/24.0;
    ske[2][2] += jac/12.0;
    ske[2][3] += jac/24.0;
    ske[3][0] += jac/24.0;
    ske[3][1] += jac/24.0;
    ske[3][2] += jac/24.0;
    ske[3][3] += jac/12.0;

    return;
}

// The following functions are added by Salman Ahmad. 

//double B11(int const ial[10], double x, double  y,double  z, const std::vector<double> &u_old_m){
	//return (1-x-y-z)*u_old_m.at(ial[9])-4*(1-x-y-z);
	//}

//A00
void CalcElem_Navier_Stokes_A00(int const ial[4], double const xc[], double ske[4][4], double fe[4], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
                 
    const double r0_n = r_old_n.at(ial[0]), r0_m = r_old_m.at(ial[0]),
                 r1_n = r_old_n.at(ial[1]), r1_m = r_old_m.at(ial[1]),
                 r2_n = r_old_n.at(ial[2]), r2_m = r_old_m.at(ial[2]),
                 r3_n = r_old_n.at(ial[3]), r3_m = r_old_m.at(ial[3]);
     
    const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
		
	double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;

    ske[0][0] = (jac/dt)*0.0167+jac*((a1+a2+a3)*(2.8969e-07)-0.0167*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(-(b1+b2+b3)*(-2.8969e-07)-0.0167*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(-(c1+c2+c3)*(-2.8969e-07)-0.0167*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*((a1+a2+a3)*(0.0028)+a1*0.0028)*u_old_m.at(ial[1])+jac*((b1+b2+b3)*(0.0028)+b1*0.0028)*v_old_m.at(ial[1])+jac*((c1+c2+c3)*(0.0028)+c1*0.0028)*w_old_m.at(ial[1])
                               +jac*((a1+a2+a3)*(0.0028)-a2*0.0028)*u_old_m.at(ial[2])+jac*((b1+b2+b3)*(0.0028)-b2*0.0028)*v_old_m.at(ial[2])+jac*((c1+c2+c3)*(0.0028)-c2*0.0028)*w_old_m.at(ial[2])
                               +jac*((a1+a2+a3)*(0.0028)-a3*0.0028)*u_old_m.at(ial[3])+jac*((b1+b2+b3)*(0.0028)-c3*0.0028)*v_old_m.at(ial[3])+jac*((c1+c2+c3)*(0.0028)-c3*0.0028)*w_old_m.at(ial[3])
                               +jac*(-(a1+a2+a3)*0.0111+5.6703e-19*a1-0.0111*a2-0.0111*a3)*u_old_m.at(ial[4])+jac*(-(b1+b2+b3)*0.0111+5.6703e-19*b1-0.0111*b2-0.0111*b3)*v_old_m.at(ial[4])+jac*(-(c1+c2+c3)*0.0111+5.6703e-19*c1-0.0111*c2-0.0111*c3)*w_old_m.at(ial[4])
                               +jac*(-(a1+a2+a3)*0.0056+0.0111*a1+0.0111*a2)*u_old_m.at(ial[5])+jac*(-(b1+b2+b3)*0.0056+0.0111*b1+0.0111*b2)*v_old_m.at(ial[5])+jac*(-(c1+c2+c3)*0.0056+0.0111*c1+0.0111*c2)*w_old_m.at(ial[5])
                               +jac*(-(a1+a2+a3)*0.0111-0.0111*a1+0.0222*a2-0.0111*a3)*u_old_m.at(ial[6])+jac*(-(b1+b2+b3)*0.0111-0.0111*b1+0.0222*b2-0.0111*b3)*v_old_m.at(ial[6])+jac*(-(c1+c2+c3)*0.0111-0.0111*c1+0.0222*c2-0.0111*c3)*w_old_m.at(ial[6])
                               +jac*(-(a1+a2+a3)*0.0111-0.0111*a1-0.0111*a2-0.0222*a3)*u_old_m.at(ial[7])+jac*(-(b1+b2+b3)*0.0111-0.0111*b1-0.0111*b2-0.0222*b3)*v_old_m.at(ial[7])+jac*(-(c1+c2+c3)*0.0111-0.0111*c1-0.0111*c2-0.0222*c3)*w_old_m.at(ial[7])
                               +jac*(-(a1+a2+a3)*0.0056+0.0111*a1+0.0111*a3)*u_old_m.at(ial[8])+jac*(-(b1+b2+b3)*0.0056+0.0111*b1+0.0111*b3)*v_old_m.at(ial[8])+jac*(-(c1+c2+c3)*0.0056+0.0111*c1+0.0111*c3)*w_old_m.at(ial[8])
                               +jac*(-(a1+a2+a3)*0.0056+0.0111*a2+0.0111*a3)*u_old_m.at(ial[9])+jac*(-(b1+b2+b3)*0.0056+0.0111*b2+0.0111*b3)*v_old_m.at(ial[9])+jac*(-(c1+c2+c3)*0.0056+0.0111*c2+0.0111*c3)*w_old_m.at(ial[9]);
                               
    ske[0][1] = (jac/dt)*0.0083+jac*(a1*(-2.8969e-07)-0.0028*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(b1*(-2.8969e-07)-0.0028*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(c1*(-2.8969e-07)-0.0028*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(a1*(-0.0028)-a1*0.0056)*u_old_m.at(ial[1])+jac*(b1*(-0.0028)-b1*0.0056)*v_old_m.at(ial[1])+jac*(c1*(-0.0028)-c1*0.0056)*w_old_m.at(ial[1])
                               +jac*(a1*(-0.0028)-a2*0.0056)*u_old_m.at(ial[2])+jac*(b1*(-0.0028)-b2*0.0056)*v_old_m.at(ial[2])+jac*(c1*(-0.0028)-c2*0.0056)*w_old_m.at(ial[2])
                               +jac*(a1*(-0.0028)-a3*0.0056)*u_old_m.at(ial[3])+jac*(b1*(-0.0028)-c3*0.0056)*v_old_m.at(ial[3])+jac*(c1*(-0.0028)-c3*0.0056)*w_old_m.at(ial[3])
                               +jac*(a1*0.0111+0.0222*a1-0.0111*a2-0.0111*a3)*u_old_m.at(ial[4])+jac*(b1*0.0111+0.0222*b1-0.0111*b2-0.0111*b3)*v_old_m.at(ial[4])+jac*(c1*0.0111+0.0222*c1-0.0111*c2-0.0111*c3)*w_old_m.at(ial[4])
                               +jac*(a1*0.0056+0.0056*a1+0.0111*a2)*u_old_m.at(ial[5])+jac*(b1*0.0056+0.0056*b1+0.0111*b2)*v_old_m.at(ial[5])+jac*(c1*0.0056+0.0056*c1+0.0111*c2)*w_old_m.at(ial[5])
                               +jac*(a1*0.0111-0.0056*a1+0.0056*a2-0.0056*a3)*u_old_m.at(ial[6])+jac*(b1*0.0111-0.0056*b1+0.0056*b2-0.0056*b3)*v_old_m.at(ial[6])+jac*(c1*0.0111-0.0056*c1+0.0056*c2-0.0056*c3)*w_old_m.at(ial[6])
                               +jac*(a1*0.0111-0.0056*a1-0.0056*a2+0.0056*a3)*u_old_m.at(ial[7])+jac*(b1*0.0111-0.0056*b1-0.0056*b2+0.0056*b3)*v_old_m.at(ial[7])+jac*(c1*0.0111-0.0056*c1-0.0056*c2+0.0056*c3)*w_old_m.at(ial[7])
                               +jac*(a1*0.0056+0.0056*a1+0.0111*a3)*u_old_m.at(ial[8])+jac*(b1*0.0056+0.0056*b1+0.0111*b3)*v_old_m.at(ial[8])+jac*(c1*0.0056+0.0056*c1+0.0111*c3)*w_old_m.at(ial[8])
                               +jac*(a1*0.0056+0.0056*a2+0.0056*a3)*u_old_m.at(ial[9])+jac*(b1*0.0056+0.0056*b2+0.0056*b3)*v_old_m.at(ial[9])+jac*(c1*0.0056+0.0056*c2+0.0056*c3)*w_old_m.at(ial[9]);
                               
    ske[0][2] = (jac/dt)*0.0083+jac*(a2*(-2.8969e-07)-0.0028*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(b2*(-2.8969e-07)-0.0028*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(c2*(-2.8969e-07)-0.0028*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(a2*(-0.0028)-a1*0.0028)*u_old_m.at(ial[1])+jac*(b2*(-0.0028)-b1*0.0028)*v_old_m.at(ial[1])+jac*(c2*(-0.0028)-c1*0.0028)*w_old_m.at(ial[1])
                               +jac*(a2*(-0.0028)+a2*0.0028)*u_old_m.at(ial[2])+jac*(b2*(-0.0028)+b2*0.0028)*v_old_m.at(ial[2])+jac*(c2*(-0.0028)+c2*0.0028)*w_old_m.at(ial[2])
                               +jac*(a2*(-0.0028)-a3*0.0028)*u_old_m.at(ial[3])+jac*(b2*(-0.0028)-c3*0.0028)*v_old_m.at(ial[3])+jac*(c2*(-0.0028)-c3*0.0028)*w_old_m.at(ial[3])
                               +jac*(a2*0.0111+0.0056*a1-0.0056*a2-0.0056*a3)*u_old_m.at(ial[4])+jac*(b2*0.0111+0.0056*b1-0.0056*b2-0.0056*b3)*v_old_m.at(ial[4])+jac*(c2*0.0111+0.0056*c1-0.0056*c2-0.0056*c3)*w_old_m.at(ial[4])
                               +jac*(a2*0.0056+0.0111*a1+0.0056*a2)*u_old_m.at(ial[5])+jac*(b2*0.0056+0.0111*b1+0.0056*b2)*v_old_m.at(ial[5])+jac*(c2*0.0056+0.0111*c1+0.0056*c2)*w_old_m.at(ial[5])
                               +jac*(a2*0.0111-0.0111*a1+5.6703e-19*a2-0.0111*a3)*u_old_m.at(ial[6])+jac*(b2*0.0111-0.0111*b1+5.6703e-19*b2-0.0111*b3)*v_old_m.at(ial[6])+jac*(c2*0.0111-0.0111*c1+5.6703e-19*c2-0.0111*c3)*w_old_m.at(ial[6])
                               +jac*(a2*0.0111-0.0056*a1-0.0056*a2+0.0056*a3)*u_old_m.at(ial[7])+jac*(b2*0.0111-0.0056*b1-0.0056*b2+0.0056*b3)*v_old_m.at(ial[7])+jac*(c2*0.0111-0.0056*c1-0.0056*c2+0.0056*c3)*w_old_m.at(ial[7])
                               +jac*(a2*0.0056+0.0056*a1+0.0056*a3)*u_old_m.at(ial[8])+jac*(b2*0.0056+0.0056*b1+0.0056*b3)*v_old_m.at(ial[8])+jac*(c2*0.0056+0.0056*c1+0.0056*c3)*w_old_m.at(ial[8])
                               +jac*(a2*0.0056+0.0056*a2+0.0111*a3)*u_old_m.at(ial[9])+jac*(b2*0.0056+0.0056*b2+0.0111*b3)*v_old_m.at(ial[9])+jac*(c2*0.0056+0.0056*c2+0.0111*c3)*w_old_m.at(ial[9]);
                               
    ske[0][3] = (jac/dt)*0.0083+jac*(a3*(-2.8969e-07)-0.0028*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(b3*(-2.8969e-07)-0.0028*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(c3*(-2.8969e-07)-0.0028*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(a3*(-0.0028)-a1*0.0028)*u_old_m.at(ial[1])+jac*(b3*(-0.0028)-b1*0.0028)*v_old_m.at(ial[1])+jac*(c3*(-0.0028)-c1*0.0028)*w_old_m.at(ial[1])
                               +jac*(a3*(-0.0028)-a2*0.0028)*u_old_m.at(ial[2])+jac*(b3*(-0.0028)-b2*0.0028)*v_old_m.at(ial[2])+jac*(c3*(-0.0028)-c2*0.0028)*w_old_m.at(ial[2])
                               +jac*(a3*(-0.0028)+a3*0.0028)*u_old_m.at(ial[3])+jac*(b3*(-0.0028)+c3*0.0028)*v_old_m.at(ial[3])+jac*(c3*(-0.0028)+c3*0.0028)*w_old_m.at(ial[3])
                               +jac*(a3*0.0111+0.0056*a1-0.0056*a2-0.0056*a3)*u_old_m.at(ial[4])+jac*(b3*0.0111+0.0056*b1-0.0056*b2-0.0056*b3)*v_old_m.at(ial[4])+jac*(c3*0.0111+0.0056*c1-0.0056*c2-0.0056*c3)*w_old_m.at(ial[4])
                               +jac*(a3*0.0056+0.0056*a1+0.0056*a2)*u_old_m.at(ial[5])+jac*(b3*0.0056+0.0056*b1+0.0056*b2)*v_old_m.at(ial[5])+jac*(c3*0.0056+0.0056*c1+0.0056*c2)*w_old_m.at(ial[5])
                               +jac*(a3*0.0111-0.0056*a1+0.0056*a2-0.0056*a3)*u_old_m.at(ial[6])+jac*(b3*0.0111-0.0056*b1+0.0056*b2-0.0056*b3)*v_old_m.at(ial[6])+jac*(c3*0.0111-0.0056*c1+0.0056*c2-0.0056*c3)*w_old_m.at(ial[6])
                               +jac*(a3*0.0111-0.0111*a1-0.0111*a2+5.6703e-19*a3)*u_old_m.at(ial[7])+jac*(b3*0.0111-0.0111*b1-0.0111*b2+5.6703e-19*b3)*v_old_m.at(ial[7])+jac*(c3*0.0111-0.0111*c1-0.0111*c2+5.6703e-19*c3)*w_old_m.at(ial[7])
                               +jac*(a3*0.0056+0.0111*a1+0.0056*a3)*u_old_m.at(ial[8])+jac*(b3*0.0056+0.0111*b1+0.0056*b3)*v_old_m.at(ial[8])+jac*(c3*0.0056+0.0111*c1+0.0056*c3)*w_old_m.at(ial[8])
                               +jac*(a3*0.0056+0.0111*a2+0.0056*a3)*u_old_m.at(ial[9])+jac*(b3*0.0056+0.0111*b2+0.0056*b3)*v_old_m.at(ial[9])+jac*(c3*0.0056+0.0111*c2+0.0056*c3)*w_old_m.at(ial[9]);
                               
    ske[1][0] = 0; (jac/dt)*0.0083+jac*(-(a1+a2+a3)*(-0.0028)+a1*(-0.0028)+a2*(-0.0028)+a3*(-0.0028))*u_old_m.at(ial[0])+jac*(-(b1+b2+b3)*(-0.0028)+b1*(-0.0028)+b2*(-0.0028)+b3*(-0.0028))*v_old_m.at(ial[0])+jac*(-(c1+c2+c3)*(-0.0028)+c1*(-0.0028)+c2*(-0.0028)+c3*(-0.0028))*w_old_m.at(ial[0])
                               +jac*(-(a1+a2+a3)*(-2.8969e-07)+a1*0.0028)*u_old_m.at(ial[1])+jac*(-(b1+b2+b3)*(-2.8969e-07)+b1*0.0028)*v_old_m.at(ial[1])+jac*(-(c1+c2+c3)*(-2.8969e-07)+c1*0.0028)*w_old_m.at(ial[1])
                               +jac*(-(a1+a2+a3)*(-0.0028)-a2*0.0028)*u_old_m.at(ial[2])+jac*(-(b1+b2+b3)*(-0.0028)-b2*0.0028)*v_old_m.at(ial[2])+jac*(-(c1+c2+c3)*(-0.0028)-c2*0.0028)*w_old_m.at(ial[2])
                               +jac*(-(a1+a2+a3)*(-0.0028)-a3*0.0028)*u_old_m.at(ial[3])+jac*(-(b1+b2+b3)*(-0.0028)-c3*0.0028)*v_old_m.at(ial[3])+jac*(-(c1+c2+c3)*(-0.0028)-c3*0.0028)*w_old_m.at(ial[3])
                               +jac*(-(a1+a2+a3)*0.0111-0.0222*a1-0.0333*a2-0.0333*a3)*u_old_m.at(ial[4])+jac*(-(b1+b2+b3)*0.0111-0.0222*b1-0.0333*b2-0.0333*b3)*v_old_m.at(ial[4])+jac*(-(c1+c2+c3)*0.0111-0.0222*c1-0.0333*c2-0.0333*c3)*w_old_m.at(ial[4])
                               +jac*(-(a1+a2+a3)*0.0111+0.0056*a1+0.0111*a2)*u_old_m.at(ial[5])+jac*(-(b1+b2+b3)*0.0111+0.0056*b1+0.0111*b2)*v_old_m.at(ial[5])+jac*(-(c1+c2+c3)*0.0111+0.0056*c1+0.0111*c2)*w_old_m.at(ial[5])
                               +jac*(-(a1+a2+a3)*0.0056-0.0111*a1+0.0056*a2-0.0111*a3)*u_old_m.at(ial[6])+jac*(-(b1+b2+b3)*0.0056-0.0111*b1+0.0056*b2-0.0111*b3)*v_old_m.at(ial[6])+jac*(-(c1+c2+c3)*0.0056-0.0111*c1+0.0056*c2-0.0111*c3)*w_old_m.at(ial[6])
                               +jac*(-(a1+a2+a3)*0.0056-0.0056*a1-0.0056*a2-0.0056*a3)*u_old_m.at(ial[7])+jac*(-(b1+b2+b3)*0.0056-0.0056*b1+0.0111*b2-0.0056*b3)*v_old_m.at(ial[7])+jac*(-(c1+c2+c3)*0.0056-0.0056*c1+0.0111*c2-0.0056*c3)*w_old_m.at(ial[7])
                               +jac*(-(a1+a2+a3)*0.0111+0.0056*a1+0.0111*a3)*u_old_m.at(ial[8])+jac*(-(b1+b2+b3)*0.0111+0.0056*b1+0.0111*b3)*v_old_m.at(ial[8])+jac*(-(c1+c2+c3)*0.0111+0.0056*c1+0.0111*c3)*w_old_m.at(ial[8])
                               +jac*(-(a1+a2+a3)*0.0056+0.0056*a2+0.0056*a3)*u_old_m.at(ial[9])+jac*(-(b1+b2+b3)*0.0056+0.0056*b2+0.0056*b3)*v_old_m.at(ial[9])+jac*(-(c1+c2+c3)*0.0056+0.0056*c2+0.0056*c3)*w_old_m.at(ial[9]);
                               
    ske[1][1] = (jac/dt)*0.0167+jac*(a1*(-0.0028)+0.0056*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(b1*(-0.0028)+0.0056*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(c1*(-0.0028)+0.0056*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(a1*(-2.8969e-07)+a1*0.0167)*u_old_m.at(ial[1])+jac*(b1*(-2.8969e-07)+b1*0.0167)*v_old_m.at(ial[1])+jac*(c1*(-2.8969e-07)+c1*0.0167)*w_old_m.at(ial[1])
                               +jac*(a1*(-0.0028)-a2*0.0056)*u_old_m.at(ial[2])+jac*(b1*(-0.0028)-b2*0.0056)*v_old_m.at(ial[2])+jac*(c1*(-0.0028)-c2*0.0056)*w_old_m.at(ial[2])
                               +jac*(a1*(-0.0028)-a3*0.0056)*u_old_m.at(ial[3])+jac*(b1*(-0.0028)-c3*0.0056)*v_old_m.at(ial[3])+jac*(c1*(-0.0028)-c3*0.0056)*w_old_m.at(ial[3])
                               +jac*(a1*0.0111+0.0222*a1-0.0111*a2-0.0111*a3)*u_old_m.at(ial[4])+jac*(b1*0.0111+0.0222*b1-0.0111*b2-0.0111*b3)*v_old_m.at(ial[4])+jac*(c1*0.0111+0.0222*c1-0.0111*c2-0.0111*c3)*w_old_m.at(ial[4])
                               +jac*(a1*0.0111+0.0111*a1+0.0333*a2)*u_old_m.at(ial[5])+jac*(b1*0.0111+0.0111*b1+0.0333*b2)*v_old_m.at(ial[5])+jac*(c1*0.0111+0.0111*c1+0.0333*c2)*w_old_m.at(ial[5])
                               +jac*(a1*0.0056-0.0111*a1+1.5999e-18*a2-0.0111*a3)*u_old_m.at(ial[6])+jac*(b1*0.0056-0.0111*b1+1.5999e-18*b2-0.0111*b3)*v_old_m.at(ial[6])+jac*(c1*0.0056-0.0111*c1+1.5999e-18*c2-0.0111*c3)*w_old_m.at(ial[6])
                               +jac*(a1*0.0056-0.0111*a1-0.0111*a2-0.0222*a3)*u_old_m.at(ial[7])+jac*(b1*0.0056-0.0111*b1-0.0111*b2-0.0222*b3)*v_old_m.at(ial[7])+jac*(c1*0.0056-0.0111*c1-0.0111*c2-0.0222*c3)*w_old_m.at(ial[7])
                               +jac*(a1*0.0111+0.0111*a1+0.0333*a3)*u_old_m.at(ial[8])+jac*(b1*0.0111+0.0111*b1+0.0333*b3)*v_old_m.at(ial[8])+jac*(c1*0.0111+0.0111*c1+0.0333*c3)*w_old_m.at(ial[8])
                               +jac*(a1*0.0056+0.0111*a2+0.0111*a3)*u_old_m.at(ial[9])+jac*(b1*0.0056+0.0111*b2+0.0111*b3)*v_old_m.at(ial[9])+jac*(c1*0.0056+0.0111*c2+0.0111*c3)*w_old_m.at(ial[9]);
                               
    ske[1][2] = (jac/dt)*0.0083+jac*(a2*(-0.0028)+0.0028*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(b2*(-0.0028)+0.0028*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(c2*(-0.0028)+0.0028*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(a2*(-2.8969e-07)+a1*0.0028)*u_old_m.at(ial[1])+jac*(b2*(-2.8969e-07)+b1*0.0028)*v_old_m.at(ial[1])+jac*(c2*(-2.8969e-07)+c1*0.0028)*w_old_m.at(ial[1])
                               +jac*(a2*(-0.0028)+a2*0.0028)*u_old_m.at(ial[2])+jac*(b2*(-0.0028)+b2*0.0028)*v_old_m.at(ial[2])+jac*(c2*(-0.0028)+c2*0.0028)*w_old_m.at(ial[2])
                               +jac*(a2*(-0.0028)-a3*0.0028)*u_old_m.at(ial[3])+jac*(b2*(-0.0028)-c3*0.0028)*v_old_m.at(ial[3])+jac*(c2*(-0.0028)-c3*0.0028)*w_old_m.at(ial[3])
                               +jac*(a2*0.0111-0.0056*a1-0.0111*a2-0.0111*a3)*u_old_m.at(ial[4])+jac*(b2*0.0111-0.0056*b1-0.0111*b2-0.0111*b3)*v_old_m.at(ial[4])+jac*(c2*0.0111-0.0056*c1-0.0111*c2-0.0111*c3)*w_old_m.at(ial[4])
                               +jac*(a2*0.0111+0.0111*a1+0.0111*a2)*u_old_m.at(ial[5])+jac*(b2*0.0111+0.0111*b1+0.0111*b2)*v_old_m.at(ial[5])+jac*(c2*0.0111+0.0111*c1+0.0111*c2)*w_old_m.at(ial[5])
                               +jac*(a2*0.0056-0.0111*a1-0.0056*a2-0.0111*a3)*u_old_m.at(ial[6])+jac*(b2*0.0056-0.0111*b1-0.0056*b2-0.0111*b3)*v_old_m.at(ial[6])+jac*(c2*0.0056-0.0111*c1-0.0056*c2-0.0111*c3)*w_old_m.at(ial[6])
                               +jac*(a2*0.0056-0.0056*a1-0.0056*a2+1.0835e-18*a3)*u_old_m.at(ial[7])+jac*(b2*0.0056-0.0056*b1-0.0056*b2+1.0835e-18*b3)*v_old_m.at(ial[7])+jac*(c2*0.0056-0.0056*c1-0.0056*c2+1.0835e-18*c3)*w_old_m.at(ial[7])
                               +jac*(a2*0.0111+0.0056*a1+0.0111*a3)*u_old_m.at(ial[8])+jac*(b2*0.0111+0.0056*b1+0.0111*b3)*v_old_m.at(ial[8])+jac*(c2*0.0111+0.0056*c1+0.0111*c3)*w_old_m.at(ial[8])
                               +jac*(a2*0.0056+0.0056*a2+0.0111*a3)*u_old_m.at(ial[9])+jac*(b2*0.0056+0.0056*b2+0.0111*b3)*v_old_m.at(ial[9])+jac*(c2*0.0056+0.0056*c2+0.0111*c3)*w_old_m.at(ial[9]);
                               
    ske[1][3] = (jac/dt)*0.0083+jac*(a3*(-0.0028)+0.0028*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(b3*(-0.0028)+0.0028*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(c3*(-0.0028)+0.0028*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(a3*(-2.8969e-07)+a1*0.0028)*u_old_m.at(ial[1])+jac*(b3*(-2.8969e-07)+b1*0.0028)*v_old_m.at(ial[1])+jac*(c3*(-2.8969e-07)+c1*0.0028)*w_old_m.at(ial[1])
                               +jac*(a3*(-0.0028)-a2*0.0028)*u_old_m.at(ial[2])+jac*(b3*(-0.0028)-b2*0.0028)*v_old_m.at(ial[2])+jac*(c3*(-0.0028)-c2*0.0028)*w_old_m.at(ial[2])
                               +jac*(a3*(-0.0028)+a3*0.0028)*u_old_m.at(ial[3])+jac*(b3*(-0.0028)+c3*0.0028)*v_old_m.at(ial[3])+jac*(c3*(-0.0028)+c3*0.0028)*w_old_m.at(ial[3])
                               +jac*(a3*0.0111-0.0056*a1-0.0111*a2-0.0111*a3)*u_old_m.at(ial[4])+jac*(b3*0.0111-0.0056*b1-0.0111*b2-0.0111*b3)*v_old_m.at(ial[4])+jac*(c3*0.0111-0.0056*c1-0.0111*c2-0.0111*c3)*w_old_m.at(ial[4])
                               +jac*(a3*0.0111+0.0056*a1+0.0111*a2)*u_old_m.at(ial[5])+jac*(b3*0.0111+0.0056*b1+0.0111*b2)*v_old_m.at(ial[5])+jac*(c3*0.0111+0.0056*c1+0.0111*c2)*w_old_m.at(ial[5])
                               +jac*(a3*0.0056-0.0056*a1+1.0835e-18*a2-0.0056*a3)*u_old_m.at(ial[6])+jac*(b3*0.0056-0.0056*b1+1.0835e-18*b2-0.0056*b3)*v_old_m.at(ial[6])+jac*(c3*0.0056-0.0056*c1+1.0835e-18*c2-0.0056*c3)*w_old_m.at(ial[6])
                               +jac*(a3*0.0056-0.0111*a1-0.0111*a2-0.0056*a3)*u_old_m.at(ial[7])+jac*(b3*0.0056-0.0111*b1-0.0111*b2-0.0056*b3)*v_old_m.at(ial[7])+jac*(c3*0.0056-0.0111*c1-0.0111*c2-0.0056*c3)*w_old_m.at(ial[7])
                               +jac*(a3*0.0111+0.0111*a1+0.0111*a3)*u_old_m.at(ial[8])+jac*(b3*0.0111+0.0111*b1+0.0111*b3)*v_old_m.at(ial[8])+jac*(c3*0.0111+0.0111*c1+0.0111*c3)*w_old_m.at(ial[8])
                               +jac*(a3*0.0056+0.0111*a2+0.0056*a3)*u_old_m.at(ial[9])+jac*(b3*0.0056+0.0111*b2+0.0056*b3)*v_old_m.at(ial[9])+jac*(c3*0.0056+0.0111*c2+0.0056*c3)*w_old_m.at(ial[9]);
                               
    ske[2][0] = 0; (jac/dt)*0.0083+jac*(-(a1+a2+a3)*(-0.0028)-0.0028*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(-(b1+b2+b3)*(-0.0028)-0.0028*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(-(c1+c2+c3)*(-0.0028)-0.0028*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(-(a1+a2+a3)*(-0.0028)-a1*0.0028)*u_old_m.at(ial[1])+jac*(-(b1+b2+b3)*(-0.0028)-b1*0.0028)*v_old_m.at(ial[1])+jac*(-(c1+c2+c3)*(-0.0028)-c1*0.0028)*w_old_m.at(ial[1])
                               +jac*(-(a1+a2+a3)*(-2.8969e-07)+a2*0.0028)*u_old_m.at(ial[2])+jac*(-(b1+b2+b3)*(-2.8969e-07)+b2*0.0028)*v_old_m.at(ial[2])+jac*(-(c1+c2+c3)*(-2.8969e-07)+c2*0.0028)*w_old_m.at(ial[2])
                               +jac*(-(a1+a2+a3)*(-0.0028)-c3*0.0056)*u_old_m.at(ial[3])+jac*(-(b1+b2+b3)*(-0.0028)-c3*0.0056)*v_old_m.at(ial[3])+jac*(-(c1+c2+c3)*(-0.0028)-c3*0.0056)*w_old_m.at(ial[3])
                               +jac*(-(a1+a2+a3)*0.0056+0.0056*a1-0.0056*a2-0.0056*a3)*u_old_m.at(ial[4])+jac*(-(b1+b2+b3)*0.0056+0.0056*b1-0.0056*b2-0.0056*b3)*v_old_m.at(ial[4])+jac*(-(c1+c2+c3)*0.0056+0.0056*c1-0.0056*c2-0.0056*c3)*w_old_m.at(ial[4])
                               +jac*(-(a1+a2+a3)*0.0111+0.0111*a1+0.0056*a2)*u_old_m.at(ial[5])+jac*(-(b1+b2+b3)*0.0111+0.0111*b1+0.0056*b2)*v_old_m.at(ial[5])+jac*(-(c1+c2+c3)*0.0111+0.0111*c1+0.0056*c2)*w_old_m.at(ial[5])
                               +jac*(-(a1+a2+a3)*0.0111-0.0111*a1+5.6703e-19*a2-0.0111*a3)*u_old_m.at(ial[6])+jac*(-(b1+b2+b3)*0.0111-0.0111*b1+5.6703e-19*b2-0.0111*b3)*v_old_m.at(ial[6])+jac*(-(c1+c2+c3)*0.0111-0.0111*c1+5.6703e-19*c2-0.0111*c3)*w_old_m.at(ial[6])
                               +jac*(-(a1+a2+a3)*0.0056-0.0056*a1-0.0056*a2+0.0056*a3)*u_old_m.at(ial[7])+jac*(-(b1+b2+b3)*0.0056-0.0056*b1-0.0056*b2+0.0056*b3)*v_old_m.at(ial[7])+jac*(-(c1+c2+c3)*0.0056-0.0056*c1-0.0056*c2+0.0056*c3)*w_old_m.at(ial[7])
                               +jac*(-(a1+a2+a3)*0.0056+0.0056*a1+0.0056*a3)*u_old_m.at(ial[8])+jac*(-(b1+b2+b3)*0.0056+0.0056*b1+0.0056*b3)*v_old_m.at(ial[8])+jac*(-(c1+c2+c3)*0.0056+0.0056*c1+0.0056*c3)*w_old_m.at(ial[8])
                               +jac*(-(a1+a2+a3)*0.0111+0.0056*a2+0.0111*a3)*u_old_m.at(ial[9])+jac*(-(b1+b2+b3)*0.0111+0.0056*b2+0.0111*b3)*v_old_m.at(ial[9])+jac*(-(c1+c2+c3)*0.0111+0.0056*c2+0.0111*c3)*w_old_m.at(ial[9]);
                               
    ske[2][1] = (jac/dt)*0.0083+jac*(a1*(-0.0028)+0.0028*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(b1*(-0.0028)+0.0028*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(c1*(-0.0028)+0.0028*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(a1*(-0.0028)+a1*0.0028)*u_old_m.at(ial[1])+jac*(b1*(-0.0028)+b1*0.0028)*v_old_m.at(ial[1])+jac*(c1*(-0.0028)+c1*0.0028)*w_old_m.at(ial[1])
                               +jac*(a1*(-2.8969e-07)+a2*0.0028)*u_old_m.at(ial[2])+jac*(b1*(-2.8969e-07)+b2*0.0028)*v_old_m.at(ial[2])+jac*(c1*(-2.8969e-07)+c2*0.0028)*w_old_m.at(ial[2])
                               +jac*(a1*(-0.0028)-a3*0.0028)*u_old_m.at(ial[3])+jac*(b1*(-0.0028)-c3*0.0028)*v_old_m.at(ial[3])+jac*(c1*(-0.0028)-c3*0.0028)*w_old_m.at(ial[3])
                               +jac*(a1*0.0056-0.0056*a1-0.0111*a2-0.0111*a3)*u_old_m.at(ial[4])+jac*(b1*0.0056-0.0056*b1-0.0111*b2-0.0111*b3)*v_old_m.at(ial[4])+jac*(c1*0.0056-0.0056*c1-0.0111*c2-0.0111*c3)*w_old_m.at(ial[4])
                               +jac*(a1*0.0111+0.0111*a1+0.0111*a2)*u_old_m.at(ial[5])+jac*(b1*0.0111+0.0111*b1+0.0111*b2)*v_old_m.at(ial[5])+jac*(c1*0.0111+0.0111*c1+0.0111*c2)*w_old_m.at(ial[5])
                               +jac*(a1*0.0111-0.0111*a1-0.0056*a2-0.0111*a3)*u_old_m.at(ial[6])+jac*(b1*0.0111-0.0111*b1-0.0056*b2-0.0111*b3)*v_old_m.at(ial[6])+jac*(c1*0.0111-0.0111*c1-0.0056*c2-0.0111*c3)*w_old_m.at(ial[6])
                               +jac*(a1*0.0056-0.0056*a1-0.0056*a2+1.0835e-18*a3)*u_old_m.at(ial[7])+jac*(b1*0.0056-0.0056*b1-0.0056*b2+1.0835e-18*b3)*v_old_m.at(ial[7])+jac*(c1*0.0056-0.0056*c1-0.0056*c2+1.0835e-18*c3)*w_old_m.at(ial[7])
                               +jac*(a1*0.0056+0.0056*a1+0.0111*a3)*u_old_m.at(ial[8])+jac*(b1*0.0056+0.0056*b1+0.0111*b3)*v_old_m.at(ial[8])+jac*(c1*0.0056+0.0056*c1+0.0111*c3)*w_old_m.at(ial[8])
                               +jac*(a1*0.0111+0.0056*a2+0.0111*a3)*u_old_m.at(ial[9])+jac*(b1*0.0111+0.0056*b2+0.0111*b3)*v_old_m.at(ial[9])+jac*(c1*0.0111+0.0056*c2+0.0111*c3)*w_old_m.at(ial[9]);
                               
    ske[2][2] = (jac/dt)*0.0167+jac*(a2*(-0.0028)+0.0056*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(b2*(-0.0028)+0.0056*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(c2*(-0.0028)+0.0056*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(a2*(-0.0028)-a1*0.0056)*u_old_m.at(ial[1])+jac*(b2*(-0.0028)-b1*0.0056)*v_old_m.at(ial[1])+jac*(c2*(-0.0028)-c1*0.0056)*w_old_m.at(ial[1])
                               +jac*(a2*(-2.8969e-07)+a2*0.0167)*u_old_m.at(ial[2])+jac*(b2*(-2.8969e-07)+b2*0.0167)*v_old_m.at(ial[2])+jac*(c2*(-2.8969e-07)+c2*0.0167)*w_old_m.at(ial[2])
                               +jac*(a2*(-0.0028)-a3*0.0056)*u_old_m.at(ial[3])+jac*(b2*(-0.0028)-c3*0.0056)*v_old_m.at(ial[3])+jac*(c2*(-0.0028)-c3*0.0056)*w_old_m.at(ial[3])
                               +jac*(a2*0.0056+1.5999e-18*a1-0.0111*a2-0.0111*a3)*u_old_m.at(ial[4])+jac*(b2*0.0056+1.5999e-18*b1-0.0111*b2-0.0111*b3)*v_old_m.at(ial[4])+jac*(c2*0.0056+1.5999e-18*c1-0.0111*c2-0.0111*c3)*w_old_m.at(ial[4])
                               +jac*(a2*0.0111+0.0333*a1+0.0111*a2)*u_old_m.at(ial[5])+jac*(b2*0.0111+0.0111*b1+0.0333*b2)*v_old_m.at(ial[5])+jac*(c2*0.0111+0.0333*c1+0.0111*c2)*w_old_m.at(ial[5])
                               +jac*(a2*0.0111-0.0333*a1-0.0222*a2-0.0333*a3)*u_old_m.at(ial[6])+jac*(b2*0.0111-0.0333*b1-0.0222*b2-0.0333*b3)*v_old_m.at(ial[6])+jac*(c2*0.0111-0.0333*c1-0.0222*c2-0.0333*c3)*w_old_m.at(ial[6])
                               +jac*(a2*0.0056-0.0111*a1-0.0111*a2+1.5999e-18*a3)*u_old_m.at(ial[7])+jac*(b2*0.0056-0.0111*b1-0.0111*b2+1.5999e-18*b3)*v_old_m.at(ial[7])+jac*(c2*0.0056-0.0111*c1-0.0111*c2+1.5999e-18*c3)*w_old_m.at(ial[7])
                               +jac*(a2*0.0056+0.0111*a1+0.0111*a3)*u_old_m.at(ial[8])+jac*(b2*0.0056+0.0111*b1+0.0111*b3)*v_old_m.at(ial[8])+jac*(c2*0.0056+0.0111*c1+0.0111*c3)*w_old_m.at(ial[8])
                               +jac*(a2*0.0111+0.0111*a2+0.0333*a3)*u_old_m.at(ial[9])+jac*(b2*0.0111+0.0111*b2+0.0333*b3)*v_old_m.at(ial[9])+jac*(c2*0.0111+0.0111*c2+0.0333*c3)*w_old_m.at(ial[9]);
                               
                               
    ske[2][3] = (jac/dt)*0.0083+jac*(a3*(-0.0028)+0.0028*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(b3*(-0.0028)+0.0028*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(c3*(-0.0028)+0.0028*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(a3*(-0.0028)-a1*0.0028)*u_old_m.at(ial[1])+jac*(b3*(-0.0028)-b1*0.0028)*v_old_m.at(ial[1])+jac*(c3*(-0.0028)-c1*0.0028)*w_old_m.at(ial[1])
                               +jac*(a3*(-2.8969e-07)+a2*0.0028)*u_old_m.at(ial[2])+jac*(b3*(-2.8969e-07)+b2*0.0028)*v_old_m.at(ial[2])+jac*(c3*(-2.8969e-07)+c2*0.0028)*w_old_m.at(ial[2])
                               +jac*(a3*(-0.0028)+a3*0.0028)*u_old_m.at(ial[3])+jac*(b3*(-0.0028)+c3*0.0028)*v_old_m.at(ial[3])+jac*(c3*(-0.0028)+c3*0.0028)*w_old_m.at(ial[3])
                               +jac*(a3*0.0056+1.0835e-18*a1-0.0056*a2-0.0056*a3)*u_old_m.at(ial[4])+jac*(b3*0.0056+1.0835e-18*b1-0.0056*b2-0.0056*b3)*v_old_m.at(ial[4])+jac*(c3*0.0056+1.0835e-18*c1-0.0056*c2-0.0056*c3)*w_old_m.at(ial[4])
                               +jac*(a3*0.0111+0.0111*a1+0.0111*a2)*u_old_m.at(ial[5])+jac*(b3*0.0111+0.0111*b1+0.0111*b2)*v_old_m.at(ial[5])+jac*(c3*0.0111+0.0111*c1+0.0111*c2)*w_old_m.at(ial[5])
                               +jac*(a3*0.0111-0.0111*a1-0.0056*a2-0.0111*a3)*u_old_m.at(ial[6])+jac*(b3*0.0111-0.0111*b1-0.0056*b2-0.0111*b3)*v_old_m.at(ial[6])+jac*(c3*0.0111-0.0111*c1-0.0056*c2-0.0111*c3)*w_old_m.at(ial[6])
                               +jac*(a3*0.0056-0.0111*a1-0.0111*a2-0.0056*a3)*u_old_m.at(ial[7])+jac*(b3*0.0056-0.0111*b1-0.0111*b2-0.0056*b3)*v_old_m.at(ial[7])+jac*(c3*0.0056-0.0111*c1-0.0111*c2-0.0056*c3)*w_old_m.at(ial[7])
                               +jac*(a3*0.0056+0.0111*a1+0.0111*a3)*u_old_m.at(ial[8])+jac*(b3*0.0056+0.0111*b1+0.0111*b3)*v_old_m.at(ial[8])+jac*(c3*0.0056+0.0111*c1+0.0056*c3)*w_old_m.at(ial[8])
                               +jac*(a3*0.0111+0.0111*a2+0.0111*a3)*u_old_m.at(ial[9])+jac*(b3*0.0111+0.0111*b2+0.0111*b3)*v_old_m.at(ial[9])+jac*(c3*0.0111+0.0111*c2+0.0111*c3)*w_old_m.at(ial[9]);
                               
    ske[3][0] = 0; (jac/dt)*0.0083+jac*(-(a1+a2+a3)*(-0.0028)-0.0028*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(-(b1+b2+b3)*(-0.0028)-0.0028*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(-(c1+c2+c3)*(-0.0028)-0.0028*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(-(a1+a2+a3)*(-0.0028)-a1*0.0028)*u_old_m.at(ial[1])+jac*(-(b1+b2+b3)*(-0.0028)-b1*0.0028)*v_old_m.at(ial[1])+jac*(-(c1+c2+c3)*(-0.0028)-c1*0.0028)*w_old_m.at(ial[1])
                               +jac*(-(a1+a2+a3)*(-0.0028)-a2*0.0028)*u_old_m.at(ial[2])+jac*(-(b1+b2+b3)*(-0.0028)-b2*0.0028)*v_old_m.at(ial[2])+jac*(-(c1+c2+c3)*(-0.0028)-c2*0.0028)*w_old_m.at(ial[2])
                               +jac*(-(a1+a2+a3)*(-2.8969e-07)+a3*0.0028)*u_old_m.at(ial[3])+jac*(-(b1+b2+b3)*(-2.8969e-07)+c3*0.0028)*v_old_m.at(ial[3])+jac*(-(c1+c2+c3)*(-2.8969e-07)+c3*0.0028)*w_old_m.at(ial[3])
                               +jac*(-(a1+a2+a3)*0.0056+0.0056*a1-0.0056*a2-0.0056*a3)*u_old_m.at(ial[4])+jac*(-(b1+b2+b3)*0.0056+0.0056*b1-0.0056*b2-0.0056*b3)*v_old_m.at(ial[4])+jac*(-(c1+c2+c3)*0.0056+0.0056*c1-0.0056*c2-0.0056*c3)*w_old_m.at(ial[4])
                               +jac*(-(a1+a2+a3)*0.0056+0.0056*a1+0.0056*a2)*u_old_m.at(ial[5])+jac*(-(b1+b2+b3)*0.0056+0.0056*b1+0.0056*b2)*v_old_m.at(ial[5])+jac*(-(c1+c2+c3)*0.0056+0.0056*c1+0.0056*c2)*w_old_m.at(ial[5])
                               +jac*(-(a1+a2+a3)*0.0056+0.0056*a1+0.0056*a2+0.0056*a3)*u_old_m.at(ial[6])+jac*(-(b1+b2+b3)*0.0056+0.0056*b1+0.0056*b2+0.0056*b3)*v_old_m.at(ial[6])+jac*(-(c1+c2+c3)*0.0056+0.0056*c1+0.0056*c2+0.0056*c3)*w_old_m.at(ial[6])
                               +jac*(-(a1+a2+a3)*0.0111-0.0111*a1-0.0111*a2+5.6703e-19*a3)*u_old_m.at(ial[7])+jac*(-(b1+b2+b3)*0.0111-0.0111*b1-0.0111*b2+5.6703e-19*b3)*v_old_m.at(ial[7])+jac*(-(c1+c2+c3)*0.0111-0.0111*c1-0.0111*c2+5.6703e-19*c3)*w_old_m.at(ial[7])
                               +jac*(-(a1+a2+a3)*0.0111+0.0111*a1+0.0056*a3)*u_old_m.at(ial[8])+jac*(-(b1+b2+b3)*0.0111+0.0111*b1+0.0056*b3)*v_old_m.at(ial[8])+jac*(-(c1+c2+c3)*0.0111+0.0111*c1+0.0056*c3)*w_old_m.at(ial[8])
                               +jac*(-(a1+a2+a3)*0.0111+0.0111*a2+0.0056*a3)*u_old_m.at(ial[9])+jac*(-(b1+b2+b3)*0.0111+0.0111*c2+0.0056*c3)*v_old_m.at(ial[9])+jac*(-(c1+c2+c3)*0.0111+0.0111*c2+0.0056*c3)*w_old_m.at(ial[9]);
                               
    ske[3][1] = (jac/dt)*0.0083+jac*(a1*(-0.0028)+0.0028*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(b1*(-0.0028)+0.0028*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(c1*(-0.0028)+0.0028*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(a1*(-0.0028)+a1*0.0028)*u_old_m.at(ial[1])+jac*(b1*(-0.0028)+b1*0.0028)*v_old_m.at(ial[1])+jac*(c1*(-0.0028)+c1*0.0028)*w_old_m.at(ial[1])
                               +jac*(a1*(-0.0028)-a2*0.0028)*u_old_m.at(ial[2])+jac*(b1*(-0.0028)-b2*0.0028)*v_old_m.at(ial[2])+jac*(c1*(-0.0028)-c2*0.0028)*w_old_m.at(ial[2])
                               +jac*(a1*(-2.8969e-07)+a3*0.0028)*u_old_m.at(ial[3])+jac*(b1*(-2.8969e-07)+c3*0.0028)*v_old_m.at(ial[3])+jac*(c1*(-2.8969e-07)+c3*0.0028)*w_old_m.at(ial[3])
                               +jac*(a1*0.0056-0.0056*a1-0.0111*a2-0.0111*a3)*u_old_m.at(ial[4])+jac*(b1*0.0056-0.0056*b1-0.0111*b2-0.0111*b3)*v_old_m.at(ial[4])+jac*(c1*0.0056-0.0056*c1-0.0111*c2-0.0111*c3)*w_old_m.at(ial[4])
                               +jac*(a1*0.0056+0.0056*a1+0.0111*a2)*u_old_m.at(ial[5])+jac*(b1*0.0056+0.0056*b1+0.0111*b2)*v_old_m.at(ial[5])+jac*(c1*0.0056+0.0056*c1+0.0111*c2)*w_old_m.at(ial[5])
                               +jac*(a1*0.0056-0.0056*a1+1.0835e-18*a2-0.0056*a3)*u_old_m.at(ial[6])+jac*(b1*0.0056-0.0056*b1+1.0835e-18*b2-0.0056*b3)*v_old_m.at(ial[6])+jac*(c1*0.0056-0.0056*c1+1.0835e-18*c2-0.0056*c3)*w_old_m.at(ial[6])
                               +jac*(a1*0.0111-0.0111*a1-0.0111*a2-0.0056*a3)*u_old_m.at(ial[7])+jac*(b1*0.0111-0.0111*b1-0.0111*b2-0.0056*b3)*v_old_m.at(ial[7])+jac*(c1*0.0111-0.0111*c1-0.0111*c2-0.0056*c3)*w_old_m.at(ial[7])
                               +jac*(a1*0.0111+0.0111*a1+0.0111*a3)*u_old_m.at(ial[8])+jac*(b1*0.0111+0.0111*b1+0.0111*b3)*v_old_m.at(ial[8])+jac*(c1*0.0111+0.0111*c1+0.0111*c3)*w_old_m.at(ial[8])
                               +jac*(a1*0.0111+0.0111*a2+0.0056*a3)*u_old_m.at(ial[9])+jac*(b1*0.0111+0.0111*b2+0.0056*b3)*v_old_m.at(ial[9])+jac*(c1*0.0111+0.0111*c2+0.0056*c3)*w_old_m.at(ial[9]);
                               
                               
    ske[3][2] = (jac/dt)*0.0083+jac*(a2*(-0.0028)+0.0056*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(b2*(-0.0028)+0.0056*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(c2*(-0.0028)+0.0056*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(a2*(-0.0028)-a1*0.0056)*u_old_m.at(ial[1])+jac*(b2*(-0.0028)-b1*0.0056)*v_old_m.at(ial[1])+jac*(c2*(-0.0028)-c1*0.0056)*w_old_m.at(ial[1])
                               +jac*(a2*(-0.0028)+a2*0.0167)*u_old_m.at(ial[2])+jac*(b2*(-0.0028)+b2*0.0167)*v_old_m.at(ial[2])+jac*(c2*(-0.0028)+c2*0.0167)*w_old_m.at(ial[2])
                               +jac*(a2*(-2.8969e-07)-a3*0.0056)*u_old_m.at(ial[3])+jac*(b2*(-2.8969e-07)-c3*0.0056)*v_old_m.at(ial[3])+jac*(c2*(-2.8969e-07)-c3*0.0056)*w_old_m.at(ial[3])
                               +jac*(a2*0.0056+1.5999e-18*a1-0.0111*a2-0.0111*a3)*u_old_m.at(ial[4])+jac*(b2*0.0056+1.5999e-18*b1-0.0111*b2-0.0111*b3)*v_old_m.at(ial[4])+jac*(c2*0.0056+1.5999e-18*c1-0.0111*c2-0.0111*c3)*w_old_m.at(ial[4])
                               +jac*(a2*0.0056+0.0333*a1+0.0111*a2)*u_old_m.at(ial[5])+jac*(b2*0.0056+0.0333*b1+0.0111*b2)*v_old_m.at(ial[5])+jac*(c2*0.0056+0.0333*c1+0.0111*c2)*w_old_m.at(ial[5])
                               +jac*(a2*0.0056-0.0333*a1-0.0222*a2-0.0333*a3)*u_old_m.at(ial[6])+jac*(b2*0.0056-0.0333*b1-0.0222*b2-0.0333*b3)*v_old_m.at(ial[6])+jac*(c2*0.0056-0.0333*c1-0.0222*c2-0.0333*c3)*w_old_m.at(ial[6])
                               +jac*(a2*0.0111-0.0111*a1-0.0111*a2+1.5999e-18*a3)*u_old_m.at(ial[7])+jac*(b2*0.0111-0.0111*b1-0.0111*b2+1.5999e-18*b3)*v_old_m.at(ial[7])+jac*(c2*0.0111-0.0111*c1-0.0111*c2+1.5999e-18*c3)*w_old_m.at(ial[7])
                               +jac*(a2*0.0111+0.0111*a1+0.0111*a3)*u_old_m.at(ial[8])+jac*(b2*0.0111+0.0111*b1+0.0111*b3)*v_old_m.at(ial[8])+jac*(c2*0.0111+0.0111*c1+0.0111*c3)*w_old_m.at(ial[8])
                               +jac*(a2*0.0111+0.0111*a2+0.0333*a3)*u_old_m.at(ial[9])+jac*(b2*0.0111+0.0111*b2+0.0333*b3)*v_old_m.at(ial[9])+jac*(c2*0.0111+0.0111*c2+0.0333*c3)*w_old_m.at(ial[9]);
                               
   ske[3][3] = (jac/dt)*0.0167+jac*(a3*(-0.0028)+0.0028*(a1+a2+a3))*u_old_m.at(ial[0])+jac*(b3*(-0.0028)+0.0028*(b1+b2+b3))*v_old_m.at(ial[0])+jac*(c3*(-0.0028)+0.0028*(c1+c2+c3))*w_old_m.at(ial[0])
                               +jac*(a3*(-0.0028)-a1*0.0028)*u_old_m.at(ial[1])+jac*(b3*(-0.0028)-b1*0.0028)*v_old_m.at(ial[1])+jac*(c3*(-0.0028)-c1*0.0028)*w_old_m.at(ial[1])
                               +jac*(a3*(-0.0028)+a2*0.0028)*u_old_m.at(ial[2])+jac*(b3*(-0.0028)+b2*0.0028)*v_old_m.at(ial[2])+jac*(c3*(-0.0028)+c2*0.0028)*w_old_m.at(ial[2])
                               +jac*(a3*(-2.8969e-07)+a3*0.0028)*u_old_m.at(ial[3])+jac*(b3*(-2.8969e-07)+c3*0.0028)*v_old_m.at(ial[3])+jac*(c3*(-2.8969e-07)+c3*0.0028)*w_old_m.at(ial[3])
                               +jac*(a3*0.0056+1.0835e-18*a1-0.0056*a2-0.0056*a3)*u_old_m.at(ial[4])+jac*(b3*0.0056+1.0835e-18*b1-0.0056*b2-0.0056*b3)*v_old_m.at(ial[4])+jac*(c3*0.0056+1.0835e-18*c1-0.0056*c2-0.0056*c3)*w_old_m.at(ial[4])
                               +jac*(a3*0.0056+0.0111*a1+0.0056*a2)*u_old_m.at(ial[5])+jac*(b3*0.0056+0.0111*b1+0.0056*b2)*v_old_m.at(ial[5])+jac*(c3*0.0056+0.0111*c1+0.0056*c2)*w_old_m.at(ial[5])
                               +jac*(a3*0.0056-0.0111*a1-0.0056*a2-0.0111*a3)*u_old_m.at(ial[6])+jac*(b3*0.0056-0.0111*b1-0.0056*b2-0.0111*b3)*v_old_m.at(ial[6])+jac*(c3*0.0056-0.0111*c1-0.0056*c2-0.0111*c3)*w_old_m.at(ial[6])
                               +jac*(a3*0.0111-0.0111*a1-0.0111*a2-0.0056*a3)*u_old_m.at(ial[7])+jac*(b3*0.0111-0.0111*b1-0.0111*b2-0.0056*b3)*v_old_m.at(ial[7])+jac*(c3*0.0111-0.0111*c1-0.0111*c2-0.0056*c3)*w_old_m.at(ial[7])
                               +jac*(a3*0.0111+0.0111*a1+0.0056*a3)*u_old_m.at(ial[8])+jac*(b3*0.0111+0.0111*b1+0.0056*b3)*v_old_m.at(ial[8])+jac*(c3*0.0111+0.0111*c1+0.0056*c3)*w_old_m.at(ial[8])
                               +jac*(a3*0.0111+0.0111*a2+0.0111*a3)*u_old_m.at(ial[9])+jac*(b3*0.0111+0.0111*b2+0.0111*b3)*v_old_m.at(ial[9])+jac*(c3*0.0111+0.0111*c2+0.0111*c3)*w_old_m.at(ial[9]);
     
    fe[0] = 0.01878132*(B0_0(a,a,a) +B0_0(b,a,a) +B0_0(a,b,a) +B0_0(a,a,b)) +0.01224884*(B0_0(c,c,c) +B0_0(d,c,c) +B0_0(c,d,c) +B0_0(c,c,d))+0.007091003*(B0_0(f,e,e) + B0_0(e,f,e) +B0_0(e,e,f) +B0_0(e,f,f)+ B0_0(f,e,f) +B0_0(f,f,e));
    fe[1] = 0.01878132*(B0_1(a,a,a) +B0_1(b,a,a) +B0_1(a,b,a) +B0_1(a,a,b)) +0.01224884*(B0_1(c,c,c) +B0_1(d,c,c) +B0_1(c,d,c) +B0_1(c,c,d))+0.007091003*(B0_1(f,e,e) + B0_1(e,f,e) +B0_1(e,e,f) +B0_1(e,f,f)+ B0_1(f,e,f) +B0_1(f,f,e));
    fe[2] = 0.01878132*(B0_2(a,a,a) +B0_2(b,a,a) +B0_2(a,b,a) +B0_2(a,a,b)) +0.01224884*(B0_2(c,c,c) +B0_2(d,c,c) +B0_2(c,d,c) +B0_2(c,c,d))+0.007091003*(B0_2(f,e,e) + B0_2(e,f,e) +B0_2(e,e,f) +B0_2(e,f,f)+ B0_2(f,e,f) +B0_2(f,f,e));
    fe[3] = 0.01878132*(B0_3(a,a,a) +B0_3(b,a,a) +B0_3(a,b,a) +B0_3(a,a,b)) +0.01224884*(B0_3(c,c,c) +B0_3(d,c,c) +B0_3(c,d,c) +B0_3(c,c,d))+0.007091003*(B0_3(f,e,e) + B0_3(e,f,e) +B0_3(e,e,f) +B0_3(e,f,f)+ B0_3(f,e,f) +B0_3(f,f,e));
}

//A01
typedef double (*A01) (double x, double y, double z);
void CalcElem_Navier_Stokes_A01(int const ial[10], double const xc[], double ske[4][10], double fe[4], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
                 
     const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
               
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
     
      A01 function01[4][10] = {{A01_00,A01_11,A01_02,A01_03,A01_04,A01_05,A01_06,A01_07,A01_08,A01_09},
		                       {A01_10,A01_11,A01_12,A01_13,A01_14,A01_15,A01_16,A01_17,A01_18,A01_19},
		                       {A01_20,A01_11,A01_22,A01_23,A01_24,A01_25,A01_26,A01_27,A01_28,A01_29},
		                       {A01_30,A01_11,A01_32,A01_33,A01_34,A01_35,A01_36,A01_37,A01_38,A01_39}};
     
             
      for(int i=0; i<=3; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function01[i][j](a,a,a) +function01[i][j](b,a,a) +function01[i][j](a,b,a) +function01[i][j](a,a,b)) +0.01224884*(function01[i][j](c,c,c) +function01[i][j](d,c,c) +function01[i][j](c,d,c) +function01[i][j](c,c,d))
			            +0.007091003*(function01[i][j](f,e,e) + function01[i][j](e,f,e) +function01[i][j](e,e,f) +function01[i][j](e,f,f)+ function01[i][j](f,e,f) +function01[i][j](f,f,e)));
		  }
	   }

    /*ske[0][0] = 0.01878132*(A01_00(a,a,a) +A01_00(b,a,a) +A01_00(a,b,a) +A01_00(a,a,b)) +0.01224884*(A01_00(c,c,c) +A01_00(d,c,c) +A01_00(c,d,c) +A01_00(c,c,d))+0.007091003*(A01_00(f,e,e) + A01_00(e,f,e) +A01_00(e,e,f) +A01_00(e,f,f)+ A01_00(f,e,f) +A01_00(f,f,e));
    ske[0][1] = 0.01878132*(A01_01(a,a,a) +A01_01(b,a,a) +A01_01(a,b,a) +A01_01(a,a,b)) +0.01224884*(A01_01(c,c,c) +A01_01(d,c,c) +A01_01(c,d,c) +A01_01(c,c,d))+0.007091003*(A01_01(f,e,e) + A01_01(e,f,e) +A01_01(e,e,f) +A01_01(e,f,f)+ A01_01(f,e,f) +A01_01(f,f,e));
    ske[0][2] = 0.01878132*(A01_02(a,a,a) +A01_02(b,a,a) +A01_02(a,b,a) +A01_02(a,a,b)) +0.01224884*(A01_02(c,c,c) +A01_02(d,c,c) +A01_02(c,d,c) +A01_02(c,c,d))+0.007091003*(A01_02(f,e,e) + A01_02(e,f,e) +A01_02(e,e,f) +A01_02(e,f,f)+ A01_02(f,e,f) +A01_02(f,f,e));
    ske[0][3] = 0.01878132*(A01_03(a,a,a) +A01_03(b,a,a) +A01_03(a,b,a) +A01_03(a,a,b)) +0.01224884*(A01_03(c,c,c) +A01_03(d,c,c) +A01_03(c,d,c) +A01_03(c,c,d))+0.007091003*(A01_03(f,e,e) + A01_03(e,f,e) +A01_03(e,e,f) +A01_03(e,f,f)+ A01_03(f,e,f) +A01_03(f,f,e));
    ske[0][4] = 0.01878132*(A01_04(a,a,a) +A01_04(b,a,a) +A01_04(a,b,a) +A01_04(a,a,b)) +0.01224884*(A01_04(c,c,c) +A01_04(d,c,c) +A01_04(c,d,c) +A01_04(c,c,d))+0.007091003*(A01_04(f,e,e) + A01_04(e,f,e) +A01_04(e,e,f) +A01_04(e,f,f)+ A01_04(f,e,f) +A01_04(f,f,e));
    ske[0][5] = 0.01878132*(A01_05(a,a,a) +A01_05(b,a,a) +A01_05(a,b,a) +A01_05(a,a,b)) +0.01224884*(A01_05(c,c,c) +A01_05(d,c,c) +A01_05(c,d,c) +A01_05(c,c,d))+0.007091003*(A01_05(f,e,e) + A01_05(e,f,e) +A01_05(e,e,f) +A01_05(e,f,f)+ A01_05(f,e,f) +A01_05(f,f,e));
    ske[0][6] = 0.01878132*(A01_06(a,a,a) +A01_06(b,a,a) +A01_06(a,b,a) +A01_06(a,a,b)) +0.01224884*(A01_06(c,c,c) +A01_06(d,c,c) +A01_06(c,d,c) +A01_06(c,c,d))+0.007091003*(A01_06(f,e,e) + A01_06(e,f,e) +A01_06(e,e,f) +A01_06(e,f,f)+ A01_06(f,e,f) +A01_06(f,f,e));
    ske[0][7] = 0.01878132*(A01_07(a,a,a) +A01_07(b,a,a) +A01_07(a,b,a) +A01_07(a,a,b)) +0.01224884*(A01_07(c,c,c) +A01_07(d,c,c) +A01_07(c,d,c) +A01_07(c,c,d))+0.007091003*(A01_07(f,e,e) + A01_07(e,f,e) +A01_07(e,e,f) +A01_07(e,f,f)+ A01_07(f,e,f) +A01_07(f,f,e));
    ske[0][8] = 0.01878132*(A01_08(a,a,a) +A01_08(b,a,a) +A01_08(a,b,a) +A01_08(a,a,b)) +0.01224884*(A01_08(c,c,c) +A01_08(d,c,c) +A01_08(c,d,c) +A01_08(c,c,d))+0.007091003*(A01_08(f,e,e) + A01_08(e,f,e) +A01_08(e,e,f) +A01_08(e,f,f)+ A01_08(f,e,f) +A01_08(f,f,e));
    ske[0][9] = 0.01878132*(A01_09(a,a,a) +A01_09(b,a,a) +A01_09(a,b,a) +A01_09(a,a,b)) +0.01224884*(A01_09(c,c,c) +A01_09(d,c,c) +A01_09(c,d,c) +A01_09(c,c,d))+0.007091003*(A01_09(f,e,e) + A01_09(e,f,e) +A01_09(e,e,f) +A01_09(e,f,f)+ A01_09(f,e,f) +A01_09(f,f,e));
    ske[1][0] = 0.01878132*(A01_10(a,a,a) +A01_10(b,a,a) +A01_10(a,b,a) +A01_10(a,a,b)) +0.01224884*(A01_10(c,c,c) +A01_10(d,c,c) +A01_10(c,d,c) +A01_10(c,c,d))+0.007091003*(A01_10(f,e,e) + A01_10(e,f,e) +A01_10(e,e,f) +A01_10(e,f,f)+ A01_10(f,e,f) +A01_10(f,f,e));
    ske[1][1] = 0.01878132*(A01_11(a,a,a) +A01_11(b,a,a) +A01_11(a,b,a) +A01_11(a,a,b)) +0.01224884*(A01_11(c,c,c) +A01_11(d,c,c) +A01_11(c,d,c) +A01_11(c,c,d))+0.007091003*(A01_11(f,e,e) + A01_11(e,f,e) +A01_11(e,e,f) +A01_11(e,f,f)+ A01_11(f,e,f) +A01_11(f,f,e));
    ske[1][2] = 0.01878132*(A01_12(a,a,a) +A01_12(b,a,a) +A01_12(a,b,a) +A01_12(a,a,b)) +0.01224884*(A01_12(c,c,c) +A01_12(d,c,c) +A01_12(c,d,c) +A01_12(c,c,d))+0.007091003*(A01_12(f,e,e) + A01_12(e,f,e) +A01_12(e,e,f) +A01_12(e,f,f)+ A01_12(f,e,f) +A01_12(f,f,e));
    ske[1][3] = 0.01878132*(A01_13(a,a,a) +A01_13(b,a,a) +A01_13(a,b,a) +A01_13(a,a,b)) +0.01224884*(A01_13(c,c,c) +A01_13(d,c,c) +A01_13(c,d,c) +A01_13(c,c,d))+0.007091003*(A01_13(f,e,e) + A01_13(e,f,e) +A01_13(e,e,f) +A01_13(e,f,f)+ A01_13(f,e,f) +A01_13(f,f,e));
    ske[1][4] = 0.01878132*(A01_14(a,a,a) +A01_14(b,a,a) +A01_14(a,b,a) +A01_14(a,a,b)) +0.01224884*(A01_14(c,c,c) +A01_14(d,c,c) +A01_14(c,d,c) +A01_14(c,c,d))+0.007091003*(A01_14(f,e,e) + A01_14(e,f,e) +A01_14(e,e,f) +A01_14(e,f,f)+ A01_14(f,e,f) +A01_14(f,f,e));
    ske[1][5] = 0.01878132*(A01_15(a,a,a) +A01_15(b,a,a) +A01_15(a,b,a) +A01_15(a,a,b)) +0.01224884*(A01_15(c,c,c) +A01_15(d,c,c) +A01_15(c,d,c) +A01_15(c,c,d))+0.007091003*(A01_15(f,e,e) + A01_15(e,f,e) +A01_15(e,e,f) +A01_15(e,f,f)+ A01_15(f,e,f) +A01_15(f,f,e));
    ske[1][6] = 0.01878132*(A01_16(a,a,a) +A01_16(b,a,a) +A01_16(a,b,a) +A01_16(a,a,b)) +0.01224884*(A01_16(c,c,c) +A01_16(d,c,c) +A01_16(c,d,c) +A01_16(c,c,d))+0.007091003*(A01_16(f,e,e) + A01_16(e,f,e) +A01_16(e,e,f) +A01_16(e,f,f)+ A01_16(f,e,f) +A01_16(f,f,e));
    ske[1][7] = 0.01878132*(A01_17(a,a,a) +A01_17(b,a,a) +A01_17(a,b,a) +A01_17(a,a,b)) +0.01224884*(A01_17(c,c,c) +A01_17(d,c,c) +A01_17(c,d,c) +A01_17(c,c,d))+0.007091003*(A01_17(f,e,e) + A01_17(e,f,e) +A01_17(e,e,f) +A01_17(e,f,f)+ A01_17(f,e,f) +A01_17(f,f,e));
    ske[1][8] = 0.01878132*(A01_18(a,a,a) +A01_18(b,a,a) +A01_18(a,b,a) +A01_18(a,a,b)) +0.01224884*(A01_18(c,c,c) +A01_18(d,c,c) +A01_18(c,d,c) +A01_18(c,c,d))+0.007091003*(A01_18(f,e,e) + A01_18(e,f,e) +A01_18(e,e,f) +A01_18(e,f,f)+ A01_18(f,e,f) +A01_18(f,f,e));
    ske[1][9] = 0.01878132*(A01_19(a,a,a) +A01_19(b,a,a) +A01_19(a,b,a) +A01_19(a,a,b)) +0.01224884*(A01_19(c,c,c) +A01_19(d,c,c) +A01_19(c,d,c) +A01_19(c,c,d))+0.007091003*(A01_19(f,e,e) + A01_19(e,f,e) +A01_19(e,e,f) +A01_19(e,f,f)+ A01_19(f,e,f) +A01_19(f,f,e));
    ske[2][0] = 0.01878132*(A01_20(a,a,a) +A01_20(b,a,a) +A01_20(a,b,a) +A01_20(a,a,b)) +0.01224884*(A01_20(c,c,c) +A01_20(d,c,c) +A01_20(c,d,c) +A01_20(c,c,d))+0.007091003*(A01_20(f,e,e) + A01_20(e,f,e) +A01_20(e,e,f) +A01_20(e,f,f)+ A01_20(f,e,f) +A01_20(f,f,e));
    ske[2][1] = 0.01878132*(A01_21(a,a,a) +A01_21(b,a,a) +A01_21(a,b,a) +A01_21(a,a,b)) +0.01224884*(A01_21(c,c,c) +A01_21(d,c,c) +A01_21(c,d,c) +A01_21(c,c,d))+0.007091003*(A01_21(f,e,e) + A01_21(e,f,e) +A01_21(e,e,f) +A01_21(e,f,f)+ A01_21(f,e,f) +A01_21(f,f,e));
    ske[2][2] = 0.01878132*(A01_22(a,a,a) +A01_22(b,a,a) +A01_22(a,b,a) +A01_22(a,a,b)) +0.01224884*(A01_22(c,c,c) +A01_22(d,c,c) +A01_22(c,d,c) +A01_22(c,c,d))+0.007091003*(A01_22(f,e,e) + A01_22(e,f,e) +A01_22(e,e,f) +A01_22(e,f,f)+ A01_22(f,e,f) +A01_22(f,f,e));
    ske[2][3] = 0.01878132*(A01_23(a,a,a) +A01_23(b,a,a) +A01_23(a,b,a) +A01_23(a,a,b)) +0.01224884*(A01_23(c,c,c) +A01_23(d,c,c) +A01_23(c,d,c) +A01_23(c,c,d))+0.007091003*(A01_23(f,e,e) + A01_23(e,f,e) +A01_23(e,e,f) +A01_23(e,f,f)+ A01_23(f,e,f) +A01_23(f,f,e));
    ske[2][4] = 0.01878132*(A01_24(a,a,a) +A01_24(b,a,a) +A01_24(a,b,a) +A01_24(a,a,b)) +0.01224884*(A01_24(c,c,c) +A01_24(d,c,c) +A01_24(c,d,c) +A01_24(c,c,d))+0.007091003*(A01_24(f,e,e) + A01_24(e,f,e) +A01_24(e,e,f) +A01_24(e,f,f)+ A01_24(f,e,f) +A01_24(f,f,e));
    ske[2][5] = 0.01878132*(A01_25(a,a,a) +A01_25(b,a,a) +A01_25(a,b,a) +A01_25(a,a,b)) +0.01224884*(A01_25(c,c,c) +A01_25(d,c,c) +A01_25(c,d,c) +A01_25(c,c,d))+0.007091003*(A01_25(f,e,e) + A01_25(e,f,e) +A01_25(e,e,f) +A01_25(e,f,f)+ A01_25(f,e,f) +A01_25(f,f,e));
    ske[2][6] = 0.01878132*(A01_26(a,a,a) +A01_26(b,a,a) +A01_26(a,b,a) +A01_26(a,a,b)) +0.01224884*(A01_26(c,c,c) +A01_26(d,c,c) +A01_26(c,d,c) +A01_26(c,c,d))+0.007091003*(A01_26(f,e,e) + A01_26(e,f,e) +A01_26(e,e,f) +A01_26(e,f,f)+ A01_26(f,e,f) +A01_26(f,f,e));
    ske[2][7] = 0.01878132*(A01_27(a,a,a) +A01_27(b,a,a) +A01_27(a,b,a) +A01_27(a,a,b)) +0.01224884*(A01_27(c,c,c) +A01_27(d,c,c) +A01_27(c,d,c) +A01_27(c,c,d))+0.007091003*(A01_27(f,e,e) + A01_27(e,f,e) +A01_27(e,e,f) +A01_27(e,f,f)+ A01_27(f,e,f) +A01_27(f,f,e));
    ske[2][8] = 0.01878132*(A01_28(a,a,a) +A01_28(b,a,a) +A01_28(a,b,a) +A01_28(a,a,b)) +0.01224884*(A01_28(c,c,c) +A01_28(d,c,c) +A01_28(c,d,c) +A01_28(c,c,d))+0.007091003*(A01_28(f,e,e) + A01_28(e,f,e) +A01_28(e,e,f) +A01_28(e,f,f)+ A01_28(f,e,f) +A01_28(f,f,e));
    ske[2][9] = 0.01878132*(A01_29(a,a,a) +A01_29(b,a,a) +A01_29(a,b,a) +A01_29(a,a,b)) +0.01224884*(A01_29(c,c,c) +A01_29(d,c,c) +A01_29(c,d,c) +A01_29(c,c,d))+0.007091003*(A01_29(f,e,e) + A01_29(e,f,e) +A01_29(e,e,f) +A01_29(e,f,f)+ A01_29(f,e,f) +A01_29(f,f,e));
    ske[3][0] = 0.01878132*(A01_30(a,a,a) +A01_30(b,a,a) +A01_30(a,b,a) +A01_30(a,a,b)) +0.01224884*(A01_30(c,c,c) +A01_30(d,c,c) +A01_30(c,d,c) +A01_30(c,c,d))+0.007091003*(A01_30(f,e,e) + A01_30(e,f,e) +A01_30(e,e,f) +A01_30(e,f,f)+ A01_30(f,e,f) +A01_30(f,f,e));
    ske[3][1] = 0.01878132*(A01_31(a,a,a) +A01_31(b,a,a) +A01_31(a,b,a) +A01_31(a,a,b)) +0.01224884*(A01_31(c,c,c) +A01_31(d,c,c) +A01_31(c,d,c) +A01_31(c,c,d))+0.007091003*(A01_31(f,e,e) + A01_31(e,f,e) +A01_31(e,e,f) +A01_31(e,f,f)+ A01_31(f,e,f) +A01_31(f,f,e));
    ske[3][2] = 0.01878132*(A01_32(a,a,a) +A01_32(b,a,a) +A01_32(a,b,a) +A01_32(a,a,b)) +0.01224884*(A01_32(c,c,c) +A01_32(d,c,c) +A01_32(c,d,c) +A01_32(c,c,d))+0.007091003*(A01_32(f,e,e) + A01_32(e,f,e) +A01_32(e,e,f) +A01_32(e,f,f)+ A01_32(f,e,f) +A01_32(f,f,e));
    ske[3][3] = 0.01878132*(A01_33(a,a,a) +A01_33(b,a,a) +A01_33(a,b,a) +A01_33(a,a,b)) +0.01224884*(A01_33(c,c,c) +A01_33(d,c,c) +A01_33(c,d,c) +A01_33(c,c,d))+0.007091003*(A01_33(f,e,e) + A01_33(e,f,e) +A01_33(e,e,f) +A01_33(e,f,f)+ A01_33(f,e,f) +A01_33(f,f,e));
    ske[3][4] = 0.01878132*(A01_34(a,a,a) +A01_34(b,a,a) +A01_34(a,b,a) +A01_34(a,a,b)) +0.01224884*(A01_34(c,c,c) +A01_34(d,c,c) +A01_34(c,d,c) +A01_34(c,c,d))+0.007091003*(A01_34(f,e,e) + A01_34(e,f,e) +A01_34(e,e,f) +A01_34(e,f,f)+ A01_34(f,e,f) +A01_34(f,f,e));
    ske[3][5] = 0.01878132*(A01_35(a,a,a) +A01_35(b,a,a) +A01_35(a,b,a) +A01_35(a,a,b)) +0.01224884*(A01_35(c,c,c) +A01_35(d,c,c) +A01_35(c,d,c) +A01_35(c,c,d))+0.007091003*(A01_35(f,e,e) + A01_35(e,f,e) +A01_35(e,e,f) +A01_35(e,f,f)+ A01_35(f,e,f) +A01_35(f,f,e));
    ske[3][6] = 0.01878132*(A01_36(a,a,a) +A01_36(b,a,a) +A01_36(a,b,a) +A01_36(a,a,b)) +0.01224884*(A01_36(c,c,c) +A01_36(d,c,c) +A01_36(c,d,c) +A01_36(c,c,d))+0.007091003*(A01_36(f,e,e) + A01_36(e,f,e) +A01_36(e,e,f) +A01_36(e,f,f)+ A01_36(f,e,f) +A01_36(f,f,e));
    ske[3][7] = 0.01878132*(A01_37(a,a,a) +A01_37(b,a,a) +A01_37(a,b,a) +A01_37(a,a,b)) +0.01224884*(A01_37(c,c,c) +A01_37(d,c,c) +A01_37(c,d,c) +A01_37(c,c,d))+0.007091003*(A01_37(f,e,e) + A01_37(e,f,e) +A01_37(e,e,f) +A01_37(e,f,f)+ A01_37(f,e,f) +A01_37(f,f,e));
    ske[3][8] = 0.01878132*(A01_38(a,a,a) +A01_38(b,a,a) +A01_38(a,b,a) +A01_38(a,a,b)) +0.01224884*(A01_38(c,c,c) +A01_38(d,c,c) +A01_38(c,d,c) +A01_38(c,c,d))+0.007091003*(A01_38(f,e,e) + A01_38(e,f,e) +A01_38(e,e,f) +A01_38(e,f,f)+ A01_38(f,e,f) +A01_38(f,f,e));
    ske[3][9] = 0.01878132*(A01_39(a,a,a) +A01_39(b,a,a) +A01_39(a,b,a) +A01_39(a,a,b)) +0.01224884*(A01_39(c,c,c) +A01_39(d,c,c) +A01_39(c,d,c) +A01_39(c,c,d))+0.007091003*(A01_39(f,e,e) + A01_39(e,f,e) +A01_39(e,e,f) +A01_39(e,f,f)+ A01_39(f,e,f) +A01_39(f,f,e));

    fe[0] = fe[1] = fe[2]= fe[3] =  0;*/
}
//A02
typedef double (*A02) (double x, double y, double z);
void CalcElem_Navier_Stokes_A02(int const ial[10], double const xc[], double ske[4][10], double fe[4], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
                 
     const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
     
     A02 function01[4][10] = { {A01_00,A01_11,A01_02,A01_03,A01_04,A01_05,A01_06,A01_07,A01_08,A01_09},
		                       {A01_10,A01_11,A01_12,A01_13,A01_14,A01_15,A01_16,A01_17,A01_18,A01_19},
		                       {A01_20,A01_11,A01_22,A01_23,A01_24,A01_25,A01_26,A01_27,A01_28,A01_29},
		                       {A01_30,A01_11,A01_32,A01_33,A01_34,A01_35,A01_36,A01_37,A01_38,A01_39}};
     
             
      for(int i=0; i<=3; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function01[i][j](a,a,a) +function01[i][j](b,a,a) +function01[i][j](a,b,a) +function01[i][j](a,a,b)) +0.01224884*(function01[i][j](c,c,c) +function01[i][j](d,c,c) +function01[i][j](c,d,c) +function01[i][j](c,c,d))
			            +0.007091003*(function01[i][j](f,e,e) + function01[i][j](e,f,e) +function01[i][j](e,e,f) +function01[i][j](e,f,f)+ function01[i][j](f,e,f) +function01[i][j](f,f,e)));
		  }
	   }

    /*ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[0][4] = 0;
    ske[0][5] = 0;
    ske[0][6] = 0;
    ske[0][7] = 0;
    ske[0][8] = 0;
    ske[0][9] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[1][4] = 0;
    ske[1][5] = 0;
    ske[1][6] = 0;
    ske[1][7] = 0;
    ske[1][8] = 0;
    ske[1][9] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[2][4] = 0;
    ske[2][5] = 0;
    ske[2][6] = 0;
    ske[2][7] = 0;
    ske[2][8] = 0;
    ske[2][9] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[3][4] = 0;
    ske[3][5] = 0;
    ske[3][6] = 0;
    ske[3][7] = 0;
    ske[3][8] = 0;
    ske[3][9] = 0;

    fe[0] = fe[1] = fe[2]= fe[3] =  0;*/
}

//A03
typedef double (*A03) (double x, double y, double z);
void CalcElem_Navier_Stokes_A03(int const ial[10], double const xc[], double ske[4][10], double fe[4], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
                 
    const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
    
    A03 function03[4][10] = {{A03_00,A03_11,A03_02,A03_03,A03_04,A03_05,A03_06,A03_07,A03_08,A03_09},
		                     {A03_10,A03_11,A03_12,A03_13,A03_14,A03_15,A03_16,A03_17,A03_18,A03_19},
		                     {A03_20,A03_11,A03_22,A03_23,A03_24,A03_25,A03_26,A03_27,A03_28,A03_29},
		                     {A03_30,A03_11,A03_32,A03_33,A03_34,A03_35,A03_36,A03_37,A03_38,A03_39}};
    
             
      for(int i=0; i<=3; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function03[i][j](a,a,a) +function03[i][j](b,a,a) +function03[i][j](a,b,a) +function03[i][j](a,a,b)) +0.01224884*(function03[i][j](c,c,c) +function03[i][j](d,c,c) +function03[i][j](c,d,c) +function03[i][j](c,c,d))
			            +0.007091003*(function03[i][j](f,e,e) + function03[i][j](e,f,e) +function03[i][j](e,e,f) +function03[i][j](e,f,f)+ function03[i][j](f,e,f) +function03[i][j](f,f,e)));
		  }
	   }

    /*ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[0][4] = 0;
    ske[0][5] = 0;
    ske[0][6] = 0;
    ske[0][7] = 0;
    ske[0][8] = 0;
    ske[0][9] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[1][4] = 0;
    ske[1][5] = 0;
    ske[1][6] = 0;
    ske[1][7] = 0;
    ske[1][8] = 0;
    ske[1][9] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[2][4] = 0;
    ske[2][5] = 0;
    ske[2][6] = 0;
    ske[2][7] = 0;
    ske[2][8] = 0;
    ske[2][9] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[3][4] = 0;
    ske[3][5] = 0;
    ske[3][6] = 0;
    ske[3][7] = 0;
    ske[3][8] = 0;
    ske[3][9] = 0;

    fe[0] = fe[1] = fe[2]= fe[3] =  0;*/
}
//A10
typedef double (*A10) (double x, double y, double z);
void CalcElem_Navier_Stokes_A10(int const ial[10], double const xc[], double ske[10][4], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
                 
     const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
     
     A03 function03[4][10] = {{A03_00,A03_11,A03_02,A03_03,A03_04,A03_05,A03_06,A03_07,A03_08,A03_09},
		                     {A03_10,A03_11,A03_12,A03_13,A03_14,A03_15,A03_16,A03_17,A03_18,A03_19},
		                     {A03_20,A03_11,A03_22,A03_23,A03_24,A03_25,A03_26,A03_27,A03_28,A03_29},
		                     {A03_30,A03_11,A03_32,A03_33,A03_34,A03_35,A03_36,A03_37,A03_38,A03_39}};
    
             
      for(int i=0; i<=3; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function03[i][j](a,a,a) +function03[i][j](b,a,a) +function03[i][j](a,b,a) +function03[i][j](a,a,b)) +0.01224884*(function03[i][j](c,c,c) +function03[i][j](d,c,c) +function03[i][j](c,d,c) +function03[i][j](c,c,d))
			            +0.007091003*(function03[i][j](f,e,e) + function03[i][j](e,f,e) +function03[i][j](e,e,f) +function03[i][j](e,f,f)+ function03[i][j](f,e,f) +function03[i][j](f,f,e)));
		  }
	   }

    ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[4][0] = 0;
    ske[4][1] = 0;
    ske[4][2] = 0;
    ske[4][3] = 0;
    ske[5][0] = 0;
    ske[5][1] = 0;
    ske[5][2] = 0;
    ske[5][3] = 0;
    ske[6][0] = 0;
    ske[6][1] = 0;
    ske[6][2] = 0;
    ske[6][3] = 0;
    ske[7][0] = 0;
    ske[7][1] = 0;
    ske[7][2] = 0;
    ske[7][3] = 0;
    ske[8][0] = 0;
    ske[8][1] = 0;
    ske[8][2] = 0;
    ske[8][3] = 0;
    ske[9][0] = 0;
    ske[9][1] = 0;
    ske[9][2] = 0;
    ske[9][3] = 0;

    fe[0] = fe[1] = fe[2] = fe[3] = fe[4] = fe[5] = fe[6] = fe[7] = fe[8] = fe[9] =  0;
}
//A11
typedef double (*A11) (double x, double y, double z);
void CalcElem_Navier_Stokes_A11(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
                 
    const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;

    ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[0][4] = 0;
    ske[0][5] = 0;
    ske[0][6] = 0;
    ske[0][7] = 0;
    ske[0][8] = 0;
    ske[0][9] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[1][4] = 0;
    ske[1][5] = 0;
    ske[1][6] = 0;
    ske[1][7] = 0;
    ske[1][8] = 0;
    ske[1][9] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[2][4] = 0;
    ske[2][5] = 0;
    ske[2][6] = 0;
    ske[2][7] = 0;
    ske[2][8] = 0;
    ske[2][9] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[3][4] = 0;
    ske[3][5] = 0;
    ske[3][6] = 0;
    ske[3][7] = 0;
    ske[3][8] = 0;
    ske[3][9] = 0;
    ske[4][0] = 0;
    ske[4][1] = 0;
    ske[4][2] = 0;
    ske[4][3] = 0;
    ske[4][4] = 0;
    ske[4][5] = 0;
    ske[4][6] = 0;
    ske[4][7] = 0;
    ske[4][8] = 0;
    ske[4][9] = 0;
    ske[5][0] = 0;
    ske[5][1] = 0;
    ske[5][2] = 0;
    ske[5][3] = 0;
    ske[5][4] = 0;
    ske[5][5] = 0;
    ske[5][6] = 0;
    ske[5][7] = 0;
    ske[5][8] = 0;
    ske[5][9] = 0;
    ske[6][0] = 0;
    ske[6][1] = 0;
    ske[6][2] = 0;
    ske[6][3] = 0;
    ske[6][4] = 0;
    ske[6][5] = 0;
    ske[6][6] = 0;
    ske[6][7] = 0;
    ske[6][8] = 0;
    ske[6][9] = 0;
    ske[7][0] = 0;
    ske[7][1] = 0;
    ske[7][2] = 0;
    ske[7][3] = 0;
    ske[7][4] = 0;
    ske[7][5] = 0;
    ske[7][6] = 0;
    ske[7][7] = 0;
    ske[7][8] = 0;
    ske[7][9] = 0;
    ske[8][0] = 0;
    ske[8][1] = 0;
    ske[8][2] = 0;
    ske[8][3] = 0;
    ske[8][4] = 0;
    ske[8][5] = 0;
    ske[8][6] = 0;
    ske[8][7] = 0;
    ske[8][8] = 0;
    ske[8][9] = 0;
    ske[9][0] = 0;
    ske[9][1] = 0;
    ske[9][2] = 0;
    ske[9][3] = 0;
    ske[9][4] = 0;
    ske[9][5] = 0;
    ske[9][6] = 0;
    ske[9][7] = 0;
    ske[9][8] = 0;
    ske[9][9] = 0;

    fe[0] = fe[1] = fe[2] = fe[3] = fe[4] = fe[5] = fe[6] = fe[7] = fe[8] = fe[9] =  0;
}
//A12
typedef double (*A12) (double x, double y, double z);
void CalcElem_Navier_Stokes_A12(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
     const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;

    ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[0][4] = 0;
    ske[0][5] = 0;
    ske[0][6] = 0;
    ske[0][7] = 0;
    ske[0][8] = 0;
    ske[0][9] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[1][4] = 0;
    ske[1][5] = 0;
    ske[1][6] = 0;
    ske[1][7] = 0;
    ske[1][8] = 0;
    ske[1][9] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[2][4] = 0;
    ske[2][5] = 0;
    ske[2][6] = 0;
    ske[2][7] = 0;
    ske[2][8] = 0;
    ske[2][9] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[3][4] = 0;
    ske[3][5] = 0;
    ske[3][6] = 0;
    ske[3][7] = 0;
    ske[3][8] = 0;
    ske[3][9] = 0;
    ske[4][0] = 0;
    ske[4][1] = 0;
    ske[4][2] = 0;
    ske[4][3] = 0;
    ske[4][4] = 0;
    ske[4][5] = 0;
    ske[4][6] = 0;
    ske[4][7] = 0;
    ske[4][8] = 0;
    ske[4][9] = 0;
    ske[5][0] = 0;
    ske[5][1] = 0;
    ske[5][2] = 0;
    ske[5][3] = 0;
    ske[5][4] = 0;
    ske[5][5] = 0;
    ske[5][6] = 0;
    ske[5][7] = 0;
    ske[5][8] = 0;
    ske[5][9] = 0;
    ske[6][0] = 0;
    ske[6][1] = 0;
    ske[6][2] = 0;
    ske[6][3] = 0;
    ske[6][4] = 0;
    ske[6][5] = 0;
    ske[6][6] = 0;
    ske[6][7] = 0;
    ske[6][8] = 0;
    ske[6][9] = 0;
    ske[7][0] = 0;
    ske[7][1] = 0;
    ske[7][2] = 0;
    ske[7][3] = 0;
    ske[7][4] = 0;
    ske[7][5] = 0;
    ske[7][6] = 0;
    ske[7][7] = 0;
    ske[7][8] = 0;
    ske[7][9] = 0;
    ske[8][0] = 0;
    ske[8][1] = 0;
    ske[8][2] = 0;
    ske[8][3] = 0;
    ske[8][4] = 0;
    ske[8][5] = 0;
    ske[8][6] = 0;
    ske[8][7] = 0;
    ske[8][8] = 0;
    ske[8][9] = 0;
    ske[9][0] = 0;
    ske[9][1] = 0;
    ske[9][2] = 0;
    ske[9][3] = 0;
    ske[9][4] = 0;
    ske[9][5] = 0;
    ske[9][6] = 0;
    ske[9][7] = 0;
    ske[9][8] = 0;
    ske[9][9] = 0;

    fe[0] = fe[1] = fe[2] = fe[3] = fe[4] = fe[5] = fe[6] = fe[7] = fe[8] = fe[9] =  0;
}
//A13
typedef double (*A13) (double x, double y, double z);
void CalcElem_Navier_Stokes_A13(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
                 
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;

    ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[0][4] = 0;
    ske[0][5] = 0;
    ske[0][6] = 0;
    ske[0][7] = 0;
    ske[0][8] = 0;
    ske[0][9] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[1][4] = 0;
    ske[1][5] = 0;
    ske[1][6] = 0;
    ske[1][7] = 0;
    ske[1][8] = 0;
    ske[1][9] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[2][4] = 0;
    ske[2][5] = 0;
    ske[2][6] = 0;
    ske[2][7] = 0;
    ske[2][8] = 0;
    ske[2][9] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[3][4] = 0;
    ske[3][5] = 0;
    ske[3][6] = 0;
    ske[3][7] = 0;
    ske[3][8] = 0;
    ske[3][9] = 0;
    ske[4][0] = 0;
    ske[4][1] = 0;
    ske[4][2] = 0;
    ske[4][3] = 0;
    ske[4][4] = 0;
    ske[4][5] = 0;
    ske[4][6] = 0;
    ske[4][7] = 0;
    ske[4][8] = 0;
    ske[4][9] = 0;
    ske[5][0] = 0;
    ske[5][1] = 0;
    ske[5][2] = 0;
    ske[5][3] = 0;
    ske[5][4] = 0;
    ske[5][5] = 0;
    ske[5][6] = 0;
    ske[5][7] = 0;
    ske[5][8] = 0;
    ske[5][9] = 0;
    ske[6][0] = 0;
    ske[6][1] = 0;
    ske[6][2] = 0;
    ske[6][3] = 0;
    ske[6][4] = 0;
    ske[6][5] = 0;
    ske[6][6] = 0;
    ske[6][7] = 0;
    ske[6][8] = 0;
    ske[6][9] = 0;
    ske[7][0] = 0;
    ske[7][1] = 0;
    ske[7][2] = 0;
    ske[7][3] = 0;
    ske[7][4] = 0;
    ske[7][5] = 0;
    ske[7][6] = 0;
    ske[7][7] = 0;
    ske[7][8] = 0;
    ske[7][9] = 0;
    ske[8][0] = 0;
    ske[8][1] = 0;
    ske[8][2] = 0;
    ske[8][3] = 0;
    ske[8][4] = 0;
    ske[8][5] = 0;
    ske[8][6] = 0;
    ske[8][7] = 0;
    ske[8][8] = 0;
    ske[8][9] = 0;
    ske[9][0] = 0;
    ske[9][1] = 0;
    ske[9][2] = 0;
    ske[9][3] = 0;
    ske[9][4] = 0;
    ske[9][5] = 0;
    ske[9][6] = 0;
    ske[9][7] = 0;
    ske[9][8] = 0;
    ske[9][9] = 0;

    fe[0] = fe[1] = fe[2] = fe[3] = fe[4] = fe[5] = fe[6] = fe[7] = fe[8] = fe[9] =  0;
}
//A20
typedef double (*A20) (double x, double y, double z);
void CalcElem_Navier_Stokes_A20(int const ial[10], double const xc[], double ske[10][4], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
     const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;

    ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[4][0] = 0;
    ske[4][1] = 0;
    ske[4][2] = 0;
    ske[4][3] = 0;
    ske[5][0] = 0;
    ske[5][1] = 0;
    ske[5][2] = 0;
    ske[5][3] = 0;
    ske[6][0] = 0;
    ske[6][1] = 0;
    ske[6][2] = 0;
    ske[6][3] = 0;
    ske[7][0] = 0;
    ske[7][1] = 0;
    ske[7][2] = 0;
    ske[7][3] = 0;
    ske[8][0] = 0;
    ske[8][1] = 0;
    ske[8][2] = 0;
    ske[8][3] = 0;
    ske[9][0] = 0;
    ske[9][1] = 0;
    ske[9][2] = 0;
    ske[9][3] = 0;

    fe[0] = fe[1] = fe[2] = fe[3] = fe[4] = fe[5] = fe[6] = fe[7] = fe[8] = fe[9] =  0;
}
//A21
typedef double (*A21) (double x, double y, double z);
void CalcElem_Navier_Stokes_A21(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
    const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;

    ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[0][4] = 0;
    ske[0][5] = 0;
    ske[0][6] = 0;
    ske[0][7] = 0;
    ske[0][8] = 0;
    ske[0][9] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[1][4] = 0;
    ske[1][5] = 0;
    ske[1][6] = 0;
    ske[1][7] = 0;
    ske[1][8] = 0;
    ske[1][9] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[2][4] = 0;
    ske[2][5] = 0;
    ske[2][6] = 0;
    ske[2][7] = 0;
    ske[2][8] = 0;
    ske[2][9] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[3][4] = 0;
    ske[3][5] = 0;
    ske[3][6] = 0;
    ske[3][7] = 0;
    ske[3][8] = 0;
    ske[3][9] = 0;
    ske[4][0] = 0;
    ske[4][1] = 0;
    ske[4][2] = 0;
    ske[4][3] = 0;
    ske[4][4] = 0;
    ske[4][5] = 0;
    ske[4][6] = 0;
    ske[4][7] = 0;
    ske[4][8] = 0;
    ske[4][9] = 0;
    ske[5][0] = 0;
    ske[5][1] = 0;
    ske[5][2] = 0;
    ske[5][3] = 0;
    ske[5][4] = 0;
    ske[5][5] = 0;
    ske[5][6] = 0;
    ske[5][7] = 0;
    ske[5][8] = 0;
    ske[5][9] = 0;
    ske[6][0] = 0;
    ske[6][1] = 0;
    ske[6][2] = 0;
    ske[6][3] = 0;
    ske[6][4] = 0;
    ske[6][5] = 0;
    ske[6][6] = 0;
    ske[6][7] = 0;
    ske[6][8] = 0;
    ske[6][9] = 0;
    ske[7][0] = 0;
    ske[7][1] = 0;
    ske[7][2] = 0;
    ske[7][3] = 0;
    ske[7][4] = 0;
    ske[7][5] = 0;
    ske[7][6] = 0;
    ske[7][7] = 0;
    ske[7][8] = 0;
    ske[7][9] = 0;
    ske[8][0] = 0;
    ske[8][1] = 0;
    ske[8][2] = 0;
    ske[8][3] = 0;
    ske[8][4] = 0;
    ske[8][5] = 0;
    ske[8][6] = 0;
    ske[8][7] = 0;
    ske[8][8] = 0;
    ske[8][9] = 0;
    ske[9][0] = 0;
    ske[9][1] = 0;
    ske[9][2] = 0;
    ske[9][3] = 0;
    ske[9][4] = 0;
    ske[9][5] = 0;
    ske[9][6] = 0;
    ske[9][7] = 0;
    ske[9][8] = 0;
    ske[9][9] = 0;

    fe[0] = fe[1] = fe[2] = fe[3] = fe[4] = fe[5] = fe[6] = fe[7] = fe[8] = fe[9] =  0;
}
//A22
typedef double (*A22) (double x, double y, double z);
void CalcElem_Navier_Stokes_A22(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
    const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;

    ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[0][4] = 0;
    ske[0][5] = 0;
    ske[0][6] = 0;
    ske[0][7] = 0;
    ske[0][8] = 0;
    ske[0][9] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[1][4] = 0;
    ske[1][5] = 0;
    ske[1][6] = 0;
    ske[1][7] = 0;
    ske[1][8] = 0;
    ske[1][9] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[2][4] = 0;
    ske[2][5] = 0;
    ske[2][6] = 0;
    ske[2][7] = 0;
    ske[2][8] = 0;
    ske[2][9] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[3][4] = 0;
    ske[3][5] = 0;
    ske[3][6] = 0;
    ske[3][7] = 0;
    ske[3][8] = 0;
    ske[3][9] = 0;
    ske[4][0] = 0;
    ske[4][1] = 0;
    ske[4][2] = 0;
    ske[4][3] = 0;
    ske[4][4] = 0;
    ske[4][5] = 0;
    ske[4][6] = 0;
    ske[4][7] = 0;
    ske[4][8] = 0;
    ske[4][9] = 0;
    ske[5][0] = 0;
    ske[5][1] = 0;
    ske[5][2] = 0;
    ske[5][3] = 0;
    ske[5][4] = 0;
    ske[5][5] = 0;
    ske[5][6] = 0;
    ske[5][7] = 0;
    ske[5][8] = 0;
    ske[5][9] = 0;
    ske[6][0] = 0;
    ske[6][1] = 0;
    ske[6][2] = 0;
    ske[6][3] = 0;
    ske[6][4] = 0;
    ske[6][5] = 0;
    ske[6][6] = 0;
    ske[6][7] = 0;
    ske[6][8] = 0;
    ske[6][9] = 0;
    ske[7][0] = 0;
    ske[7][1] = 0;
    ske[7][2] = 0;
    ske[7][3] = 0;
    ske[7][4] = 0;
    ske[7][5] = 0;
    ske[7][6] = 0;
    ske[7][7] = 0;
    ske[7][8] = 0;
    ske[7][9] = 0;
    ske[8][0] = 0;
    ske[8][1] = 0;
    ske[8][2] = 0;
    ske[8][3] = 0;
    ske[8][4] = 0;
    ske[8][5] = 0;
    ske[8][6] = 0;
    ske[8][7] = 0;
    ske[8][8] = 0;
    ske[8][9] = 0;
    ske[9][0] = 0;
    ske[9][1] = 0;
    ske[9][2] = 0;
    ske[9][3] = 0;
    ske[9][4] = 0;
    ske[9][5] = 0;
    ske[9][6] = 0;
    ske[9][7] = 0;
    ske[9][8] = 0;
    ske[9][9] = 0;

    fe[0] = fe[1] = fe[2] = fe[3] = fe[4] = fe[5] = fe[6] = fe[7] = fe[8] = fe[9] =  0;
}
//A23
typedef double (*A23) (double x, double y, double z);
void CalcElem_Navier_Stokes_A23(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
    const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;

    ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[0][4] = 0;
    ske[0][5] = 0;
    ske[0][6] = 0;
    ske[0][7] = 0;
    ske[0][8] = 0;
    ske[0][9] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[1][4] = 0;
    ske[1][5] = 0;
    ske[1][6] = 0;
    ske[1][7] = 0;
    ske[1][8] = 0;
    ske[1][9] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[2][4] = 0;
    ske[2][5] = 0;
    ske[2][6] = 0;
    ske[2][7] = 0;
    ske[2][8] = 0;
    ske[2][9] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[3][4] = 0;
    ske[3][5] = 0;
    ske[3][6] = 0;
    ske[3][7] = 0;
    ske[3][8] = 0;
    ske[3][9] = 0;
    ske[4][0] = 0;
    ske[4][1] = 0;
    ske[4][2] = 0;
    ske[4][3] = 0;
    ske[4][4] = 0;
    ske[4][5] = 0;
    ske[4][6] = 0;
    ske[4][7] = 0;
    ske[4][8] = 0;
    ske[4][9] = 0;
    ske[5][0] = 0;
    ske[5][1] = 0;
    ske[5][2] = 0;
    ske[5][3] = 0;
    ske[5][4] = 0;
    ske[5][5] = 0;
    ske[5][6] = 0;
    ske[5][7] = 0;
    ske[5][8] = 0;
    ske[5][9] = 0;
    ske[6][0] = 0;
    ske[6][1] = 0;
    ske[6][2] = 0;
    ske[6][3] = 0;
    ske[6][4] = 0;
    ske[6][5] = 0;
    ske[6][6] = 0;
    ske[6][7] = 0;
    ske[6][8] = 0;
    ske[6][9] = 0;
    ske[7][0] = 0;
    ske[7][1] = 0;
    ske[7][2] = 0;
    ske[7][3] = 0;
    ske[7][4] = 0;
    ske[7][5] = 0;
    ske[7][6] = 0;
    ske[7][7] = 0;
    ske[7][8] = 0;
    ske[7][9] = 0;
    ske[8][0] = 0;
    ske[8][1] = 0;
    ske[8][2] = 0;
    ske[8][3] = 0;
    ske[8][4] = 0;
    ske[8][5] = 0;
    ske[8][6] = 0;
    ske[8][7] = 0;
    ske[8][8] = 0;
    ske[8][9] = 0;
    ske[9][0] = 0;
    ske[9][1] = 0;
    ske[9][2] = 0;
    ske[9][3] = 0;
    ske[9][4] = 0;
    ske[9][5] = 0;
    ske[9][6] = 0;
    ske[9][7] = 0;
    ske[9][8] = 0;
    ske[9][9] = 0;

    fe[0] = fe[1] = fe[2] = fe[3] = fe[4] = fe[5] = fe[6] = fe[7] = fe[8] = fe[9] =  0;
}
//A30
typedef double (*A30) (double x, double y, double z);
void CalcElem_Navier_Stokes_A30(int const ial[10], double const xc[], double ske[10][4], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
    const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;

    ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[4][0] = 0;
    ske[4][1] = 0;
    ske[4][2] = 0;
    ske[4][3] = 0;
    ske[5][0] = 0;
    ske[5][1] = 0;
    ske[5][2] = 0;
    ske[5][3] = 0;
    ske[6][0] = 0;
    ske[6][1] = 0;
    ske[6][2] = 0;
    ske[6][3] = 0;
    ske[7][0] = 0;
    ske[7][1] = 0;
    ske[7][2] = 0;
    ske[7][3] = 0;
    ske[8][0] = 0;
    ske[8][1] = 0;
    ske[8][2] = 0;
    ske[8][3] = 0;
    ske[9][0] = 0;
    ske[9][1] = 0;
    ske[9][2] = 0;
    ske[9][3] = 0;

    fe[0] = fe[1] = fe[2] = fe[3] = fe[4] = fe[5] = fe[6] = fe[7] = fe[8] = fe[9] =  0;
}
//A31
typedef double (*A31) (double x, double y, double z);
void CalcElem_Navier_Stokes_A31(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
     const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;

    ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[0][4] = 0;
    ske[0][5] = 0;
    ske[0][6] = 0;
    ske[0][7] = 0;
    ske[0][8] = 0;
    ske[0][9] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[1][4] = 0;
    ske[1][5] = 0;
    ske[1][6] = 0;
    ske[1][7] = 0;
    ske[1][8] = 0;
    ske[1][9] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[2][4] = 0;
    ske[2][5] = 0;
    ske[2][6] = 0;
    ske[2][7] = 0;
    ske[2][8] = 0;
    ske[2][9] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[3][4] = 0;
    ske[3][5] = 0;
    ske[3][6] = 0;
    ske[3][7] = 0;
    ske[3][8] = 0;
    ske[3][9] = 0;
    ske[4][0] = 0;
    ske[4][1] = 0;
    ske[4][2] = 0;
    ske[4][3] = 0;
    ske[4][4] = 0;
    ske[4][5] = 0;
    ske[4][6] = 0;
    ske[4][7] = 0;
    ske[4][8] = 0;
    ske[4][9] = 0;
    ske[5][0] = 0;
    ske[5][1] = 0;
    ske[5][2] = 0;
    ske[5][3] = 0;
    ske[5][4] = 0;
    ske[5][5] = 0;
    ske[5][6] = 0;
    ske[5][7] = 0;
    ske[5][8] = 0;
    ske[5][9] = 0;
    ske[6][0] = 0;
    ske[6][1] = 0;
    ske[6][2] = 0;
    ske[6][3] = 0;
    ske[6][4] = 0;
    ske[6][5] = 0;
    ske[6][6] = 0;
    ske[6][7] = 0;
    ske[6][8] = 0;
    ske[6][9] = 0;
    ske[7][0] = 0;
    ske[7][1] = 0;
    ske[7][2] = 0;
    ske[7][3] = 0;
    ske[7][4] = 0;
    ske[7][5] = 0;
    ske[7][6] = 0;
    ske[7][7] = 0;
    ske[7][8] = 0;
    ske[7][9] = 0;
    ske[8][0] = 0;
    ske[8][1] = 0;
    ske[8][2] = 0;
    ske[8][3] = 0;
    ske[8][4] = 0;
    ske[8][5] = 0;
    ske[8][6] = 0;
    ske[8][7] = 0;
    ske[8][8] = 0;
    ske[8][9] = 0;
    ske[9][0] = 0;
    ske[9][1] = 0;
    ske[9][2] = 0;
    ske[9][3] = 0;
    ske[9][4] = 0;
    ske[9][5] = 0;
    ske[9][6] = 0;
    ske[9][7] = 0;
    ske[9][8] = 0;
    ske[9][9] = 0;

    fe[0] = fe[1] = fe[2] = fe[3] = fe[4] = fe[5] = fe[6] = fe[7] = fe[8] = fe[9] =  0;
}
//A32
typedef double (*A32) (double x, double y, double z);
void CalcElem_Navier_Stokes_A32(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
    
    const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;

    ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[0][4] = 0;
    ske[0][5] = 0;
    ske[0][6] = 0;
    ske[0][7] = 0;
    ske[0][8] = 0;
    ske[0][9] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[1][4] = 0;
    ske[1][5] = 0;
    ske[1][6] = 0;
    ske[1][7] = 0;
    ske[1][8] = 0;
    ske[1][9] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[2][4] = 0;
    ske[2][5] = 0;
    ske[2][6] = 0;
    ske[2][7] = 0;
    ske[2][8] = 0;
    ske[2][9] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[3][4] = 0;
    ske[3][5] = 0;
    ske[3][6] = 0;
    ske[3][7] = 0;
    ske[3][8] = 0;
    ske[3][9] = 0;
    ske[4][0] = 0;
    ske[4][1] = 0;
    ske[4][2] = 0;
    ske[4][3] = 0;
    ske[4][4] = 0;
    ske[4][5] = 0;
    ske[4][6] = 0;
    ske[4][7] = 0;
    ske[4][8] = 0;
    ske[4][9] = 0;
    ske[5][0] = 0;
    ske[5][1] = 0;
    ske[5][2] = 0;
    ske[5][3] = 0;
    ske[5][4] = 0;
    ske[5][5] = 0;
    ske[5][6] = 0;
    ske[5][7] = 0;
    ske[5][8] = 0;
    ske[5][9] = 0;
    ske[6][0] = 0;
    ske[6][1] = 0;
    ske[6][2] = 0;
    ske[6][3] = 0;
    ske[6][4] = 0;
    ske[6][5] = 0;
    ske[6][6] = 0;
    ske[6][7] = 0;
    ske[6][8] = 0;
    ske[6][9] = 0;
    ske[7][0] = 0;
    ske[7][1] = 0;
    ske[7][2] = 0;
    ske[7][3] = 0;
    ske[7][4] = 0;
    ske[7][5] = 0;
    ske[7][6] = 0;
    ske[7][7] = 0;
    ske[7][8] = 0;
    ske[7][9] = 0;
    ske[8][0] = 0;
    ske[8][1] = 0;
    ske[8][2] = 0;
    ske[8][3] = 0;
    ske[8][4] = 0;
    ske[8][5] = 0;
    ske[8][6] = 0;
    ske[8][7] = 0;
    ske[8][8] = 0;
    ske[8][9] = 0;
    ske[9][0] = 0;
    ske[9][1] = 0;
    ske[9][2] = 0;
    ske[9][3] = 0;
    ske[9][4] = 0;
    ske[9][5] = 0;
    ske[9][6] = 0;
    ske[9][7] = 0;
    ske[9][8] = 0;
    ske[9][9] = 0;

    fe[0] = fe[1] = fe[2] = fe[3] = fe[4] = fe[5] = fe[6] = fe[7] = fe[8] = fe[9] =  0;
}
//A33
typedef double (*A33) (double x, double y, double z);
void CalcElem_Navier_Stokes_A33(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3], i5 = 3 * ial[5], i6 = 3 * ial[6],
     i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
     
     const double u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
                 
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;

    ske[0][0] = 0;
    ske[0][1] = 0;
    ske[0][2] = 0;
    ske[0][3] = 0;
    ske[0][4] = 0;
    ske[0][5] = 0;
    ske[0][6] = 0;
    ske[0][7] = 0;
    ske[0][8] = 0;
    ske[0][9] = 0;
    ske[1][0] = 0;
    ske[1][1] = 0;
    ske[1][2] = 0;
    ske[1][3] = 0;
    ske[1][4] = 0;
    ske[1][5] = 0;
    ske[1][6] = 0;
    ske[1][7] = 0;
    ske[1][8] = 0;
    ske[1][9] = 0;
    ske[2][0] = 0;
    ske[2][1] = 0;
    ske[2][2] = 0;
    ske[2][3] = 0;
    ske[2][4] = 0;
    ske[2][5] = 0;
    ske[2][6] = 0;
    ske[2][7] = 0;
    ske[2][8] = 0;
    ske[2][9] = 0;
    ske[3][0] = 0;
    ske[3][1] = 0;
    ske[3][2] = 0;
    ske[3][3] = 0;
    ske[3][4] = 0;
    ske[3][5] = 0;
    ske[3][6] = 0;
    ske[3][7] = 0;
    ske[3][8] = 0;
    ske[3][9] = 0;
    ske[4][0] = 0;
    ske[4][1] = 0;
    ske[4][2] = 0;
    ske[4][3] = 0;
    ske[4][4] = 0;
    ske[4][5] = 0;
    ske[4][6] = 0;
    ske[4][7] = 0;
    ske[4][8] = 0;
    ske[4][9] = 0;
    ske[5][0] = 0;
    ske[5][1] = 0;
    ske[5][2] = 0;
    ske[5][3] = 0;
    ske[5][4] = 0;
    ske[5][5] = 0;
    ske[5][6] = 0;
    ske[5][7] = 0;
    ske[5][8] = 0;
    ske[5][9] = 0;
    ske[6][0] = 0;
    ske[6][1] = 0;
    ske[6][2] = 0;
    ske[6][3] = 0;
    ske[6][4] = 0;
    ske[6][5] = 0;
    ske[6][6] = 0;
    ske[6][7] = 0;
    ske[6][8] = 0;
    ske[6][9] = 0;
    ske[7][0] = 0;
    ske[7][1] = 0;
    ske[7][2] = 0;
    ske[7][3] = 0;
    ske[7][4] = 0;
    ske[7][5] = 0;
    ske[7][6] = 0;
    ske[7][7] = 0;
    ske[7][8] = 0;
    ske[7][9] = 0;
    ske[8][0] = 0;
    ske[8][1] = 0;
    ske[8][2] = 0;
    ske[8][3] = 0;
    ske[8][4] = 0;
    ske[8][5] = 0;
    ske[8][6] = 0;
    ske[8][7] = 0;
    ske[8][8] = 0;
    ske[8][9] = 0;
    ske[9][0] = 0;
    ske[9][1] = 0;
    ske[9][2] = 0;
    ske[9][3] = 0;
    ske[9][4] = 0;
    ske[9][5] = 0;
    ske[9][6] = 0;
    ske[9][7] = 0;
    ske[9][8] = 0;
    ske[9][9] = 0;

    fe[0] = fe[1] = fe[2] = fe[3] = fe[4] = fe[5] = fe[6] = fe[7] = fe[8] = fe[9] =  0;
}

/*
void CalcElem_heat_equation_crank_nichelson(int const ial[4], double const xc[], double ske[4][4], double fe[4], 
const vector<double> &u_old, const double dt, const double t, const double c)
{
    double tau=dt/c;

    const int  i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];
    const double x1 = xc[i2 + 0] - xc[i1 + 0],  y1 = xc[i2 + 1] - xc[i1 + 1], z1 = xc[i2 + 2] - xc[i1 + 2],
                 x2 = xc[i3 + 0] - xc[i1 + 0],  y2 = xc[i3 + 1] - xc[i1 + 1], z2 = xc[i3 + 2] - xc[i1 + 2],
                 x3 = xc[i4 + 0] - xc[i1 + 0],  y3 = xc[i4 + 1] - xc[i1 + 1], z3 = xc[i4 + 2] - xc[i1 + 2],
                 x4 = xc[i3 + 0] - xc[i2 + 0],  y4 = xc[i3 + 1] - xc[i2 + 1], z4 = xc[i3 + 2] - xc[i2 + 2],
                 x5 = xc[i2 + 0] - xc[i4 + 0],  y5 = xc[i2 + 1] - xc[i4 + 1], z5 = xc[i2 + 2] - xc[i4 + 2],
                 x6 = xc[i4 + 0] - xc[i3 + 0],  y6 = xc[i4 + 1] - xc[i3 + 1], z6 = xc[i4 + 2] - xc[i3 + 2];
    const double jac = fabs(x3*(y1*z2-y2*z1)+y3*(x2*z1-x1*z2)+z3*(x1*y2-y1*x2));

    ske[0][0] = (1/(6*jac)) * ((x4*x4+y4*y4+z4*z4)*(x5*x5+y5*y5+z5*z5)-(x4*x5+y4*y5+z4*z5)*(x4*x5+y4*y5+z4*z5));
    ske[0][1] = (1/(6*jac)) * ((x2*x6+y2*y6+z2*z6)*(x4*x6+y4*y6+z4*z6)-(x2*x4+y2*y4+z2*z4)*(x6*x6+y6*y6+z6*z6));
    ske[0][2] = (1/(6*jac)) * ((x3*x4+y3*y4+z3*z4)*(x5*x5+y5*y5+z5*z5)-(x3*x5+y3*y5+z3*z5)*(x4*x5+y4*y5+z4*z5));
    ske[0][3] = (1/(6*jac)) * ((x1*x6+y1*y6+z1*z6)*(x4*x4+y4*y4+z4*z4)-(x1*x4+y1*y4+z1*z4)*(x4*x6+y4*y6+z4*z6));
    ske[1][0] = ske[0][1];
    ske[1][1] = (1/(6*jac)) * ((x2*x2+y2*y2+z2*z2)*(x3*x3+y3*y3+z3*z3)-(x2*x3+y2*y3+z2*z3)*(x2*x3+y2*y3+z2*z3));
    ske[1][2] = (1/(6*jac)) * ((x1*x6+y1*y6+z1*z6)*(x3*x3+y3*y3+z3*z3)-(x1*x3+y1*y3+z1*z3)*(x3*x6+y3*y6+z3*z6));
    ske[1][3] = (1/(6*jac)) * ((x1*x3+y1*y3+z1*z3)*(x4*x6+y4*y6+z4*z6)-(x1*x6+y1*y6+z1*z6)*(x3*x4+y3*y4+z3*z4));
    ske[2][0] = ske[0][2];
    ske[2][1] = ske[1][2];
    ske[2][2] = (1/(6*jac)) * ((x1*x1+y1*y1+z1*z1)*(x3*x3+y3*y3+z3*z3)-(x1*x3+y1*y3+z1*z3)*(x1*x3+y1*y3+z1*z3));
    ske[2][3] = (1/(6*jac)) * ((x1*x3+y1*y3+z1*z3)*(x1*x4+y1*y4+z1*z4)-(x1*x1+y1*y1+z1*z1)*(x3*x4+y3*y4+z3*z4));
    ske[3][0] = ske[0][3];
    ske[3][1] = ske[1][3];
    ske[3][2] = ske[2][3];
    ske[3][3] = (1/(6*jac)) * ((x1*x1+y1*y1+z1*z1)*(x4*x4+y4*y4+z4*z4)-(x1*x4+y1*y4+z1*z4)*(x1*x4+y1*y4+z1*z4));
    

    const double xm    = (xc[i1 + 0] + xc[i2 + 0] + xc[i3 + 0]+ xc[i4 + 0]) / 4.0,
                 ym    = (xc[i1 + 1] + xc[i2 + 1] + xc[i3 + 1]+ xc[i4 + 1]) / 4.0,
                 zm    = (xc[i1 + 2] + xc[i2 + 2] + xc[i3 + 2]+ xc[i4 + 2]) / 4.0;
    //fe[0] = fe[1] = fe[2] = 0.5 * jac * FunctF(xm, ym) / 3.0;
    fe[0] = fe[1] = fe[2]= fe[3] =  0*(jac/6.0) * fNice(t+dt/2, xm, ym, zm) / 4.0; // jac*0.5/3 is integral of one hat function over the triangle (volume of pyramid)

    // add contributions from crank nichelson to right hand-side

    fe[0] += -ske[0][0]*tau/2*u_old.at(ial[0])-ske[0][1]*tau/2*u_old.at(ial[1])-ske[0][2]*tau/2*u_old.at(ial[2])-ske[0][3]*tau/2*u_old.at(ial[3]);
    fe[1] += -ske[1][0]*tau/2*u_old.at(ial[0])-ske[1][1]*tau/2*u_old.at(ial[1])-ske[1][2]*tau/2*u_old.at(ial[2])-ske[1][3]*tau/2*u_old.at(ial[3]);
    fe[2] += -ske[2][0]*tau/2*u_old.at(ial[0])-ske[2][1]*tau/2*u_old.at(ial[1])-ske[2][2]*tau/2*u_old.at(ial[2])-ske[2][3]*tau/2*u_old.at(ial[3]);
    fe[3] += -ske[3][0]*tau/2*u_old.at(ial[0])-ske[3][1]*tau/2*u_old.at(ial[1])-ske[3][2]*tau/2*u_old.at(ial[2])-ske[3][3]*tau/2*u_old.at(ial[3]);


    ske[0][0] = ske[0][0]*tau/2 + jac/12.0;
    ske[0][1] = ske[0][1]*tau/2 + jac/24.0;
    ske[0][2] = ske[0][2]*tau/2 + jac/24.0;
    ske[0][3] = ske[0][2]*tau/2 + jac/24.0;
    ske[1][0] = ske[1][0]*tau/2 + jac/24.0;
    ske[1][1] = ske[1][1]*tau/2 + jac/12.0;
    ske[1][2] = ske[1][2]*tau/2 + jac/24.0;
    ske[1][3] = ske[1][2]*tau/2 + jac/24.0;
    ske[2][0] = ske[2][0]*tau/2 + jac/24.0;
    ske[2][1] = ske[2][1]*tau/2 + jac/24.0;
    ske[2][2] = ske[2][2]*tau/2 + jac/12.0;
    ske[2][3] = ske[2][2]*tau/2 + jac/24.0;
    ske[3][0] = ske[2][0]*tau/2 + jac/24.0;
    ske[3][1] = ske[2][1]*tau/2 + jac/24.0;
    ske[3][2] = ske[2][2]*tau/2 + jac/24.0;
    ske[3][3] = ske[2][2]*tau/2 + jac/12.0;


    // add contributions from mass matrix/time derivative to rhs

    fe[0] += jac/12.0*u_old.at(ial[0]) + jac/24.0*u_old.at(ial[1]) + jac/24.0*u_old.at(ial[2])+ jac/24.0*u_old.at(ial[3]);
    fe[1] += jac/24.0*u_old.at(ial[0]) + jac/12.0*u_old.at(ial[1]) + jac/24.0*u_old.at(ial[2])+ jac/24.0*u_old.at(ial[3]);
    fe[2] += jac/24.0*u_old.at(ial[0]) + jac/24.0*u_old.at(ial[1]) + jac/12.0*u_old.at(ial[2])+ jac/24.0*u_old.at(ial[3]);
    fe[3] += jac/24.0*u_old.at(ial[0]) + jac/24.0*u_old.at(ial[1]) + jac/24.0*u_old.at(ial[2])+ jac/12.0*u_old.at(ial[3]);


}
*/
/*
// generalization to 3D for P2-P1 by Salman Ahmad

// generalization to 3D for P2 (quardatic) polynamial
void CalcElem_heat_equation_crank_nichelson_P2(int const ial[4], double const xc[], double ske[4][4], double fe[4], 
const vector<double> &u_old, const double dt, const double t, const double c)
{
    double tau=dt/c;

    const int  i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];
    const double x0 = xc[i1 + 0],  y0 = xc[i1 + 1], z0 = xc[i1 + 2],
                 x1 = xc[i2 + 0],  y1 = xc[i2 + 1], z1 = xc[i2 + 2],
                 x2 = xc[i3 + 0],  y2 = xc[i3 + 1], z2 = xc[i3 + 2],
                 x3 = xc[i4 + 0],  y3 = xc[i4 + 1], z3 = xc[i4 + 2],
                 x4 = (xc[i1 + 0]+xc[i2 + 0])/2,  y4 = (xc[i1 + 1]+xc[i2 + 1])/2, z4 = (xc[i1 + 2]+xc[i2 + 2])/2,
                 x5 = (xc[i2 + 0]+xc[i4 + 0])/2,  y5 = (xc[i2 + 1]+xc[i4 + 1])/2, z5 = (xc[i2 + 2]+xc[i4 + 2])/2,
                 x6 = (xc[i2 + 0]+xc[i3 + 0])/2,  y6 = (xc[i2 + 1]+xc[i3 + 1])/2, z6 = (xc[i2 + 2]+xc[i3 + 2])/2,
                 x7 = (xc[i3 + 0]+xc[i4 + 0])/2,  y7 = (xc[i3 + 1]+xc[i4 + 1])/2, z7 = (xc[i3 + 2]+xc[i4 + 2])/2,
                 x8 = (xc[i1 + 0]+xc[i3 + 0])/2,  y8 = (xc[i1 + 1]+xc[i3 + 1])/2, z8 = (xc[i1 + 2]+xc[i3 + 2])/2,
                 x9 = (xc[i1 + 0]+xc[i4 + 0])/2,  y9 = (xc[i1 + 1]+xc[i4 + 1])/2, z9 = (xc[i1 + 2]+xc[i4 + 2])/2;
   
    
    const double jac = fabs(x3*(y1*z2-y2*z1)+y3*(x2*z1-x1*z2)+z3*(x1*y2-y1*x2));
    
    const double  a0, a1, a2, a3, a4, a5, a6, a7, a8, a9,
                  b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, 
                  c0, c1, c2, c3, c4, c5, c6, c7, c8, c9,
                  d0, d1, d2, d3, d4, d5, d6, d7, d8, d9,
                  e0, e1, e2, e3, e4, e5, e6, e7, e8, e9,
                  f0, f1, f2, f3, f4, f5, f6, f7, f8, f9,
                  g0, g1, g2, g3, g4, g5, g6, g7, g8, g9,
                  h0, h1, h2, h3, h4, h5, h6, h7, h8, h9,
                  k0, k1, k2, k3, k4, k5, k6, k7, k8, k9,
                  r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
      a0 = ()/jac;
      a1 = -()/jac;
      a2 = ()/jac;
      a3 = -()/jac;
      a4 = ()/jac;
      a5 = -()/jac;
      a6 = ()/jac;
      a7 = -()/jac;
      a8 = ()/jac;
      a9 = -()/jac;
      
      b0 = -()/jac;
      b1 = ()/jac;
      b2 = -()/jac;
      b3 = ()/jac;
      b4 = -()/jac;
      b5 = ()/jac;
      b6 = -()/jac;
      b7 = ()/jac;
      b8 = -()/jac;
      b9 = ()/jac;
      
      c0 = ()/jac;
      c1 = -()/jac;
      c2 = ()/jac;
      c3 = -()/jac;
      c4 = ()/jac;
      c5 = -()/jac;
      c6 = ()/jac;
      c7 = -()/jac;
      c8 = ()/jac;
      c9 = -()/jac;
      
      d0 = -()/jac;
      d1 = ()/jac;
      d2 = -()/jac;
      d3 = ()/jac;
      d4 = -()/jac;
      d5 = ()/jac;
      d6 = -()/jac;
      d7 = ()/jac;
      d8 = -()/jac;
      d9 = ()/jac;
      
      e0 = ()/jac;
      e1 = -()/jac;
      e2 = ()/jac;
      e3 = -()/jac;
      e4 = ()/jac;
      e5 = -()/jac;
      e6 = ()/jac;
      e7 = -()/jac;
      e8 = ()/jac;
      e9 = -()/jac;
      
      f0 = -()/jac;
      f1 = ()/jac;
      f2 = -()/jac;
      f3 = ()/jac;
      f4 = -()/jac;
      f5 = ()/jac;
      f6 = -()/jac;
      f7 = ()/jac;
      f8 = -()/jac;
      f9 = ()/jac;
      
      g0 = ()/jac;
      g1 = -()/jac;
      g2 = ()/jac;
      g3 = -()/jac;
      g4 = ()/jac;
      g5 = -()/jac;
      g6 = ()/jac;
      g7 = -()/jac;
      g8 = ()/jac;
      g9 = -()/jac;
      
      h0 = -()/jac;
      h1 = ()/jac;
      h2 = -()/jac;
      h3 = ()/jac;
      h4 = -()/jac;
      h5 = ()/jac;
      h6 = -()/jac;
      h7 = ()/jac;
      h8 = -()/jac;
      h9 = ()/jac;
      
      k0 = ()/jac;
      k1 = -()/jac;
      k2 = ()/jac;
      k3 = -()/jac;
      k4 = ()/jac;
      k5 = -()/jac;
      k6 = ()/jac;
      k7 = -()/jac;
      k8 = ()/jac;
      k9 = -()/jac;
      
      r0 = -()/jac;
      r1 = ()/jac;
      r2 = -()/jac;
      r3 = ()/jac;
      r4 = -()/jac;
      r5 = ()/jac;
      r6 = -()/jac;
      r7 = ()/jac;
      r8 = -()/jac;
      r9 = ()/jac;
      
      
   
    
    

    ske[0][0] = (1/(6*jac)) * ((x4*x4+y4*y4+z4*z4)*(x5*x5+y5*y5+z5*z5)-(x4*x5+y4*y5+z4*z5)*(x4*x5+y4*y5+z4*z5));
    ske[0][1] = (1/(6*jac)) * ((x2*x6+y2*y6+z2*z6)*(x4*x6+y4*y6+z4*z6)-(x2*x4+y2*y4+z2*z4)*(x6*x6+y6*y6+z6*z6));
    ske[0][2] = (1/(6*jac)) * ((x3*x4+y3*y4+z3*z4)*(x5*x5+y5*y5+z5*z5)-(x3*x5+y3*y5+z3*z5)*(x4*x5+y4*y5+z4*z5));
    ske[0][3] = (1/(6*jac)) * ((x1*x6+y1*y6+z1*z6)*(x4*x4+y4*y4+z4*z4)-(x1*x4+y1*y4+z1*z4)*(x4*x6+y4*y6+z4*z6));
    ske[1][0] = ske[0][1];
    ske[1][1] = (1/(6*jac)) * ((x2*x2+y2*y2+z2*z2)*(x3*x3+y3*y3+z3*z3)-(x2*x3+y2*y3+z2*z3)*(x2*x3+y2*y3+z2*z3));
    ske[1][2] = (1/(6*jac)) * ((x1*x6+y1*y6+z1*z6)*(x3*x3+y3*y3+z3*z3)-(x1*x3+y1*y3+z1*z3)*(x3*x6+y3*y6+z3*z6));
    ske[1][3] = (1/(6*jac)) * ((x1*x3+y1*y3+z1*z3)*(x4*x6+y4*y6+z4*z6)-(x1*x6+y1*y6+z1*z6)*(x3*x4+y3*y4+z3*z4));
    ske[2][0] = ske[0][2];
    ske[2][1] = ske[1][2];
    ske[2][2] = (1/(6*jac)) * ((x1*x1+y1*y1+z1*z1)*(x3*x3+y3*y3+z3*z3)-(x1*x3+y1*y3+z1*z3)*(x1*x3+y1*y3+z1*z3));
    ske[2][3] = (1/(6*jac)) * ((x1*x3+y1*y3+z1*z3)*(x1*x4+y1*y4+z1*z4)-(x1*x1+y1*y1+z1*z1)*(x3*x4+y3*y4+z3*z4));
    ske[3][0] = ske[0][3];
    ske[3][1] = ske[1][3];
    ske[3][2] = ske[2][3];
    ske[3][3] = (1/(6*jac)) * ((x1*x1+y1*y1+z1*z1)*(x4*x4+y4*y4+z4*z4)-(x1*x4+y1*y4+z1*z4)*(x1*x4+y1*y4+z1*z4));
    

    const double xm    = (xc[i1 + 0] + xc[i2 + 0] + xc[i3 + 0]+ xc[i4 + 0]) / 4.0,
                 ym    = (xc[i1 + 1] + xc[i2 + 1] + xc[i3 + 1]+ xc[i4 + 1]) / 4.0,
                 zm    = (xc[i1 + 2] + xc[i2 + 2] + xc[i3 + 2]+ xc[i4 + 2]) / 4.0;
    //fe[0] = fe[1] = fe[2] = 0.5 * jac * FunctF(xm, ym) / 3.0;
    fe[0] = fe[1] = fe[2]= fe[3] =  0*(jac/6.0) * fNice(t+dt/2, xm, ym, zm) / 4.0; // jac*0.5/3 is integral of one hat function over the triangle (volume of pyramid)

    // add contributions from crank nichelson to right hand-side

    fe[0] += -ske[0][0]*tau/2*u_old.at(ial[0])-ske[0][1]*tau/2*u_old.at(ial[1])-ske[0][2]*tau/2*u_old.at(ial[2])-ske[0][3]*tau/2*u_old.at(ial[3]);
    fe[1] += -ske[1][0]*tau/2*u_old.at(ial[0])-ske[1][1]*tau/2*u_old.at(ial[1])-ske[1][2]*tau/2*u_old.at(ial[2])-ske[1][3]*tau/2*u_old.at(ial[3]);
    fe[2] += -ske[2][0]*tau/2*u_old.at(ial[0])-ske[2][1]*tau/2*u_old.at(ial[1])-ske[2][2]*tau/2*u_old.at(ial[2])-ske[2][3]*tau/2*u_old.at(ial[3]);
    fe[3] += -ske[3][0]*tau/2*u_old.at(ial[0])-ske[3][1]*tau/2*u_old.at(ial[1])-ske[3][2]*tau/2*u_old.at(ial[2])-ske[3][3]*tau/2*u_old.at(ial[3]);


    ske[0][0] = ske[0][0]*tau/2 + jac/12.0;
    ske[0][1] = ske[0][1]*tau/2 + jac/24.0;
    ske[0][2] = ske[0][2]*tau/2 + jac/24.0;
    ske[0][3] = ske[0][2]*tau/2 + jac/24.0;
    ske[1][0] = ske[1][0]*tau/2 + jac/24.0;
    ske[1][1] = ske[1][1]*tau/2 + jac/12.0;
    ske[1][2] = ske[1][2]*tau/2 + jac/24.0;
    ske[1][3] = ske[1][2]*tau/2 + jac/24.0;
    ske[2][0] = ske[2][0]*tau/2 + jac/24.0;
    ske[2][1] = ske[2][1]*tau/2 + jac/24.0;
    ske[2][2] = ske[2][2]*tau/2 + jac/12.0;
    ske[2][3] = ske[2][2]*tau/2 + jac/24.0;
    ske[3][0] = ske[2][0]*tau/2 + jac/24.0;
    ske[3][1] = ske[2][1]*tau/2 + jac/24.0;
    ske[3][2] = ske[2][2]*tau/2 + jac/24.0;
    ske[3][3] = ske[2][2]*tau/2 + jac/12.0;


    // add contributions from mass matrix/time derivative to rhs

    fe[0] += jac/12.0*u_old.at(ial[0]) + jac/24.0*u_old.at(ial[1]) + jac/24.0*u_old.at(ial[2])+ jac/24.0*u_old.at(ial[3]);
    fe[1] += jac/24.0*u_old.at(ial[0]) + jac/12.0*u_old.at(ial[1]) + jac/24.0*u_old.at(ial[2])+ jac/24.0*u_old.at(ial[3]);
    fe[2] += jac/24.0*u_old.at(ial[0]) + jac/24.0*u_old.at(ial[1]) + jac/12.0*u_old.at(ial[2])+ jac/24.0*u_old.at(ial[3]);
    fe[3] += jac/24.0*u_old.at(ial[0]) + jac/24.0*u_old.at(ial[1]) + jac/24.0*u_old.at(ial[2])+ jac/12.0*u_old.at(ial[3]);


}

*/

// generalization to 3D for P1 (linear) polynamial

/*

void CalcElem_heat_equation_crank_nichelson_P1(int const ial[4], double const xc[], double ske[4][4], double fe[4], 
const vector<double> &u_old, const double dt, const double t, const double c)
{
    double tau=dt/c;

    const int  i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];
    const double x0 = xc[i1 + 0],  y0 = xc[i1 + 1], z0 = xc[i1 + 2],
                 x1 = xc[i2 + 0],  y1 = xc[i2 + 1], z1 = xc[i2 + 2],
                 x2 = xc[i3 + 0],  y2 = xc[i3 + 1], z2 = xc[i3 + 2],
                 x3 = xc[i4 + 0],  y3 = xc[i4 + 1], z3 = xc[i4 + 2];
                 
    
    const double jac = fabs(x1*y3*z2 - x1*y2*z3 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1 + x1*y2*z3 - x1*y4*z2 - x2*y1*z4 + x2*y4*z1 + x4*y1*z2 
    - x4*y2*z1 - x1*y3*z4 + x1*y4*z3 + x3*y1*z4 - x3*y4*z1 - x4*y1*z3 + x4*y3*z1 + x2*y3*z4 - x2*y3*z3 - x3*y2*z4 + x3*y4*z2 + x4*y2*z3 - x4*y3*z2);
    
    const double  a0, a1, a2, a3,
                  b0, b1, b2, b3, 
                  c0, c1, c2, c3, 
                  d0, d1, d2, d3;
      a0 = x2*y3*z4 - x2*y4*z3 - x3*y2*z4 + x3*y4*z2 + x4*y2*z3 - x4*y3*z2;
      a1 = x1*y4*z3 - x1*y3*z4 + x3*y1*z4 - x3*y4*z1 - x4*y1*z3 + x4*y3*z1;
      a2 = x1*y2*z4 - x1*y4*z2 - x2*y1*z4 + x2*y4*z1 + x4*y1*z2 - x4*y2*z1;
      a3 = x1*y3*z2 - x1*y2*z3 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1;
      
      
      b0 = y3*z2 - y2*z3 + y2*z4 - y4*z2 - y3*z4 + y4*z3;
      b1 = y1*z3 - y3*z1 - y1*z4 + y4*z1 + y3*z4 - y4*z3;
      b2 = y2*z1 - y1*z2 + y1*z4 - y4*z1 - y2*z4 + y4*z2;
      b3 = y1*z2 - y2*z1 - y1*z3 + y3*z1 + y2*z3 - y3*z2;
      
      
      
      c0 = x2*z3 - x3*z2 - x2*z4 + x4*z2 + x3*z4 - x4*z3;
      c1 = x3*z1 - x1*z3 + x1*z4 - x4*z1 - x3*z4 + x4*z3;
      c2 = x1*z2 - x2*z1 - x1*z4 + x4*z1 + x2*z4 - x4*z2;
      c3 = x2*z1 - x1*z2 + x1*z3 - x3*z1 - x2*z3 + x3*z2;
      
      
      d0 = x3*y2 - x2*y3 + x2*y4 - x4*y2 - x3*y4 + x4*y3;
      d1 = x1*y3 - x3*y1 - x1*y4 + x4*y1 + x3*y4 - x4*y3;
      d2 = x2*y1 - x1*y2 + x1*y4 - x4*y1 - x2*y4 + x4*y2;
      d3 = x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2;
      
      
    
    

    ske[0][0] = (1/(6*jac)) * ((x4*x4+y4*y4+z4*z4)*(x5*x5+y5*y5+z5*z5)-(x4*x5+y4*y5+z4*z5)*(x4*x5+y4*y5+z4*z5));
    ske[0][1] = (1/(6*jac)) * ((x2*x6+y2*y6+z2*z6)*(x4*x6+y4*y6+z4*z6)-(x2*x4+y2*y4+z2*z4)*(x6*x6+y6*y6+z6*z6));
    ske[0][2] = (1/(6*jac)) * ((x3*x4+y3*y4+z3*z4)*(x5*x5+y5*y5+z5*z5)-(x3*x5+y3*y5+z3*z5)*(x4*x5+y4*y5+z4*z5));
    ske[0][3] = (1/(6*jac)) * ((x1*x6+y1*y6+z1*z6)*(x4*x4+y4*y4+z4*z4)-(x1*x4+y1*y4+z1*z4)*(x4*x6+y4*y6+z4*z6));
    ske[1][0] = ske[0][1];
    ske[1][1] = (1/(6*jac)) * ((x2*x2+y2*y2+z2*z2)*(x3*x3+y3*y3+z3*z3)-(x2*x3+y2*y3+z2*z3)*(x2*x3+y2*y3+z2*z3));
    ske[1][2] = (1/(6*jac)) * ((x1*x6+y1*y6+z1*z6)*(x3*x3+y3*y3+z3*z3)-(x1*x3+y1*y3+z1*z3)*(x3*x6+y3*y6+z3*z6));
    ske[1][3] = (1/(6*jac)) * ((x1*x3+y1*y3+z1*z3)*(x4*x6+y4*y6+z4*z6)-(x1*x6+y1*y6+z1*z6)*(x3*x4+y3*y4+z3*z4));
    ske[2][0] = ske[0][2];
    ske[2][1] = ske[1][2];
    ske[2][2] = (1/(6*jac)) * ((x1*x1+y1*y1+z1*z1)*(x3*x3+y3*y3+z3*z3)-(x1*x3+y1*y3+z1*z3)*(x1*x3+y1*y3+z1*z3));
    ske[2][3] = (1/(6*jac)) * ((x1*x3+y1*y3+z1*z3)*(x1*x4+y1*y4+z1*z4)-(x1*x1+y1*y1+z1*z1)*(x3*x4+y3*y4+z3*z4));
    ske[3][0] = ske[0][3];
    ske[3][1] = ske[1][3];
    ske[3][2] = ske[2][3];
    ske[3][3] = (1/(6*jac)) * ((x1*x1+y1*y1+z1*z1)*(x4*x4+y4*y4+z4*z4)-(x1*x4+y1*y4+z1*z4)*(x1*x4+y1*y4+z1*z4));
    

    const double xm    = (xc[i1 + 0] + xc[i2 + 0] + xc[i3 + 0]+ xc[i4 + 0]) / 4.0,
                 ym    = (xc[i1 + 1] + xc[i2 + 1] + xc[i3 + 1]+ xc[i4 + 1]) / 4.0,
                 zm    = (xc[i1 + 2] + xc[i2 + 2] + xc[i3 + 2]+ xc[i4 + 2]) / 4.0;
    //fe[0] = fe[1] = fe[2] = 0.5 * jac * FunctF(xm, ym) / 3.0;
    fe[0] = fe[1] = fe[2]= fe[3] =  0*(jac/6.0) * fNice(t+dt/2, xm, ym, zm) / 4.0; // jac*0.5/3 is integral of one hat function over the triangle (volume of pyramid)

    // add contributions from crank nichelson to right hand-side

    fe[0] += -ske[0][0]*tau/2*u_old.at(ial[0])-ske[0][1]*tau/2*u_old.at(ial[1])-ske[0][2]*tau/2*u_old.at(ial[2])-ske[0][3]*tau/2*u_old.at(ial[3]);
    fe[1] += -ske[1][0]*tau/2*u_old.at(ial[0])-ske[1][1]*tau/2*u_old.at(ial[1])-ske[1][2]*tau/2*u_old.at(ial[2])-ske[1][3]*tau/2*u_old.at(ial[3]);
    fe[2] += -ske[2][0]*tau/2*u_old.at(ial[0])-ske[2][1]*tau/2*u_old.at(ial[1])-ske[2][2]*tau/2*u_old.at(ial[2])-ske[2][3]*tau/2*u_old.at(ial[3]);
    fe[3] += -ske[3][0]*tau/2*u_old.at(ial[0])-ske[3][1]*tau/2*u_old.at(ial[1])-ske[3][2]*tau/2*u_old.at(ial[2])-ske[3][3]*tau/2*u_old.at(ial[3]);


    ske[0][0] = ske[0][0]*tau/2 + jac/12.0;
    ske[0][1] = ske[0][1]*tau/2 + jac/24.0;
    ske[0][2] = ske[0][2]*tau/2 + jac/24.0;
    ske[0][3] = ske[0][2]*tau/2 + jac/24.0;
    ske[1][0] = ske[1][0]*tau/2 + jac/24.0;
    ske[1][1] = ske[1][1]*tau/2 + jac/12.0;
    ske[1][2] = ske[1][2]*tau/2 + jac/24.0;
    ske[1][3] = ske[1][2]*tau/2 + jac/24.0;
    ske[2][0] = ske[2][0]*tau/2 + jac/24.0;
    ske[2][1] = ske[2][1]*tau/2 + jac/24.0;
    ske[2][2] = ske[2][2]*tau/2 + jac/12.0;
    ske[2][3] = ske[2][2]*tau/2 + jac/24.0;
    ske[3][0] = ske[2][0]*tau/2 + jac/24.0;
    ske[3][1] = ske[2][1]*tau/2 + jac/24.0;
    ske[3][2] = ske[2][2]*tau/2 + jac/24.0;
    ske[3][3] = ske[2][2]*tau/2 + jac/12.0;


    // add contributions from mass matrix/time derivative to rhs

    fe[0] += jac/12.0*u_old.at(ial[0]) + jac/24.0*u_old.at(ial[1]) + jac/24.0*u_old.at(ial[2])+ jac/24.0*u_old.at(ial[3]);
    fe[1] += jac/24.0*u_old.at(ial[0]) + jac/12.0*u_old.at(ial[1]) + jac/24.0*u_old.at(ial[2])+ jac/24.0*u_old.at(ial[3]);
    fe[2] += jac/24.0*u_old.at(ial[0]) + jac/24.0*u_old.at(ial[1]) + jac/12.0*u_old.at(ial[2])+ jac/24.0*u_old.at(ial[3]);
    fe[3] += jac/24.0*u_old.at(ial[0]) + jac/24.0*u_old.at(ial[1]) + jac/24.0*u_old.at(ial[2])+ jac/12.0*u_old.at(ial[3]);


}
*/
// ############################################################################################################################################################################################################


void AddElem(int const ial[4], double const ske[4][4], double const fe[4],
             int const id[], int const ik[], double sk[], double f[])
{
    for (int i = 0; i < 4; ++i)
    {
		for(int ki=0; ki<4; ++ki){
        const int ii  = 4*ial[i]+ki,           // row ii (global index)
                  id1 = id[ii],           // start and
                  id2 = id[ii + 1];       // end of row ii in matrix
        int ip  = id1;
        for (int j = 0; j < 4; ++j)       // no symmetry assumed
        {
			for(int kj=0; kj<4; ++kj){
            const int jj = 4*ial[j]+kj;
            bool not_found = true;
            do       // find entry jj (global index) in row ii
            {
                not_found = (ik[ip] != jj);
                ++ip;
            }
            while (not_found && ip < id2);

#ifndef NDEBUG                 // compiler option -DNDEBUG switches off the check
            if (not_found)     // no entry found !!
            {
                cout << "Error in AddElem: (" << ii << "," << jj << ") ["
                     << ial[0] << "," << ial[1] << "," << ial[2] << "," << ial[4] << "]\n";
                assert(!not_found);
            }
#endif

#pragma omp atomic
            sk[ip - 1] += ske[4*i+ki][4*j+kj];
        }
#pragma omp atomic
        f[ii] += fe[4*i+ki];
	}
	}
    }
}


// ----------------------------------------------------------------------------




// ############################################################################################################################################################################################################



CRS_Matrix::CRS_Matrix(Mesh const & mesh, int ndof_v)
 : _mesh(mesh), _nrows(0), _nnz(0), _id(0), _ik(0), _sk(0), _ndof_v(ndof_v)
{
    Derive_Matrix_Pattern();
    Skalar2VectorMatrix(ndof_v);
    return;
}

void CRS_Matrix::Derive_Matrix_Pattern()
{
    int const nelem(_mesh.Nelems());
    int const ndof_e(_mesh.NdofsElement());
    auto const &ia(_mesh.GetConnectivity());
//  Determine the number of matrix rows
    _nrows = *max_element(ia.cbegin(), ia.cbegin() + ndof_e * nelem);
    ++_nrows;                                 // node numberng: 0 ... nnode-1
    assert(*min_element(ia.cbegin(), ia.cbegin() + ndof_e * nelem) == 0); // numbering starts with 0 ?

//  Collect for each node those nodes it is connected to (multiple entries)
//  Detect the neighboring nodes
    vector< list<int> > cc(_nrows);             //  cc[i] is the  list of nodes a node i is connected to
    for (int i = 0; i < nelem; ++i)
    {
        int const idx = ndof_e * i;
        for (int k = 0; k < ndof_e; ++k)
        {
            list<int> &cck = cc.at(ia[idx + k]);
            cck.insert( cck.end(), ia.cbegin() + idx, ia.cbegin() + idx + ndof_e );
        }
    }
//  Delete the multiple entries
    _nnz = 0;
    for (auto &it : cc)
    {
        it.sort();
        it.unique();
        _nnz += it.size();
        // cout << it.size() << " :: "; copy(it->begin(),it->end(), ostream_iterator<int,char>(cout,"  ")); cout << endl;
    }

// CSR data allocation
    _id.resize(_nrows + 1);                  // Allocate memory for CSR row pointer
    _ik.resize(_nnz);                        // Allocate memory for CSR column index vector

//  copy CSR data
    _id[0] = 0;                                 // begin of first row
    for (size_t i = 0; i < cc.size(); ++i)
    {
        //cout << i << "   " << nid.at(i) << endl;;
        const list<int> &ci = cc.at(i);
        const auto nci = static_cast<int>(ci.size());
        _id[i + 1] = _id[i] + nci; // begin of next line
        copy(ci.begin(), ci.end(), _ik.begin() + _id[i] );
    }

    assert(_nnz == _id[_nrows]);
    _sk.resize(_nnz);                      // Allocate memory for CSR column index vector
    return;
}

void CRS_Matrix::Skalar2VectorMatrix(int ndof_v)
{
    this->Debug();
    cout << "\n########################\n";
    if (1 == ndof_v) return;
    assert(2 == ndof_v);

    auto old_id = _id;
    auto old_ik = _ik;

    _sk.resize(ndof_v * ndof_v * _sk.size(), -1.0);
    _id.resize(ndof_v * (_id.size() - 1) + 1);
    _ik.resize(ndof_v * ndof_v * _ik.size(), -7);

    _id[0] = 0;
    for (int kold = 0; kold < Nrows(); ++kold) {
        int nr = old_id[kold + 1] - old_id[kold];
        
        for (int ii=1; ii<=ndof_v; ++ii){
			_id[ndof_v * kold + ii] = _id[ndof_v * kold + ii-1] + ndof_v * nr;
		}
  
    }

    for (int newrow = 0; newrow < ndof_v * Nrows(); ++newrow ) {
        int oldrow = newrow / ndof_v;
        int idx = _id[newrow];
        //cout << " ("<< newrow<<")  " << idx;
        for (int oid = old_id[oldrow]; oid < old_id[oldrow + 1]; ++oid) {
            int oldcol = old_ik[oid];
            
            for (int ii=0; ii < ndof_v; ++ii){
				_ik[idx + ii] = ndof_v * oldcol+ii;
				}
         
            idx += ndof_v;
        }
    }

    _nrows =  _id.size() - 1;
    _nnz   =  _ik.size();

    return;
}

void CRS_Matrix::Debug() const
{
//  ID points to first entry of row
//  no symmetry assumed
    cout << "\nMatrix  (nnz = " << _id[_nrows] << ")\n";

    for (int row = 0; row < _nrows; ++row)
    {
        cout << "Row " << row << " : ";
        int const id1 = _id[row];
        int const id2 = _id[row + 1];
        for (int j = id1; j < id2; ++j)
        {
            cout.setf(ios::right, ios::adjustfield);
            cout << "[" << setw(2) << _ik[j] << "]  " << setw(4) << _sk[j] << "  ";
        }
        cout << endl;
    }
    return;
}

void CRS_Matrix::CalculateLaplace(vector<double> &f)
{
    assert(_mesh.NdofsElement() == 4);               // only for triangular, linear elements
    //cout << _nnz << " vs. " << _id[_nrows] << "  " << _nrows<< endl;
    assert(_nnz == _id[_nrows]);

    for (int k = 0; k < _nrows; ++k)
    {
        _sk[k] = 0.0;
    }
    for (int k = 0; k < _nrows; ++k)
    {
        f[k] = 0.0;
    }

    //double ske[3][3], fe[3];               // move inside loop (==> thread private)
    //  Loop over all elements
    auto const nelem = _mesh.Nelems();
    auto const &ia   = _mesh.GetConnectivity();
    auto const &xc   = _mesh.GetCoords();

#pragma omp parallel for
    for (int i = 0; i < nelem; ++i)
    {
        double ske[4][4], fe[4];             // OpenMP: Thread private
        CalcElem(ia.data()+4 * i, xc.data(), ske, fe);
        //AddElem(ia.data()+3 * i, ske, fe, _id.data(), _ik.data(), _sk.data(), f.data()); // GH: deprecated
        AddElem_3(ia.data()+4 * i, ske, fe, f);
    }

    //Debug();

    return;
}

/*
void FEM_Matrix::ApplyDirichletBC(std::vector<double> const &u, std::vector<double> &f)
{
    double const PENALTY = 1e6;
    auto const idx = _mesh.Index_DirichletNodes();
    int const nidx = idx.size();

    for (int row=0; row<nidx; ++row)
    {
        int const k = idx[row];
		int const id1 = fetch(k, k); // Find diagonal entry of row
		assert(id1 >= 0);
		_sk[id1] += PENALTY;		// matrix weighted scaling feasible
        f[k] += PENALTY * u[k];
    }

    return;
}
*/
/*
void FEM_Matrix::ApplyDirichletBC_Box(std::vector<double> const &u, std::vector<double> &f,
            double xl, double xh, double yl, double yh,double zl, double zh )
{
    double const PENALTY = 1e6;
    // auto const idx = _mesh.Index_DirichletNodes();
//    auto const idx = _mesh.Index_DirichletNodes_Box(xl, xh, yl, yh, zl, zh);
    auto lmesh=dynamic_cast<const Mesh_3d_4_matlab&>(_mesh);  // GH: Dirty hack
    auto const idx = lmesh.Index_DirichletNodes_Box(xl, xh, yl, yh, zl, zh);

    int const nidx = idx.size();

    for (int row=0; row<nidx; ++row)
    {
        int const k = idx[row];
		int const id1 = fetch(k, k); // Find diagonal entry of row
		assert(id1 >= 0);
		_sk[id1] += PENALTY;		// matrix weighted scaling feasible
        f[k] += PENALTY * u[k];
    }

    return;
}
*/
void CRS_Matrix::GetDiag(vector<double> &d) const
{
    assert( _nrows==static_cast<int>(d.size()) );

    for (int row = 0; row < _nrows; ++row)
    {
        const int ia = fetch(row, row); // Find diagonal entry of row
        assert(ia >= 0);
        d[row] = _sk[ia];
    }
    return;
}

bool CRS_Matrix::Compare2Old(int nnode, int const id[], int const ik[], double const sk[]) const
{
    bool bn = (nnode==_nrows);       // number of rows
    if (!bn)
    {
        cout << "#########   Error: " << "number of rows" << endl;
    }

    bool bz = (id[nnode]==_nnz);     // number of non zero elements
    if (!bz)
    {
        cout << "#########   Error: " << "number of non zero elements" << endl;
    }

    bool bd = equal(id,id+nnode+1,_id.cbegin());  // row starts
    if (!bd)
    {
        cout << "#########   Error: " << "row starts" << endl;
    }

    bool bk = equal(ik,ik+id[nnode],_ik.cbegin());  // column indices
    if (!bk)
    {
        cout << "#########   Error: " << "column indices" << endl;
    }

    bool bv = equal(sk,sk+id[nnode],_sk.cbegin());  // values
    if (!bv)
    {
        cout << "#########   Error: " << "values" << endl;
    }

    return bn && bz && bd && bk && bv;
}


void CRS_Matrix::Mult(vector<double> &w, vector<double> const &u) const
{
    assert( _nrows==static_cast<int>(w.size()) );
    assert( w.size()==u.size() );

#pragma omp parallel for
    for (int row = 0; row < _nrows; ++row)
    {
        double wi = 0.0;
        for (int ij = _id[row]; ij < _id[row + 1]; ++ij)
        {
            wi += _sk[ij] * u[ _ik[ij] ];
        }
        w[row] = wi;
    }
    return;
}

void CRS_Matrix::Defect(vector<double> &w,
                        vector<double> const &f, vector<double> const &u) const
{
    assert( _nrows==static_cast<int>(w.size()) );
    assert( w.size()==u.size() && u.size()==f.size() );

//#pragma omp parallel for default(none) shared(f,w,u)// class members are all shared
#pragma omp parallel for                          // works too because all data are shared
    for (int row = 0; row < _nrows; ++row)
    {
        double wi = f[row];
        for (int ij = _id[row]; ij < _id[row + 1]; ++ij)
        {
            wi -= _sk[ij] * u[ _ik[ij] ];
        }
        w[row] = wi;
    }
    return;
}

int CRS_Matrix::fetch(int const row, int const col) const
{
    int const id2 = _id[row + 1];    // end   and
    int       ip  = _id[row];        // start of recent row (global index)

    while (ip < id2 && _ik[ip] != col)  // find index col (global index)
    {
        ++ip;
    }
    if (ip >= id2)
    {
        ip = -1;
#ifndef NDEBUG                 // compiler option -DNDEBUG switches off the check
        cout << "No column  " << col << "  in row  " << row << endl;
        assert(ip >= id2);
#endif
    }
    return ip;
}


void CRS_Matrix::AddElem_3(int const ial[4], double const ske[4][4], double const fe[4], vector<double> & f)
{
    for (int i = 0; i < 4; ++i)
    {
        const int ii  = ial[i];           // row ii (global index)
        for (int j = 0; j < 4; ++j)       // no symmetry assumed
        {
            const int jj = ial[j];        // column jj (global index)
            int ip = fetch(ii,jj);        // find column entry jj in row ii
#ifndef NDEBUG                 // compiler option -DNDEBUG switches off the check
            if (ip<0)          // no entry found !!
            {
                cout << "Error in AddElem: (" << ii << "," << jj << ") ["
                     << ial[0] << "," << ial[1] << "," << ial[2] << "," << ial[3] << "]\n";
                assert(ip>=0);
            }
#endif
#pragma omp atomic
            _sk[ip] += ske[i][j];
        }
#pragma omp atomic
        f[ii] += fe[i];
    }
}
// ############################################################################################################################################################################################################


Matrix::Matrix(int const nrows, int const ncols)
    : _nrows(nrows), _ncols(ncols), _dd(0)
{}


Matrix::~Matrix()
{}




CRS_Matrix1::CRS_Matrix1()
    : Matrix(0, 0), _nnz(0), _id(0), _ik(0), _sk(0)
{}

CRS_Matrix1::CRS_Matrix1(const std::string &file) : Matrix(0, 0), _nnz(0), _id(0), _ik(0), _sk(0)
{
    readBinary(file);
    _nrows = static_cast<int>(size(_id) - 1);
    _ncols = _nrows;
}


CRS_Matrix1::~CRS_Matrix1()
{}

void CRS_Matrix1::Mult(vector<double> &w, vector<double> const &u) const
{
    assert( _ncols == static_cast<int>(u.size()) ); // compatibility of inner dimensions
    assert( _nrows == static_cast<int>(w.size()) ); // compatibility of outer dimensions

    #pragma omp parallel for
    for (int row = 0; row < _nrows; ++row) {
        double wi = 0.0;
        for (int ij = _id[row]; ij < _id[row + 1]; ++ij) {
            wi += _sk[ij] * u[ _ik[ij] ];
        }
        w[row] = wi;
    }
    return;
}

void CRS_Matrix1::Defect(vector<double> &w,
                        vector<double> const &f, vector<double> const &u) const
{
    assert( _ncols == static_cast<int>(u.size()) ); // compatibility of inner dimensions
    assert( _nrows == static_cast<int>(w.size()) ); // compatibility of outer dimensions
    assert( w.size() == f.size() );

    #pragma omp parallel for
    for (int row = 0; row < _nrows; ++row) {
        double wi = f[row];
        for (int ij = _id[row]; ij < _id[row + 1]; ++ij) {
            wi -= _sk[ij] * u[ _ik[ij] ];
        }
        w[row] = wi;
    }
    return;
}


void CRS_Matrix1::JacobiSmoother(std::vector<double> const &f, std::vector<double> &u,
                    std::vector<double> &r, int nsmooth, double const omega, bool zero) const
{
    // ToDO: ensure compatible dimensions
    assert(_ncols==_nrows);
    assert( _ncols == static_cast<int>(u.size()) ); // compatibility of inner dimensions
    assert( _nrows == static_cast<int>(r.size()) ); // compatibility of outer dimensions
    assert( r.size() == f.size() );
        
    int const nnodes = static_cast<int>(u.size());
    auto const &D = Matrix::GetDiag();        // accumulated diagonal of matrix @p SK.
    if (zero) {            // assumes initial solution is zero
        #pragma omp parallel for    
        for (int k = 0; k < nnodes; ++k) {
            // u := u + om*D^{-1}*f
            u[k] = omega*f[k] / D[k]; // MPI: distributed to accumulated vector needed
        }
        --nsmooth;                           // first smoothing sweep done
    }
    //cout << zero << endl;

    for (int ns = 1; ns <= nsmooth; ++ns) {
        //Defect(r, f, u);                  //  r := f - K*u    
        #pragma omp parallel for
        for (int row = 0; row < _nrows; ++row) {
            double wi = f[row];
            for (int ij = _id[row]; ij < _id[row + 1]; ++ij) {
                wi -= _sk[ij] * u[ _ik[ij] ];
            }
            r[row] = wi;
        }        
        #pragma omp parallel for    
        for (int k = 0; k < _nrows; ++k) {
            // u := u + om*D^{-1}*r
            u[k] = u[k] + omega * r[k] / D[k]; // MPI: distributed to accumulated vector needed
        }
    }

    return;
}

void CRS_Matrix1::GetDiag(vector<double> &d) const
{
    // be carefull when using a rectangular matrix
    int const nm = min(_nrows, _ncols);

    assert( nm == static_cast<int>(d.size()) ); // instead of stopping we could resize d and warn the user

    #pragma omp parallel for
    for (int row = 0; row < nm; ++row) {
        const int ia = fetch(row, row); // Find diagonal entry of row
        assert(ia >= 0);
        d[row] = _sk[ia];
    }
    cout << ">>>>> CRS_Matrix::GetDiag  <<<<<" << endl;
    return;
}

inline
int CRS_Matrix1::fetch(int const row, int const col) const
{
    int const id2 = _id[row + 1];    // end   and
    int       ip  = _id[row];        // start of recent row (global index)

    while (ip < id2 && _ik[ip] != col) { // find index col (global index)
        ++ip;
    }
    if (ip >= id2) {
        ip = -1;
#ifndef NDEBUG                 // compiler option -DNDEBUG switches off the check
        cout << "No column  " << col << "  in row  " << row << endl;
        assert(ip >= id2);
#endif
    }
    return ip;
}

void CRS_Matrix1::Debug() const
{
//  ID points to first entry of row
//  no symmetry assumed
    cout << "\nMatrix  (" << _nrows << " x " << _ncols << "  with  nnz = " << _id[_nrows] << ")\n";

    for (int row = 0; row < _nrows; ++row) {
        cout << "Row " << row << " : ";
        int const id1 = _id[row];
        int const id2 = _id[row + 1];
        for (int j = id1; j < id2; ++j) {
            cout.setf(ios::right, ios::adjustfield);
            cout << "[" << setw(2) << _ik[j] << "]  " << setw(4) << _sk[j] << "  ";
        }
        cout << endl;
    }
    return;
}



bool CRS_Matrix1::Compare2Old(int nnode, int const id[], int const ik[], double const sk[]) const
{
    bool bn = (nnode == _nrows);     // number of rows
    if (!bn) {
        cout << "#########   Error: " << "number of rows" << endl;
    }

    bool bz = (id[nnode] == _nnz);   // number of non zero elements
    if (!bz) {
        cout << "#########   Error: " << "number of non zero elements" << endl;
    }

    bool bd = equal(id, id + nnode + 1, _id.cbegin()); // row starts
    if (!bd) {
        cout << "#########   Error: " << "row starts" << endl;
    }

    bool bk = equal(ik, ik + id[nnode], _ik.cbegin()); // column indices
    if (!bk) {
        cout << "#########   Error: " << "column indices" << endl;
    }

    bool bv = equal(sk, sk + id[nnode], _sk.cbegin()); // values
    if (!bv) {
        cout << "#########   Error: " << "values" << endl;
    }

    return bn && bz && bd && bk && bv;
}


void CRS_Matrix1::writeBinary(const std::string &file)
{
    vector<int> cnt(size(_id) - 1);
    for (size_t k = 0; k < size(cnt); ++k) {
        cnt[k] = _id[k + 1] - _id[k];
    }
    //adjacent_difference( cbegin(_id)+1, cend(_id), cnt );
    write_binMatrix(file, cnt, _ik, _sk);
}

void CRS_Matrix1::readBinary(const std::string &file)
{
    vector<int> cnt;
    write_binMatrix(file, cnt, _ik, _sk);
    _id.resize(size(cnt) + 1);
    _id[0] = 0;
    for (size_t k = 0; k < size(cnt); ++k) {
        _id[k + 1] = _id[k] + cnt[k];
    }
    //partial_sum( cbegin(cnt), cend(cnt), begin(_id)+1 );
}



















// ############################################################################################################################################################################################################


FEM_Matrix::FEM_Matrix(Mesh const &mesh, int ndof_v)
    : CRS_Matrix1(), _mesh(mesh), _ndof_v(ndof_v)
{
    Derive_Matrix_Pattern();
    Skalar2VectorMatrix(ndof_v);
    return;
}

FEM_Matrix::~FEM_Matrix()
{}


void FEM_Matrix::Derive_Matrix_Pattern_fast()
{
    cout << "\n############   FEM_Matrix::Derive_Matrix_Pattern ";
    auto tstart = clock();
    int const nelem(_mesh.Nelems());
    int const ndof_e(_mesh.NdofsElement());
    auto const &ia(_mesh.GetConnectivity());
//  Determine the number of matrix rows
    _nrows = *max_element(ia.cbegin(), ia.cbegin() + ndof_e * nelem);
    ++_nrows;                                 // node numberng: 0 ... nnode-1
    assert(*min_element(ia.cbegin(), ia.cbegin() + ndof_e * nelem) == 0); // numbering starts with 0 ?

// CSR data allocation
    _id.resize(_nrows + 1);                  // Allocate memory for CSR row pointer
//##########################################################################
    auto const v2v = _mesh.Node2NodeGraph();
    _nnz = 0;                         // number of connections
    _id[0] = 0;                       // start of matrix row zero
    for (size_t v = 0; v < v2v.size(); ++v ) {
        _id[v + 1] = _id[v] + v2v[v].size();
        _nnz += v2v[v].size();
    }
    assert(_nnz == _id[_nrows]);
    _sk.resize(_nnz);                        // Allocate memory for CSR column index vector

// CSR data allocation
    _ik.resize(_nnz);                        // Allocate memory for CSR column index vector
// Copy column indices
    int kk = 0;
    for (size_t v = 0; v < v2v.size(); ++v ) {
        for (size_t vi = 0; vi < v2v[v].size(); ++vi) {
            _ik[kk] = v2v[v][vi];
            ++kk;
        }
    }
    _ncols = *max_element(_ik.cbegin(), _ik.cend());  // maximal column number
    ++_ncols;                                         // node numbering: 0 ... nnode-1
    //cout << _nrows << "  " << _ncols << endl;
    assert(_ncols == _nrows);

    double duration = static_cast<double>(clock() - tstart) / CLOCKS_PER_SEC;
    cout << "finished in  " <<  duration  << " sec.    ########\n";

    return;
}


void FEM_Matrix::Derive_Matrix_Pattern_slow()
{
    cout << "\n############   FEM_Matrix::Derive_Matrix_Pattern slow ";
    auto tstart = clock();
    int const nelem(_mesh.Nelems());
    int const ndof_e(_mesh.NdofsElement());
    auto const &ia(_mesh.GetConnectivity());
//  Determine the number of matrix rows
    _nrows = *max_element(ia.cbegin(), ia.cbegin() + ndof_e * nelem);
    ++_nrows;                                 // node numberng: 0 ... nnode-1
    assert(*min_element(ia.cbegin(), ia.cbegin() + ndof_e * nelem) == 0); // numbering starts with 0 ?

//  Collect for each node those nodes it is connected to (multiple entries)
//  Detect the neighboring nodes
    vector< list<int> > cc(_nrows);             //  cc[i] is the  list of nodes a node i is connected to
    for (int i = 0; i < nelem; ++i) {
        int const idx = ndof_e * i;
        for (int k = 0; k < ndof_e; ++k) {
            list<int> &cck = cc[ia[idx + k]];
            cck.insert( cck.end(), ia.cbegin() + idx, ia.cbegin() + idx + ndof_e );
        }
    }
//  Delete the multiple entries
    _nnz = 0;
    for (auto &it : cc) {
        it.sort();
        it.unique();
        _nnz += it.size();
        // cout << it.size() << " :: "; copy(it->begin(),it->end(), ostream_iterator<int,char>(cout,"  ")); cout << endl;
    }

// CSR data allocation
    _id.resize(_nrows + 1);                  // Allocate memory for CSR row pointer
    _ik.resize(_nnz);                        // Allocate memory for CSR column index vector

//  copy CSR data
    _id[0] = 0;                              // begin of first row
    for (size_t i = 0; i < cc.size(); ++i) {
        //cout << i << "   " << nid.at(i) << endl;;
        const list<int> &ci = cc[i];
        const auto nci = static_cast<int>(ci.size());
        _id[i + 1] = _id[i] + nci; // begin of next line
        copy(ci.begin(), ci.end(), _ik.begin() + _id[i] );
    }

    assert(_nnz == _id[_nrows]);
    _sk.resize(_nnz);                        // Allocate memory for CSR column index vector

    _ncols = *max_element(_ik.cbegin(), _ik.cend());  // maximal column number
    ++_ncols;                                 // node numbering: 0 ... nnode-1
    //cout << _nrows << "  " << _ncols << endl;
    assert(_ncols == _nrows);

    double duration = static_cast<double>(clock() - tstart) / CLOCKS_PER_SEC;
    cout << "finished in  " <<  duration  << " sec.    ########\n";

    return;
}

void FEM_Matrix::Skalar2VectorMatrix(int ndof_v)
{
    this->Debug();
    cout << "\n########################\n";
    if (1 == ndof_v) return;
    assert(4 == ndof_v);

    auto old_id = _id;
    auto old_ik = _ik;

    _sk.resize(ndof_v * ndof_v * _sk.size(), 0.0);
    _id.resize(ndof_v * (_id.size() - 1) + 1);
    _ik.resize(ndof_v * ndof_v * _ik.size(), -7);

    _id[0] = 0;
    for (int kold = 0; kold < Nrows(); ++kold) {
        int nr = old_id[kold + 1] - old_id[kold];
        
        for (int ii=1; ii<=ndof_v; ++ii){
			_id[ndof_v * kold + ii] = _id[ndof_v * kold + ii-1] + ndof_v * nr;
		}
  
    }

    for (int newrow = 0; newrow < ndof_v * Nrows(); ++newrow ) {
        int oldrow = newrow / ndof_v;
        int idx = _id[newrow];
        //cout << " ("<< newrow<<")  " << idx;
        for (int oid = old_id[oldrow]; oid < old_id[oldrow + 1]; ++oid) {
            int oldcol = old_ik[oid];
            
            for (int ii=0; ii < ndof_v; ++ii){
				_ik[idx + ii] = ndof_v * oldcol+ii;
				}
         
            idx += ndof_v;
        }
    }

    _nrows =  _id.size() - 1;
    _nnz   =  _ik.size();

    return;
}
/*
void CRS_Matrix::Skalar2VectorMatrix(int ndof_v)
{
    this->Debug();
    cout << "\n########################\n";
    if (1 == ndof_v) return;
    assert(4 == ndof_v);

    auto old_id = _id;
    auto old_ik = _ik;

    _sk.resize(ndof_v * ndof_v * _sk.size(), 0.0);
    _id.resize(ndof_v * (_id.size() - 1) + 1);
    _ik.resize(ndof_v * ndof_v * _ik.size(), -7);

    _id[0] = 0;
    for (int kold = 0; kold < Nrows(); ++kold) {
        int nr = old_id[kold + 1] - old_id[kold];
        
        for (int ii=1; ii<=ndof_v; ++ii){
			_id[ndof_v * kold + ii] = _id[ndof_v * kold + ii-1] + ndof_v * nr;
		}
  
    }

    for (int newrow = 0; newrow < ndof_v * Nrows(); ++newrow ) {
        int oldrow = newrow / ndof_v;
        int idx = _id[newrow];
        //cout << " ("<< newrow<<")  " << idx;
        for (int oid = old_id[oldrow]; oid < old_id[oldrow + 1]; ++oid) {
            int oldcol = old_ik[oid];
            
            for (int ii=0; ii < ndof_v; ++ii){
				_ik[idx + ii] = ndof_v * oldcol+ii;
				}
         
            idx += ndof_v;
        }
    }

    _nrows =  _id.size() - 1;
    _nnz   =  _ik.size();

    return;
}



void FEM_Matrix::Skalar2VectorMatrix(int ndof_v)
{
	if (1==ndof_v) return;
	
	auto old_id = _id;
	auto old_ik = _ik;
	
	_sk.resize(ndof_v*ndof_v*_sk.size());
	_id.resize(ndof_v*_id.size());
	_ik.resize(ndof_v*ndof_v*_ik.size());
    
    for(int kold = 0; kold < 14; ++kold){
    _id[4*ndof_v*kold+0] = ndof_v*old_id[kold]*ndof_v;
    _id[4*kold+1] = _id[4*kold]+ndof_v*(old_id[kold+1]-old_id[kold]);
    _id[4*kold+2] = _id[4*kold+1]+ndof_v*(old_id[kold+1]-old_id[kold]);
    _id[4*kold+3] = _id[4*kold+2]+ndof_v*(old_id[kold+2]-old_id[kold]);
     }
     
     for(int kold = 0; kold < 14; ++kold){
    _ik[16*kold+0] = old_id[kold]*ndof_v+0;
    _ik[16*kold+1] = old_id[kold]*ndof_v+1;
    _ik[16*kold+2] = old_id[kold]*ndof_v+2;
    _ik[16*kold+3] = old_id[kold]*ndof_v+3;
    
    _ik[16*kold+4] = _ik[16*kold+0];
    _ik[16*kold+5] = _ik[16*kold+1];
    _ik[16*kold+6] = _ik[16*kold+2];
    _ik[16*kold+7] = _ik[16*kold+3];
    
    _ik[16*kold+8] = _ik[16*kold+0];
    _ik[16*kold+9] = _ik[16*kold+1];
    _ik[16*kold+10] = _ik[16*kold+2];
    _ik[16*kold+11] = _ik[16*kold+3];
    
    _ik[16*kold+12] = _ik[16*kold+0];
    _ik[16*kold+13] = _ik[16*kold+1];
    _ik[16*kold+14] = _ik[16*kold+2];
    _ik[16*kold+15] = _ik[16*kold+3];
    
    
     }
    
	return;
}
*/


void FEM_Matrix::CalculateLaplace(vector<double> &f)
{
    cout << "\n############   FEM_Matrix::CalculateLaplace ";
    //double tstart = clock();
    double tstart = omp_get_wtime();                  // OpenMP
    //cout << _nnz << " vs. " << _id[_nrows] << "  " << _nrows<< endl;
    assert(_nnz == _id[_nrows]);

    for (int k = 0; k < _nrows; ++k) {
        _sk[k] = 0.0;
    }
    for (int k = 0; k < _nrows; ++k) {
        f[k] = 0.0;
    }

    //  Loop over all elements
    auto const nelem = _mesh.Nelems();
    auto const &ia   = _mesh.GetConnectivity();
    auto const &xc   = _mesh.GetCoords();

    assert(_mesh.NdofsElement() == 4);               // only for triangular, linear elements
    double ske[4][4], fe[4];
    #pragma omp parallel for private(ske,fe)
    for (int i = 0; i < nelem; ++i) {
        CalcElem(ia.data() + 4 * i, xc.data(), ske, fe);
        AddElem_3(ia.data() + 4 * i, ske, fe, f);
    }

    //double duration = (clock() - tstart) / CLOCKS_PER_SEC;
    double duration = omp_get_wtime() - tstart;             // OpenMP
    cout << "finished in  " <<  duration  << " sec.    ########\n";
    //Debug();

    return;
}

void FEM_Matrix::CalculateLaplace_heat_equation(vector<double> &f, const vector<double> &u_old, const double dt, const double t, const double c)
{
    cout << "\n############   FEM_Matrix::CalculateLaplace ";
    //double tstart = clock();
    double tstart = omp_get_wtime();                  // OpenMP
    assert(_mesh.NdofsElement() == 4);               // only for triangular, linear elements
    //cout << _nnz << " vs. " << _id[_nrows] << "  " << _nrows<< endl;
    assert(_nnz == _id[_nrows]);

    for (int k = 0; k < _nrows; ++k)
    {
        _sk[k] = 0.0;
    }
    for (int k = 0; k < _nrows; ++k)
    {
        f[k] = 0.0;
    }

    double ske[4][4], fe[4];
    //  Loop over all elements
    auto const nelem = _mesh.Nelems();
    auto const &ia   = _mesh.GetConnectivity();
    auto const &xc   = _mesh.GetCoords();

#pragma omp parallel for private(ske,fe)
    for (int i = 0; i < nelem; ++i)
    {
        CalcElem_heat_equation_crank_nichelson(ia.data()+4 * i, xc.data(), ske, fe, u_old, dt, t, c);
        //AddElem(ia.data()+3 * i, ske, fe, _id.data(), _ik.data(), _sk.data(), f.data()); // GH: deprecated
        AddElem_3(ia.data()+4 * i, ske, fe, f);
    }

    //double duration = (clock() - tstart) / CLOCKS_PER_SEC;
    double duration = omp_get_wtime() - tstart;             // OpenMP
    cout << "finished in  " <<  duration  << " sec.    ########\n";
    //Debug();

    return;
}
/*
void FEM_Matrix::ApplyDirichletBC_Box(std::vector<double> const &u, std::vector<double> &f,
            double xl, double xh, double yl, double yh,double zl, double zh )
 double const PENALTY = 1e6;
auto const idx = _mesh.Index_DirichletNodes();
int const nidx = idx.size();

for (int i=0; i<nidx; ++i)
{
int const k = idx[i];
int const id1 = fetch(k, k); // Find diagonal entry of k
assert(id1 >= 0);
_sk[id1] += PENALTY;		// matrix weighted scaling feasible
f[k] += PENALTY * u[k];
}

return;
}*/

void FEM_Matrix::ApplyDirichletBC_Box(std::vector<double> const &u, std::vector<double> &f,
            double xl, double xh, double yl, double yh,double zl, double zh )
{
    double const PENALTY = 1e6;
    // auto const idx = _mesh.Index_DirichletNodes();
//    auto const idx = _mesh.Index_DirichletNodes_Box(xl, xh, yl, yh, zl, zh);
    auto lmesh=dynamic_cast<const Mesh_3d_4_matlab&>(_mesh);  // GH: Dirty hack
    auto const idx = lmesh.Index_DirichletNodes_Box(xl, xh, yl, yh, zl, zh);

    int const nidx = idx.size();

    for (int row=0; row<nidx; ++row)
    {
        int const k = idx[row];
		int const id1 = fetch(k, k); // Find diagonal entry of row
		assert(id1 >= 0);
		_sk[id1] += PENALTY;		// matrix weighted scaling feasible
        f[k] += PENALTY * u[k];
    }

    return;
}


void FEM_Matrix::ApplyDirichletBC(std::vector<double> const &u, std::vector<double> &f)
{
    auto const idx = _mesh.Index_DirichletNodes();
    int const nidx = idx.size();

    for (int i = 0; i < nidx; ++i) {
        int const row = idx[i];
        for (int ij = _id[row]; ij < _id[row + 1]; ++ij) {
            int const col = _ik[ij];
            if (col == row) {
                _sk[ij] = 1.0;
                f[row]  = u[row];
            }
            else {
                int const id1 = fetch(col, row); // Find entry (col,row)
                assert(id1 >= 0);
                f[col] -= _sk[id1] * u[row];
                _sk[id1] = 0.0;
                _sk[ij]  = 0.0;
            }
        }
    }

    return;
}



void FEM_Matrix::AddElem_3(int const ial[4], double const ske[4][4], double const fe[4], vector<double> &f)
{
    for (int i = 0; i < 4; ++i) {
		
		for (int ki = 0; ki < 4; ++ki) {
		
        const int ii  = 4*ial[i]+ki;           // row ii (global index)
        for (int j = 0; j < 4; ++j) {     // no symmetry assumed
			
			for (int kj = 0; kj < 4; ++kj) {
				
            const int jj = 4*ial[j]+kj;        // column jj (global index)
            const int ip = fetch(ii, jj);       // find column entry jj in row ii
#ifndef NDEBUG                 // compiler option -DNDEBUG switches off the check
            if (ip < 0) {      // no entry found !!
                cout << "Error in AddElem: (" << ii << "," << jj << ") ["
                     << ial[0] << "," << ial[1] << "," << ial[2] << "," << ial[3] << "," << ial[4] <<"]\n";
                assert(ip >= 0);
            }
#endif
            #pragma omp atomic
            _sk[ip] += ske[4*i+ki][4*j+kj];
        }
        #pragma omp atomic
        f[ii] += fe[4*i+ki];
	}
    }
    }
}

bool FEM_Matrix::CheckSymmetry() const
{
    cout << "+++  Check matrix symmetry  +++" << endl;
    bool bs{true};
    for (int row = 0; row < Nrows(); ++row) {
        for (int ij = _id[row]; ij < _id[row + 1]; ++ij) {
            const int col = _ik[ij];       // column col (global index)
            const int ip = fetch(col, row);  // find column entry row in row col
            if (ip < 0) {      // no entry found !!
                cout << "Matrix has non-symmetric pattern at (" << row << "," << col << ")" << endl;
                bs = false;
                //assert(ip >= 0);
            }
            if ( std::abs(_sk[ij] - _sk[ip]) > 1e-13) {
                cout << "Matrix has non-symmetric entries at (" << row << "," << col << ")" << endl;
                bs = false;
            }
        }
    }
    return bs;
}


bool FEM_Matrix::CheckRowSum() const
{
    cout << "+++  Check row sum  +++" << endl;
    vector<double> rhs(_ncols, 1.0);  //replace Ncols() by _ncols
    vector<double> res(_nrows); //replace Nrows() by _nrows

    Mult(res, rhs);

    bool bb{true};
    for (size_t k = 0; k < res.size(); ++k) {
        //if (std::abs(res[k]) != 0.0)
        if (std::abs(res[k]) > 1e-14) {
            cout << "!! Nonzero row " << k << " : sum = " << res[k] << endl;
            bb = false;
        }
    }
    return bb;
}

bool FEM_Matrix::CheckMproperty() const
{
    cout << "+++  Check M property  +++" << endl;
    bool bm{true};
    for (int row = 0; row < Nrows(); ++row) {
        for (int ij = _id[row]; ij < _id[row + 1]; ++ij) {
            bool b_diag{true}, b_off{true};
            if (_ik[ij] == row) {
                b_diag = _sk[ij] > 0.0;
                if (!b_diag) {
                    cout << "## negative diag in row " << row << " : " << _sk[ij] << endl;
                    bm = false;
                }
            }
            else {
                b_off = _sk[ij] <= 0.0;
                if (!b_off) {
                    cout << "!! positive off-diag [" << row << "," << _ik[ij] << "] : " << _sk[ij] << endl;
                    bm = false;
                }
            }
        }
    }
    return bm;
}

bool FEM_Matrix::ForceMproperty()
{
    cout << "+++  Force M property  +++" << endl;
    bool bm{false};
    for (int row = 0; row < Nrows(); ++row) {
        double corr{0.0};
        int idiag = {-1};
        for (int ij = _id[row]; ij < _id[row + 1]; ++ij) {
            if (_ik[ij] != row &&  _sk[ij] > 0.0) {
                corr   += _sk[ij];
                _sk[ij] = 0.0;
                bm = true;
            }
            if (_ik[ij] == row) {
                idiag = ij;
            }
        }
        assert(idiag >= 0);
        _sk[idiag] += corr;
    }
    return bm;
}

bool FEM_Matrix::CheckMatrix() const
{
    bool b1 = CheckRowSum();
    if (!b1) {
        cout << " !!!!  R O W   S U M   E R R O R" << endl;
    }

    bool b2 = CheckMproperty();
    if (!b2) {
        cout << " !!!!  N O   M - M A T R I X" << endl;
    }

    return b1 && b2;
}

void FEM_Matrix::GetDiag_M(vector<double> &d) const
{
    // be carefull when using a rectangular matrix
    int const nm = min(_nrows, _ncols);

    assert( nm == static_cast<int>(d.size()) ); // instead of stopping we could resize d and warn the user

    for (int row = 0; row < Nrows(); ++row) {
        d[row] = 0.0;
        double v_ii{-1.0};
        for (int ij = _id[row]; ij < _id[row + 1]; ++ij) {
            if (_ik[ij] != row) {
                d[row] += std::abs(_sk[ij]);
            }
            else {
                v_ii = _sk[ij];
            }
        }
        if ( d[row] < v_ii ) {
            d[row] = v_ii;
        }
    }
    cout << "<<<<<<<  GetDiag_M (finished)   >>>>>>>>>" << endl;
    return;
}




// #####################################################################

BisectInterpolation::BisectInterpolation()
    : Matrix( 0, 0 ), _iv(), _vv()
{
}

BisectInterpolation::BisectInterpolation(std::vector<int> const &fathers)
    : Matrix( static_cast<int>(fathers.size()) / 2, 1 + * max_element(fathers.cbegin(), fathers.cend()) ),
      _iv(fathers), _vv(fathers.size(), 0.5)
{
}

BisectInterpolation::~BisectInterpolation()
{}

void BisectInterpolation::GetDiag(vector<double> &d) const
{
    assert( Nrows() == static_cast<int>(d.size()) );

    for (int k = 0; k < Nrows(); ++k) {
        if ( _iv[2 * k] == _iv[2 * k + 1] ) {
            d[k] = 1.0;
        }
        else {
            d[k] = 0.0;
        }
    }
    return;
}

void BisectInterpolation::Mult(vector<double> &wf, vector<double> const &uc) const
{
    assert( Nrows() == static_cast<int>(wf.size()) );
    assert( Ncols() == static_cast<int>(uc.size()) );

    #pragma omp parallel for
    for (int k = 0; k < Nrows(); ++k) {
        wf[k] = _vv[2 * k] * uc[_iv[2 * k]] + _vv[2 * k + 1] * uc[_iv[2 * k + 1]];
    }
    return;
}

// old version
//void BisectInterpolation::MultT(vector<double> const &wf, vector<double> &uc) const
//{
    //assert( Nrows() == static_cast<int>(wf.size()) );
    //assert( Ncols() == static_cast<int>(uc.size()) );
//// GH: atomic slows down the code ==> use different storage for MultT operation (CRS-matrix?)
//////#pragma omp parallel for
    //for (int k = 0; k < Ncols(); ++k)  uc[k] = 0.0;
    ////#pragma omp parallel for
    //for (int k = 0; k < Nrows(); ++k) {
        //if (_iv[2 * k] != _iv[2 * k + 1]) {
            ////#pragma omp atomic
            //uc[_iv[2 * k]  ] += _vv[2 * k  ] * wf[k];
            ////#pragma omp atomic
            //uc[_iv[2 * k + 1]] += _vv[2 * k + 1] * wf[k];
        //}
        //else {
            ////#pragma omp atomic
            //uc[_iv[2 * k]  ] +=  2.0*_vv[2 * k  ] * wf[k]; // uses a property of class BisectInterpolation
            ////uc[_iv[2 * k]  ] +=  _vv[2 * k  ] * wf[k]; // uses a property of class BisectInterpolation
        //}
    //}
    //return;
//}

void BisectInterpolation::MultT(vector<double> const &wf, vector<double> &uc) const
{
    assert( Nrows() == static_cast<int>(wf.size()) );
    assert( Ncols() == static_cast<int>(uc.size()) );
    #pragma omp parallel for
    for (int k = 0; k < Ncols(); ++k)  uc[k] = 0.0;
    
// GH: atomic slows down the code ==> use different storage for MultT operation (CRS-matrix?)
    #pragma omp parallel for
    for (int k = 0; k < Nrows(); ++k) {
        #pragma omp atomic
        uc[_iv[2 * k]  ] += _vv[2 * k  ] * wf[k];
        #pragma omp atomic
        uc[_iv[2 * k + 1]] += _vv[2 * k + 1] * wf[k];
    }
    return;
}


void BisectInterpolation::MultT_Full(vector<double> const &wf, vector<double> &uc) const
{
    assert( Nrows() == static_cast<int>(wf.size()) );
    assert( Ncols() == static_cast<int>(uc.size()) );
// GH: atomic slows down the code ==> use different storage for MultT operation (CRS-matrix?)
////#pragma omp parallel for
    for (int k = 0; k < Ncols(); ++k)  uc[k] = 0.0;
    vector<double> full(uc.size(),0.0);
    //#pragma omp parallel for
    for (int k = 0; k < Nrows(); ++k) {
        if (_iv[2 * k] != _iv[2 * k + 1]) {
            //#pragma omp atomic
            uc[_iv[2 * k]  ] += _vv[2 * k  ] * wf[k];
            //#pragma omp atomic
            uc[_iv[2 * k + 1]] += _vv[2 * k + 1] * wf[k];
            full[_iv[2 * k    ]] += _vv[2 * k  ];
            full[_iv[2 * k + 1]] += _vv[2 * k + 1];
        }
        else {
            //#pragma omp atomic
            uc[_iv[2 * k]  ] +=  2.0*_vv[2 * k  ] * wf[k]; // uses a property of class BisectInterpolation
            //uc[_iv[2 * k]  ] +=  _vv[2 * k  ] * wf[k]; // uses a property of class BisectInterpolation
            full[_iv[2 * k] ] += 2.0*_vv[2 * k  ];
        }
    }
    for (size_t k=0; k<uc.size(); ++k)  uc[k] /= full[k];
    return;
}

void BisectInterpolation::Defect(vector<double> &w,
                                 vector<double> const &f, vector<double> const &u) const
{
    assert( Nrows() == static_cast<int>(w.size()) );
    assert( Ncols() == static_cast<int>(u.size()) );
    assert( w.size() == f.size() );

    for (int k = 0; k < Nrows(); ++k) {
        w[k] = f[k] - _vv[2 * k] * u[_iv[2 * k]] + _vv[2 * k + 1] * u[_iv[2 * k + 1]];
    }
    return;
}

void BisectInterpolation::Debug() const
{
    for (int k = 0; k < Nrows(); ++k) {
        cout << k << " : fathers(" << _iv[2 * k] << "," << _iv[2 * k + 1] << ")    ";
        cout << "weights(" << _vv[2 * k] << "," << _vv[2 * k + 1] << endl;
    }
    cout << endl;
    return;
}

int BisectInterpolation::fetch(int row, int col) const
{
    int idx(-1);
    if (_iv[2 * row  ] == col) idx = 2 * row;
    if (_iv[2 * row + 1] == col) idx = 2 * row + 1;
    assert(idx >= 0);
    return idx;
}

// #####################################################################

//BisectIntDirichlet::BisectIntDirichlet(std::vector<int> const &fathers, std::vector<int> const &idxc_dir)
    //: BisectInterpolation(fathers)
//{
    //vector<bool> bdir(Ncols(), false);               // Indicator for Dirichlet coarse nodes
    //for (size_t kc = 0; kc < idxc_dir.size(); ++kc) {
        //bdir.at(idxc_dir[kc]) = true;                // Mark Dirichlet node from coarse mesh
    //}

    //for (size_t j = 0; j < _iv.size(); ++j) {
        //if ( bdir.at(_iv[j]) )   _vv[j] = 0.0;       // set weight to zero iff (at least) one father is Dirichlet node
    //}
    //return;
//}

BisectIntDirichlet::BisectIntDirichlet(std::vector<int> const &fathers, std::vector<int> const &idxc_dir)
    : BisectInterpolation(fathers), _idxDir(idxc_dir)
{
    //vector<bool> bdir(Ncols(), false);               // Indicator for Dirichlet coarse nodes
    //for (size_t kc = 0; kc < idxc_dir.size(); ++kc) {
        //bdir.at(idxc_dir[kc]) = true;                // Mark Dirichlet node from coarse mesh
    //}

    //for (size_t j = 0; j < _iv.size(); ++j) {
        //if ( bdir.at(_iv[j]) )   _vv[j] = 0.0;       // set weight to zero iff (at least) one father is Dirichlet node
    //}
    return;
}

BisectIntDirichlet::~BisectIntDirichlet()
{}


void BisectIntDirichlet::MultT(vector<double> const &wf, vector<double> &uc) const
{
    BisectInterpolation::MultT(wf, uc);
    for (size_t kc = 0; kc < _idxDir.size(); ++kc) {
        uc.at(_idxDir[kc]) = 0.0;                // Set Dirichlet node on coarse mesh to Zero
    }

    return;
}

// #####################################################################

void DefectRestrict(CRS_Matrix1 const &SK1, BisectInterpolation const &P,
                    vector<double> &fc, vector<double> &ff, vector<double> &uf)
{
    assert( P.Nrows() == static_cast<int>(ff.size()) );
    assert( P.Ncols() == static_cast<int>(fc.size()) );
    assert( ff.size() == uf.size() );
    assert( P.Nrows() == SK1.Nrows() );

//#pragma omp parallel for
    for (int k = 0; k < P.Ncols(); ++k)  fc[k] = 0.0;

// GH: atomic slows down the code ==> use different storage for MultT operation (CRS-matrix?)
    #pragma omp parallel for
    for (int row = 0; row < SK1._nrows; ++row) {
        double wi = ff[row];
        for (int ij = SK1._id[row]; ij < SK1._id[row + 1]; ++ij) {
            wi -= SK1._sk[ij] * uf[ SK1._ik[ij] ];
        }

        const int i1 = P._iv[2 * row];
        const int i2 = P._iv[2 * row + 1];
        if (i1 != i2) {
            #pragma omp atomic
            fc[i1] += P._vv[2 * row  ] * wi;
            #pragma omp atomic
            fc[i2] += P._vv[2 * row + 1] * wi;
        }
        else {
            #pragma omp atomic
            fc[i1] += 2.0 * P._vv[2 * row  ] * wi; // uses a property of class BisectInterpolation
        }
    }
    return;
}

