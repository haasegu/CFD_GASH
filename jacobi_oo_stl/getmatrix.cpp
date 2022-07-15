#include "getmatrix.h"
#include "userset.h"
#include "binaryIO.h"
#include "elements.h"              // GH
#include "sa_elements.h"
#include "utils.h"
#include "omp.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>
using namespace std;

/*
P1_2vec_Const::P1_2vec_Const(int const ial[4], std::vector<double> const &xc, 
      std::vector<double> const &u_old_n, std::vector<double> const &u_old_m)
// no member intitialization here
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
     a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
                 
    r0_n = r_old_n.at(ial[0]), r0_m = r_old_m.at(ial[0]),
                 r1_n = r_old_n.at(ial[1]), r1_m = r_old_m.at(ial[1]),
                 r2_n = r_old_n.at(ial[2]), r2_m = r_old_m.at(ial[2]),
                 r3_n = r_old_n.at(ial[3]), r3_m = r_old_m.at(ial[3]);
     
    u0_n = u_old_n.at(ial[0]), u0_m = u_old_m.at(ial[0]), v0_n = v_old_n.at(ial[0]), v0_m = v_old_m.at(ial[0]), w0_n = w_old_n.at(ial[0]), w0_m = w_old_m.at(ial[0]),
                 u1_n = u_old_n.at(ial[1]), u1_m = u_old_m.at(ial[1]), v1_n = v_old_n.at(ial[1]), v1_m = v_old_m.at(ial[1]), w1_n = w_old_n.at(ial[1]), w1_m = w_old_m.at(ial[1]),
                 u2_n = u_old_n.at(ial[2]), u2_m = u_old_m.at(ial[2]), v2_n = v_old_n.at(ial[2]), v2_m = v_old_m.at(ial[2]), w2_n = w_old_n.at(ial[2]), w2_m = w_old_m.at(ial[2]),
                 u3_n = u_old_n.at(ial[3]), u3_m = u_old_m.at(ial[3]), v3_n = v_old_n.at(ial[3]), v3_m = v_old_m.at(ial[3]), w3_n = w_old_n.at(ial[3]), w3_m = w_old_m.at(ial[3]),
                 u4_n = u_old_n.at(ial[4]), u4_m = u_old_m.at(ial[4]), v4_n = v_old_n.at(ial[4]), v4_m = v_old_m.at(ial[4]), w4_n = w_old_n.at(ial[4]), w4_m = w_old_m.at(ial[4]),
                 u5_n = u_old_n.at(ial[5]), u5_m = u_old_m.at(ial[5]), v5_n = v_old_n.at(ial[5]), v5_m = v_old_m.at(ial[5]), w5_n = w_old_n.at(ial[5]), w5_m = w_old_m.at(ial[5]),
                 u6_n = u_old_n.at(ial[6]), u6_m = u_old_m.at(ial[6]), v6_n = v_old_n.at(ial[6]), v6_m = v_old_m.at(ial[6]), w6_n = w_old_n.at(ial[6]), w6_m = w_old_m.at(ial[6]),
                 u7_n = u_old_n.at(ial[7]), u7_m = u_old_m.at(ial[7]), v7_n = v_old_n.at(ial[7]), v7_m = v_old_m.at(ial[7]), w7_n = w_old_n.at(ial[7]), w7_m = w_old_m.at(ial[7]),
                 u8_n = u_old_n.at(ial[8]), u8_m = u_old_m.at(ial[8]), v8_n = v_old_n.at(ial[8]), v8_m = v_old_m.at(ial[8]), w8_n = w_old_n.at(ial[8]), w8_m = w_old_m.at(ial[8]),
                 u9_n = u_old_n.at(ial[9]), u9_m = u_old_m.at(ial[9]), v9_n = v_old_n.at(ial[9]), v9_m = v_old_m.at(ial[9]), w9_n = w_old_n.at(ial[9]), w9_m = w_old_m.at(ial[9]);
		
//	double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;	
	
	
    f0x = -a1-a2-a3, f1x = a1, f2x = a2, f3x = a3,
    f0y = -b1-b2-b3, f1y = b1, f2y = b2, f3y = b3,
    f0z = -c1-c2-c3, f1z = c1, f2z = c2, f3z = c3;	
}*/


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
typedef double (*A00) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, 
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m);
typedef double (*B0) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double r0_n, double r1_n, double r2_n, double r3_n, double r0_m, double r1_m, double r2_m, double r3_m, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, 
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m);
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
                 
    const double r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double //u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
     
	double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
	
	A00 function00[4][10] = {{A00_00,A00_01,A00_02,A00_03},
		                     {A00_10,A00_11,A00_12,A00_13},
		                     {A00_20,A00_21,A00_22,A00_23},
		                     {A00_30,A00_31,A00_32,A00_33}};
     
       B0 function0[4] = {B0_0,B0_1,B0_2,B0_3};
             
      for(int i=0; i<=3; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function00[i][j](a,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function00[i][j](b,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function00[i][j](a,b,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function00[i][j](a,a,b,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)) 
			           +0.01224884*(function00[i][j](c,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function00[i][j](d,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function00[i][j](c,d,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function00[i][j](c,c,d,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m))
			           +0.007091003*(function00[i][j](f,e,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           + function00[i][j](e,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function00[i][j](e,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function00[i][j](e,f,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			           + function00[i][j](f,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function00[i][j](f,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m, v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)));
		  }
	   }
	
	for(int i=0; i<=9; ++i)
      {
			  fe[i]=jac*(0.01878132*(function0[i](a,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function0[i](b,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function0[i](a,b,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function0[i](a,a,b,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)) 
			  +0.01224884*(function0[i](c,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function0[i](d,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function0[i](c,d,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function0[i](c,c,d,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m))
			  +0.007091003*(function0[i](f,e,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  + function0[i](e,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function0[i](e,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function0[i](e,f,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			  + function0[i](f,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function0[i](f,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,r0_n,r1_n,r2_n,r3_n,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)));
	   }
	
}

//A01
typedef double (*A01) (double x, double y, double z, double a1, double a2, double a3, double r0_m, double r1_m, double r2_m, double r3_m);
void CalcElem_Navier_Stokes_A01(int const ial[10], double const xc[], double ske[4][10], double fe[4], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
                 
     const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     //const double u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  //u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  //v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  //w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
        
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
     
      A01 function01[4][10] = {{A01_00,A01_01,A01_02,A01_03,A01_04,A01_05,A01_06,A01_07,A01_08,A01_09},
		                       {A01_10,A01_11,A01_12,A01_13,A01_14,A01_15,A01_16,A01_17,A01_18,A01_19},
		                       {A01_20,A01_21,A01_22,A01_23,A01_24,A01_25,A01_26,A01_27,A01_28,A01_29},
		                       {A01_30,A01_31,A01_32,A01_33,A01_34,A01_35,A01_36,A01_37,A01_38,A01_39}};
     
             
      for(int i=0; i<=3; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function01[i][j](a,a,a,a1,a2,a3,r0_m,r1_m, r2_m,r3_m) 
			           +function01[i][j](b,a,a,a1,a2,a3,r0_m,r1_m, r2_m,r3_m) 
			           +function01[i][j](a,b,a,a1,a2,a3,r0_m,r1_m, r2_m,r3_m) 
			           +function01[i][j](a,a,b,a1,a2,a3,r0_m,r1_m, r2_m,r3_m)) 
			           +0.01224884*(function01[i][j](c,c,c,a1,a2,a3,r0_m,r1_m, r2_m,r3_m) 
			           +function01[i][j](d,c,c,a1,a2,a3,r0_m,r1_m, r2_m,r3_m) 
			           +function01[i][j](c,d,c,a1,a2,a3,r0_m,r1_m, r2_m,r3_m) 
			           +function01[i][j](c,c,d,a1,a2,a3,r0_m,r1_m, r2_m,r3_m))
			            +0.007091003*(function01[i][j](f,e,e,a1,a2,a3,r0_m,r1_m, r2_m,r3_m) 
			            + function01[i][j](e,f,e,a1,a2,a3,r0_m,r1_m, r2_m,r3_m) 
			            +function01[i][j](e,e,f,a1,a2,a3,r0_m,r1_m, r2_m,r3_m) 
			            +function01[i][j](e,f,f,a1,a2,a3,r0_m,r1_m, r2_m,r3_m)
			            + function01[i][j](f,e,f,a1,a2,a3,r0_m,r1_m, r2_m,r3_m) 
			            +function01[i][j](f,f,e,a1,a2,a3,r0_m,r1_m, r2_m,r3_m)));
		  }
	   }

}
//A02
typedef double (*A02) (double x, double y, double z, double b1, double b2, double b3, double r0_m, double r1_m, double r2_m, double r3_m);
void CalcElem_Navier_Stokes_A02(int const ial[10], double const xc[], double ske[4][10], double fe[4], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
                 
      const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     //const double u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  //u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  //v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  //w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
      
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
     
     A02 function01[4][10] = { {A01_00,A01_01,A01_02,A01_03,A01_04,A01_05,A01_06,A01_07,A01_08,A01_09},
		                       {A01_10,A01_11,A01_12,A01_13,A01_14,A01_15,A01_16,A01_17,A01_18,A01_19},
		                       {A01_20,A01_21,A01_22,A01_23,A01_24,A01_25,A01_26,A01_27,A01_28,A01_29},
		                       {A01_30,A01_31,A01_32,A01_33,A01_34,A01_35,A01_36,A01_37,A01_38,A01_39}};
     
             
      for(int i=0; i<=3; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function01[i][j](a,a,a,b1,b2,b3,r0_m,r1_m, r2_m,r3_m) 
			           +function01[i][j](b,a,a,b1,b2,b3,r0_m,r1_m, r2_m,r3_m) 
			           +function01[i][j](a,b,a,b1,b2,b3,r0_m,r1_m, r2_m,r3_m) 
			           +function01[i][j](a,a,b,b1,b2,b3,r0_m,r1_m, r2_m,r3_m)) 
			           +0.01224884*(function01[i][j](c,c,c,b1,b2,b3,r0_m,r1_m, r2_m,r3_m) 
			           +function01[i][j](d,c,c,b1,b2,b3,r0_m,r1_m, r2_m,r3_m) 
			           +function01[i][j](c,d,c,b1,b2,b3,r0_m,r1_m, r2_m,r3_m) 
			           +function01[i][j](c,c,d,b1,b2,b3,r0_m,r1_m, r2_m,r3_m))
			            +0.007091003*(function01[i][j](f,e,e,b1,b2,b3,r0_m,r1_m, r2_m,r3_m) 
			            + function01[i][j](e,f,e,b1,b2,b3,r0_m,r1_m, r2_m,r3_m) 
			            +function01[i][j](e,e,f,b1,b2,b3,r0_m,r1_m, r2_m,r3_m) 
			            +function01[i][j](e,f,f,b1,b2,b3,r0_m,r1_m, r2_m,r3_m)
			            + function01[i][j](f,e,f,b1,b2,b3,r0_m,r1_m, r2_m,r3_m) 
			            +function01[i][j](f,f,e,b1,b2,b3,r0_m,r1_m, r2_m,r3_m)));
		  }
	   }

    
}

//A03
typedef double (*A03) (double x, double y, double z, double c1, double c2, double c3, double r0_m, double r1_m, double r2_m, double r3_m);
void CalcElem_Navier_Stokes_A03(int const ial[10], double const xc[], double ske[4][10], double fe[4], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double //a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 //b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
     
     const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     //const double u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  //u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                 // v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                 // w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                 // w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
      
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
    
    A03 function03[4][10] = {{A03_00,A03_01,A03_02,A03_03,A03_04,A03_05,A03_06,A03_07,A03_08,A03_09},
		                     {A03_10,A03_11,A03_12,A03_13,A03_14,A03_15,A03_16,A03_17,A03_18,A03_19},
		                     {A03_20,A03_21,A03_22,A03_23,A03_24,A03_25,A03_26,A03_27,A03_28,A03_29},
		                     {A03_30,A03_31,A03_32,A03_33,A03_34,A03_35,A03_36,A03_37,A03_38,A03_39}};
    
             
      for(int i=0; i<=3; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function03[i][j](a,a,a,c1,c2,c3,r0_m,r1_m, r2_m,r3_m) 
			  +function03[i][j](b,a,a,c1,c2,c3,r0_m,r1_m, r2_m,r3_m) 
			  +function03[i][j](a,b,a,c1,c2,c3,r0_m,r1_m, r2_m,r3_m) 
			  +function03[i][j](a,a,b,c1,c2,c3,r0_m,r1_m, r2_m,r3_m)) 
			  +0.01224884*(function03[i][j](c,c,c,c1,c2,c3,r0_m,r1_m, r2_m,r3_m) 
			  +function03[i][j](d,c,c,c1,c2,c3,r0_m,r1_m, r2_m,r3_m) 
			  +function03[i][j](c,d,c,c1,c2,c3,r0_m,r1_m, r2_m,r3_m) 
			  +function03[i][j](c,c,d,c1,c2,c3,r0_m,r1_m, r2_m,r3_m))
			  +0.007091003*(function03[i][j](f,e,e,c1,c2,c3,r0_m,r1_m, r2_m,r3_m) 
			  + function03[i][j](e,f,e,c1,c2,c3,r0_m,r1_m, r2_m,r3_m) 
			  +function03[i][j](e,e,f,c1,c2,c3,r0_m,r1_m, r2_m,r3_m) 
			  +function03[i][j](e,f,f,c1,c2,c3,r0_m,r1_m, r2_m,r3_m)
			  + function03[i][j](f,e,f,c1,c2,c3,r0_m,r1_m, r2_m,r3_m) 
			  +function03[i][j](f,f,e,c1,c2,c3,r0_m,r1_m, r2_m,r3_m)));
		  }
	   }

    
}
//A10
typedef double (*A10) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, double u6_n, double u7_n, double u8_n, double u9_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, 
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m);
typedef double (*B1) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double mu, double lambda, double kp, double r0_m, double r1_m, double r2_m, double r3_m, 
double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, double u6_n, double u7_n, double u8_n, double u9_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, 
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m);
void CalcElem_Navier_Stokes_A10(int const ial[10], double const xc[], double ske[10][4], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
      
      const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
        
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
     
     A10 function10[10][4] = {{A10_00,A10_01,A10_02,A10_03},
		                      {A10_10,A10_11,A10_12,A10_13},
		                      {A10_20,A10_21,A10_22,A10_23},
		                      {A10_30,A10_31,A10_32,A10_33},
		                      {A10_40,A10_41,A10_42,A10_43},
		                      {A10_50,A10_51,A10_52,A10_53},
		                      {A10_60,A10_61,A10_62,A10_63},
		                      {A10_70,A10_71,A10_72,A10_73},
		                      {A10_80,A10_81,A10_82,A10_83},
		                      {A10_90,A10_91,A10_92,A10_93}};
		                      
	 B1 function1[10] = {B1_0,B1_1,B1_2,B1_3,B1_4,B1_5,B1_6,B1_7,B1_8,B1_9};
	 
	 for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=3; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function10[i][j](a,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			          +function10[i][j](b,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			          +function10[i][j](a,b,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			          +function10[i][j](a,a,b,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)) 
			          +0.01224884*(function10[i][j](c,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			          +function10[i][j](d,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			          +function10[i][j](c,d,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			          +function10[i][j](c,c,d,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m))
			          +0.007091003*(function10[i][j](f,e,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			          + function10[i][j](e,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			          +function10[i][j](e,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			          +function10[i][j](e,f,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			          + function10[i][j](f,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			          +function10[i][j](f,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)));
		  }
	   }
    
             
      for(int i=0; i<=9; ++i)
      {
			  fe[i]=jac*(0.01878132*(function1[i](a,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function1[i](b,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function1[i](a,b,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function1[i](a,a,b,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)) 
			       +0.01224884*(function1[i](c,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function1[i](d,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function1[i](c,d,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function1[i](c,c,d,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m))
			       +0.007091003*(function1[i](f,e,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       + function1[i](e,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function1[i](e,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function1[i](e,f,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			       + function1[i](f,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function1[i](f,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u6_n,u7_n,u8_n,u9_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)));
	   }

    
}
//A11
typedef double (*A11) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double mu, double lambda, double kp, double r0_m, double r1_m, double r2_m, double r3_m, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, 
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m );
void CalcElem_Navier_Stokes_A11(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
     
     const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double //u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
           
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
    
    A11 function11[10][10] = {{A11_00,A11_01,A11_02,A11_03,A11_04,A11_05,A11_06,A11_07,A11_08,A11_09},
		                      {A11_10,A11_11,A11_12,A11_13,A11_14,A11_15,A11_16,A11_17,A11_18,A11_19},
		                      {A11_20,A11_21,A11_22,A11_23,A11_24,A11_25,A11_26,A11_27,A11_28,A11_29},
		                      {A11_30,A11_31,A11_32,A11_33,A11_34,A11_35,A11_36,A11_37,A11_38,A11_39},
		                      {A11_40,A11_41,A11_42,A11_43,A11_44,A11_45,A11_46,A11_47,A11_48,A11_49},
		                      {A11_50,A11_51,A11_52,A11_53,A11_54,A11_55,A11_56,A11_57,A11_58,A11_59},
		                      {A11_60,A11_61,A11_62,A11_63,A11_64,A11_65,A11_66,A11_67,A11_68,A11_69},
		                      {A11_70,A11_71,A11_72,A11_73,A11_74,A11_75,A11_76,A11_77,A11_78,A11_79},
		                      {A11_80,A11_81,A11_82,A11_83,A11_84,A11_85,A11_86,A11_87,A11_88,A11_89},
		                      {A11_90,A11_91,A11_92,A11_93,A11_94,A11_95,A11_96,A11_97,A11_98,A11_99}};
  
             
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function11[i][j](a,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function11[i][j](b,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function11[i][j](a,b,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function11[i][j](a,a,b,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)) 
			           +0.01224884*(function11[i][j](c,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function11[i][j](d,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function11[i][j](c,d,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function11[i][j](c,c,d,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m))
			           +0.007091003*(function11[i][j](f,e,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           + function11[i][j](e,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function11[i][j](e,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function11[i][j](e,f,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			           + function11[i][j](f,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function11[i][j](f,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)));
		  }
	   }

  
}
//A12
typedef double (*A12) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double mu, double lambda, double r0_m, double r1_m, double r2_m, double r3_m, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m);
void CalcElem_Navier_Stokes_A12(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
     
    const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double //u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]);
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  //v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  //w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
     
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
     
     A12 function12[10][10] = {{A12_00,A12_01,A12_02,A12_03,A12_04,A12_05,A12_06,A12_07,A12_08,A12_09},
		                       {A12_10,A12_11,A12_12,A12_13,A12_14,A12_15,A12_16,A12_17,A12_18,A12_19},
		                       {A12_20,A12_21,A12_22,A12_23,A12_24,A12_25,A12_26,A12_27,A12_28,A12_29},
		                       {A12_30,A12_31,A12_32,A12_33,A12_34,A12_35,A12_36,A12_37,A12_38,A12_39},
		                       {A12_40,A12_41,A12_42,A12_43,A12_44,A12_45,A12_46,A12_47,A12_48,A12_49},
		                       {A12_50,A12_51,A12_52,A12_53,A12_54,A12_55,A12_56,A12_57,A12_58,A12_59},
		                       {A12_60,A12_61,A12_62,A12_63,A12_64,A12_65,A12_66,A12_67,A12_68,A12_69},
		                       {A12_70,A12_71,A12_72,A12_73,A12_74,A12_75,A12_76,A12_77,A12_78,A12_79},
		                       {A12_80,A12_81,A12_82,A12_83,A12_84,A12_85,A12_86,A12_87,A12_88,A12_89},
		                       {A12_90,A12_91,A12_92,A12_93,A12_94,A12_95,A12_96,A12_97,A12_98,A12_99}};
             
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function12[i][j](a,a,a,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function12[i][j](b,a,a,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function12[i][j](a,b,a,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function12[i][j](a,a,b,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m)) 
			           +0.01224884*(function12[i][j](c,c,c,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function12[i][j](d,c,c,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function12[i][j](c,d,c,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function12[i][j](c,c,d,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m))
			           +0.007091003*(function12[i][j](f,e,e,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           + function12[i][j](e,f,e,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function12[i][j](e,e,f,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function12[i][j](e,f,f,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m)
			           + function12[i][j](f,e,f,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function12[i][j](f,f,e,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m)));
		  }
	   }

    
}
//A13
typedef double (*A13) (double x, double y, double z, double a1, double a2, double a3, double c1, double c2, double c3, double mu, double lambda, double r0_m, double r1_m, double r2_m, double r3_m, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m);
void CalcElem_Navier_Stokes_A13(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 //b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
      
      const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double //u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]);
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  //v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  //w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
                     
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
     
     A13 function13[10][10] = {{A13_00,A13_01,A13_02,A13_03,A13_04,A13_05,A13_06,A13_07,A13_08,A13_09},
		                       {A13_10,A13_11,A13_12,A13_13,A13_14,A13_15,A13_16,A13_17,A13_18,A13_19},
		                       {A13_20,A13_21,A13_22,A13_23,A13_24,A13_25,A13_26,A13_27,A13_28,A13_29},
		                       {A13_30,A13_31,A13_32,A13_33,A13_34,A13_35,A13_36,A13_37,A13_38,A13_39},
		                       {A13_40,A13_41,A13_42,A13_43,A13_44,A13_45,A13_46,A13_47,A13_48,A13_49},
		                       {A13_50,A13_51,A13_52,A13_53,A13_54,A13_55,A13_56,A13_57,A13_58,A13_59},
		                       {A13_60,A13_61,A13_62,A13_63,A13_64,A13_65,A13_66,A13_67,A13_68,A13_69},
		                       {A13_70,A13_71,A13_72,A13_73,A13_74,A13_75,A13_76,A13_77,A13_78,A13_79},
		                       {A13_80,A13_81,A13_82,A13_83,A13_84,A13_85,A13_86,A13_87,A13_88,A13_89},
		                       {A13_90,A13_91,A13_92,A13_93,A13_94,A13_95,A13_96,A13_97,A13_98,A13_99}};
            
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function13[i][j](a,a,a,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function13[i][j](b,a,a,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function13[i][j](a,b,a,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function13[i][j](a,a,b,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m)) 
			           +0.01224884*(function13[i][j](c,c,c,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function13[i][j](d,c,c,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function13[i][j](c,d,c,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function13[i][j](c,c,d,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m))
			           +0.007091003*(function13[i][j](f,e,e,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           + function13[i][j](e,f,e,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function13[i][j](e,e,f,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function13[i][j](e,f,f,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m)
			           + function13[i][j](f,e,f,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m) 
			           +function13[i][j](f,f,e,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m)));
		  }
	   }

    
}
//A20
typedef double (*A20) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, 
double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, double v6_n, double v7_n, double v8_n, double v9_n, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, 
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m );
typedef double (*B2) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double mu, double lambda, double kp, double r0_m, double r1_m, double r2_m, double r3_m, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, 
double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, double v6_n, double v7_n, double v8_n, double v9_n, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, 
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m);
void CalcElem_Navier_Stokes_A20(int const ial[10], double const xc[], double ske[10][4], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
     
     const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double //u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
           
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
     
     A20 function20[10][4] = {{A20_00,A20_01,A20_02,A20_03},
		                      {A20_10,A20_11,A20_12,A20_13},
		                      {A20_20,A20_21,A20_22,A20_23},
		                      {A20_30,A20_31,A20_32,A20_33},
		                      {A20_40,A20_41,A20_42,A20_43},
		                      {A20_50,A20_51,A20_52,A20_53},
		                      {A20_60,A20_61,A20_62,A20_63},
		                      {A20_70,A20_71,A20_72,A20_73},
		                      {A20_80,A20_81,A20_82,A20_83},
		                      {A20_90,A20_91,A20_92,A20_93}};
		                      
	 B2 function2[10] = {B2_0,B2_1,B2_2,B2_3,B2_4,B2_5,B2_6,B2_7,B2_8,B2_9};
	 
	 for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=3; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function20[i][j](a,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function20[i][j](b,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function20[i][j](a,b,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function20[i][j](a,a,b,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)) 
			  +0.01224884*(function20[i][j](c,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function20[i][j](d,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function20[i][j](c,d,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function20[i][j](c,c,d,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m))
			  +0.007091003*(function20[i][j](f,e,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  + function20[i][j](e,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function20[i][j](e,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function20[i][j](e,f,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			  + function20[i][j](f,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			  +function20[i][j](f,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)));
		  }
	   }
    
             
      for(int i=0; i<=9; ++i)
      {
			  fe[i]=jac*(0.01878132*(function2[i](a,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function2[i](b,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function2[i](a,b,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function2[i](a,a,b,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)) 
			       +0.01224884*(function2[i](c,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function2[i](d,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function2[i](c,d,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function2[i](c,c,d,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m))
			       +0.007091003*(function2[i](f,e,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function2[i](e,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function2[i](e,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function2[i](e,f,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			       +function2[i](f,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			       +function2[i](f,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,v6_n,v7_n,v8_n,v9_n,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)));
	   }

   
}
//A21
typedef double (*A21) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double mu, double lambda, double r0_m, double r1_m, double r2_m, double r3_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m);
void CalcElem_Navier_Stokes_A21(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac;
                 //c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
     
     const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double //u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  //u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]);
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  //w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
           
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
    
    A21 function21[10][10] = {{A21_00,A21_01,A21_02,A21_03,A21_04,A21_05,A21_06,A21_07,A21_08,A21_09},
		                      {A21_10,A21_11,A21_12,A21_13,A21_14,A21_15,A21_16,A21_17,A21_18,A21_19},
		                      {A21_20,A21_21,A21_22,A21_23,A21_24,A21_25,A21_26,A21_27,A21_28,A21_29},
		                      {A21_30,A21_31,A21_32,A21_33,A21_34,A21_35,A21_36,A21_37,A21_38,A21_39},
		                      {A21_40,A21_41,A21_42,A21_43,A21_44,A21_45,A21_46,A21_47,A21_48,A21_49},
		                      {A21_50,A21_51,A21_52,A21_53,A21_54,A21_55,A21_56,A21_57,A21_58,A21_59},
		                      {A21_60,A21_61,A21_62,A21_63,A21_64,A21_65,A21_66,A21_67,A21_68,A21_69},
		                      {A21_70,A21_71,A21_72,A21_73,A21_74,A21_75,A21_76,A21_77,A21_78,A21_79},
		                      {A21_80,A21_81,A21_82,A21_83,A21_84,A21_85,A21_86,A21_87,A21_88,A21_89},
		                      {A21_90,A21_91,A21_92,A21_93,A21_94,A21_95,A21_96,A21_97,A21_98,A21_99}};
		                                
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function21[i][j](a,a,a,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function21[i][j](b,a,a,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function21[i][j](a,b,a,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function21[i][j](a,a,b,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m)) 
			  +0.01224884*(function21[i][j](c,c,c,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function21[i][j](d,c,c,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function21[i][j](c,d,c,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function21[i][j](c,c,d,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m))
			  +0.007091003*(function21[i][j](f,e,e,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  + function21[i][j](e,f,e,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function21[i][j](e,e,f,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function21[i][j](e,f,f,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m)
			  + function21[i][j](f,e,f,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function21[i][j](f,f,e,a1,a2,a3,b1,b2,b3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m)));
		  }
	   }

  
}
//A22
typedef double (*A22) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double mu, double lambda, double kp, double r0_m, double r1_m, double r2_m, double r3_m, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, 
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m );
void CalcElem_Navier_Stokes_A22(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
     
     const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double //u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
        
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
    
    A22 function22[10][10] = {{A22_00,A22_01,A22_02,A22_03,A22_04,A22_05,A22_06,A22_07,A22_08,A22_09},
		                      {A22_10,A22_11,A22_12,A22_13,A22_14,A22_15,A22_16,A22_17,A22_18,A22_19},
		                      {A22_20,A22_21,A22_22,A22_23,A22_24,A22_25,A22_26,A22_27,A22_28,A22_29},
		                      {A22_30,A22_31,A22_32,A22_33,A22_34,A22_35,A22_36,A22_37,A22_38,A22_39},
		                      {A22_40,A22_41,A22_42,A22_43,A22_44,A22_45,A22_46,A22_47,A22_48,A22_49},
		                      {A22_50,A22_51,A22_52,A22_53,A22_54,A22_55,A22_56,A22_57,A22_58,A22_59},
		                      {A22_60,A22_61,A22_62,A22_63,A22_64,A22_65,A22_66,A22_67,A22_68,A22_69},
		                      {A22_70,A22_71,A22_72,A22_73,A22_74,A22_75,A22_76,A22_77,A22_78,A22_79},
		                      {A22_80,A22_81,A22_82,A22_83,A22_84,A22_85,A22_86,A22_87,A22_88,A22_89},
		                      {A22_90,A22_91,A22_92,A22_93,A22_94,A22_95,A22_96,A22_97,A22_98,A22_99}};
		                               
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function22[i][j](a,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			           +function22[i][j](b,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function22[i][j](a,b,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function22[i][j](a,a,b,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)) 
			           +0.01224884*(function22[i][j](c,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function22[i][j](d,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function22[i][j](c,d,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function22[i][j](c,c,d,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m))
			           +0.007091003*(function22[i][j](f,e,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           + function22[i][j](e,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function22[i][j](e,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function22[i][j](e,f,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			           +function22[i][j](f,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function22[i][j](f,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)));
		  }
	   }

   
}
//A23
typedef double (*A23) (double x, double y, double z, double b1, double b2, double b3, double c1, double c2, double c3, double mu, double lambda, double r0_m, double r1_m, double r2_m, double r3_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m);
void CalcElem_Navier_Stokes_A23(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double //a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
    
    const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double //u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  //u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]);
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  //w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
          
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
    
    A23 function23[10][10] = {{A23_00,A23_01,A23_02,A23_03,A23_04,A23_05,A23_06,A23_07,A23_08,A23_09},
		                      {A23_10,A23_11,A23_12,A23_13,A23_14,A23_15,A23_16,A23_17,A23_18,A23_19},
		                      {A23_20,A23_21,A23_22,A23_23,A23_24,A23_25,A23_26,A23_27,A23_28,A23_29},
		                      {A23_30,A23_31,A23_32,A23_33,A23_34,A23_35,A23_36,A23_37,A23_38,A23_39},
		                      {A23_40,A23_41,A23_42,A23_43,A23_44,A23_45,A23_46,A23_47,A23_48,A23_49},
		                      {A23_50,A23_51,A23_52,A23_53,A23_54,A23_55,A23_56,A23_57,A23_58,A23_59},
		                      {A23_60,A23_61,A23_62,A23_63,A23_64,A23_65,A23_66,A23_67,A23_68,A23_69},
		                      {A23_70,A23_71,A23_72,A23_73,A23_74,A23_75,A23_76,A23_77,A23_78,A23_79},
		                      {A23_80,A23_81,A23_82,A23_83,A23_84,A23_85,A23_86,A23_87,A23_88,A23_89},
		                      {A23_90,A23_91,A23_92,A23_93,A23_94,A23_95,A23_96,A23_97,A23_98,A23_99}};
		                              
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function23[i][j](a,a,a,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function23[i][j](b,a,a,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function23[i][j](a,b,a,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function23[i][j](a,a,b,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m)) 
			  +0.01224884*(function23[i][j](c,c,c,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function23[i][j](d,c,c,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function23[i][j](c,d,c,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function23[i][j](c,c,d,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m))
			  +0.007091003*(function23[i][j](f,e,e,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function23[i][j](e,f,e,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function23[i][j](e,e,f,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function23[i][j](e,f,f,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m)
			  +function23[i][j](f,e,f,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			  +function23[i][j](f,f,e,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m)));
		  }
	   }

    

}
//A30
typedef double (*A30) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, 
double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, double w6_n, double w7_n, double w8_n, double w9_n, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m);
typedef double (*B3) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double mu, double lambda, double kp, double r0_m, double r1_m, double r2_m, double r3_m, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, 
double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, double w6_n, double w7_n, double w8_n, double w9_n, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m );
void CalcElem_Navier_Stokes_A30(int const ial[10], double const xc[], double ske[10][4], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
    
    const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double //u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                  w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
           
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
    
    A30 function30[10][4] = {{A30_00,A30_01,A30_02,A30_03},
		                     {A30_10,A30_11,A30_12,A30_13},
		                     {A30_20,A30_21,A30_22,A30_23},
		                     {A30_30,A30_31,A30_32,A30_33},
		                     {A30_40,A30_41,A30_42,A30_43},
		                     {A30_50,A30_51,A30_52,A30_53},
		                     {A30_60,A30_61,A30_62,A30_63},
		                     {A30_70,A30_71,A30_72,A30_73},
		                     {A30_80,A30_81,A30_82,A30_83},
		                     {A30_90,A30_91,A30_92,A30_93}};
		                      
	 B3 function3[10] = {B3_0,B3_1,B3_2,B3_3,B3_4,B3_5,B3_6,B3_7,B3_8,B3_9};
	 
	 for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=3; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function30[i][j](a,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function30[i][j](b,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function30[i][j](a,b,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function30[i][j](a,a,b,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)) 
			           +0.01224884*(function30[i][j](c,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function30[i][j](d,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function30[i][j](c,d,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function30[i][j](c,c,d,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m))
			           +0.007091003*(function30[i][j](f,e,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function30[i][j](e,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function30[i][j](e,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function30[i][j](e,f,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			           +function30[i][j](f,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function30[i][j](f,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m, w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)));
		  }
	   }
    
             
      for(int i=0; i<=9; ++i)
      {
			  fe[i]=jac*(0.01878132*(function3[i](a,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			      +function3[i](b,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			      +function3[i](a,b,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			      +function3[i](a,a,b,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)) 
			      +0.01224884*(function3[i](c,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			      +function3[i](d,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			      +function3[i](c,d,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			      +function3[i](c,c,d,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m))
			      +0.007091003*(function3[i](f,e,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			      + function3[i](e,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			      +function3[i](e,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			      +function3[i](e,f,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			      + function3[i](f,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			      +function3[i](f,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,w6_n,w7_n,w8_n,w9_n,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)));
	   }

  
}
//A31
typedef double (*A31) (double x, double y, double z, double a1, double a2, double a3, double c1, double c2, double c3, double mu, double lambda, double r0_m, double r1_m, double r2_m, double r3_m, 
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m);
void CalcElem_Navier_Stokes_A31(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];//i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 //b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
     
     const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double //u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  //u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  //v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
      
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
     
     A31 function31[10][10] = {{A31_00,A31_01,A31_02,A31_03,A31_04,A31_05,A31_06,A31_07,A31_08,A31_09},
		                       {A31_10,A31_11,A31_12,A31_13,A31_14,A31_15,A31_16,A31_17,A31_18,A31_19},
		                       {A31_20,A31_21,A31_22,A31_23,A31_24,A31_25,A31_26,A31_27,A31_28,A31_29},
		                       {A31_30,A31_31,A31_32,A31_33,A31_34,A31_35,A31_36,A31_37,A31_38,A31_39},
		                       {A31_40,A31_41,A31_42,A31_43,A31_44,A31_45,A31_46,A31_47,A31_48,A31_49},
		                       {A31_50,A31_51,A31_52,A31_53,A31_54,A31_55,A31_56,A31_57,A31_58,A31_59},
		                       {A31_60,A31_61,A31_62,A31_63,A31_64,A31_65,A31_66,A31_67,A31_68,A31_69},
		                       {A31_70,A31_71,A31_72,A31_73,A31_74,A31_75,A31_76,A31_77,A31_78,A31_79},
		                       {A31_80,A31_81,A31_82,A31_83,A31_84,A31_85,A31_86,A31_87,A31_88,A31_89},
		                       {A31_90,A31_91,A31_92,A31_93,A31_94,A31_95,A31_96,A31_97,A31_98,A31_99}};
		                                
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function31[i][j](a,a,a,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function31[i][j](b,a,a,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			           +function31[i][j](a,b,a,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function31[i][j](a,a,b,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)) 
			           +0.01224884*(function31[i][j](c,c,c,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function31[i][j](d,c,c,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function31[i][j](c,d,c,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function31[i][j](c,c,d,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m))
			           +0.007091003*(function31[i][j](f,e,e,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           + function31[i][j](e,f,e,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function31[i][j](e,e,f,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function31[i][j](e,f,f,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			           +function31[i][j](f,e,f,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function31[i][j](f,f,e,a1,a2,a3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)));
		  }
	   }

    
}
//A32
typedef double (*A32) (double x, double y, double z, double b1, double b2, double b3, double c1, double c2, double c3, double mu, double lambda, double r0_m, double r1_m, double r2_m, double r3_m, 
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m);
void CalcElem_Navier_Stokes_A32(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double //a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
    
    const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double //u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  //u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]);
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  //w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
       
    double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
    
    A32 function32[10][10] = { {A32_00,A32_01,A32_02,A32_03,A32_04,A32_05,A32_06,A32_07,A32_08,A32_09},
		                       {A32_10,A32_11,A32_12,A32_13,A32_14,A32_15,A32_16,A32_17,A32_18,A32_19},
		                       {A32_20,A32_21,A32_22,A32_23,A32_24,A32_25,A32_26,A32_27,A32_28,A32_29},
		                       {A32_30,A32_31,A32_32,A32_33,A32_34,A32_35,A32_36,A32_37,A32_38,A32_39},
		                       {A32_40,A32_41,A32_42,A32_43,A32_44,A32_45,A32_46,A32_47,A32_48,A32_49},
		                       {A32_50,A32_51,A32_52,A32_53,A32_54,A32_55,A32_56,A32_57,A32_58,A32_59},
		                       {A32_60,A32_61,A32_62,A32_63,A32_64,A32_65,A32_66,A32_67,A32_68,A32_69},
		                       {A32_70,A32_71,A32_72,A32_73,A32_74,A32_75,A32_76,A32_77,A32_78,A32_79},
		                       {A32_80,A32_81,A32_82,A32_83,A32_84,A32_85,A32_86,A32_87,A32_88,A32_89},
		                       {A32_90,A32_91,A32_92,A32_93,A32_94,A32_95,A32_96,A32_97,A32_98,A32_99}};
		                                
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function32[i][j](a,a,a,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			           +function32[i][j](b,a,a,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			           +function32[i][j](a,b,a,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			           +function32[i][j](a,a,b,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m)) 
			           +0.01224884*(function32[i][j](c,c,c,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			           +function32[i][j](d,c,c,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			           +function32[i][j](c,d,c,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			           +function32[i][j](c,c,d,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m))
			           +0.007091003*(function32[i][j](f,e,e,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			           + function32[i][j](e,f,e,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			           +function32[i][j](e,e,f,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			           +function32[i][j](e,f,f,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m)
			           +function32[i][j](f,e,f,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m) 
			           +function32[i][j](f,f,e,b1,b2,b3,c1,c2,c3,mu,lambda,r0_m,r1_m,r2_m,r3_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m)));
		  }
	   }

   
}
//A33
typedef double (*A33) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double mu, double lambda, double kp, double r0_m, double r1_m, double r2_m, double r3_m, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, 
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m);
void CalcElem_Navier_Stokes_A33(int const ial[10], double const xc[], double ske[10][10], double fe[10], const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, 
const double dt, const double t, const double mu, const double lambda, const double kp)
{
    const int    i1  = 3 * ial[0],   i2 = 3 * ial[1],   i3 = 3 * ial[2], i4 = 3 * ial[3];// i5 = 3 * ial[5], i6 = 3 * ial[6], i7 = 3 * ial[7], i8 = 3 * ial[8], i9 = 3 * ial[9], i10 = 3 * ial[10];
    const double x1 = xc[i1 + 0] - xc[i4 + 0],  y1 = xc[i1 + 1] - xc[i4 + 1], z1 = xc[i1 + 2] - xc[i4 + 2],
                 x2 = xc[i2 + 0] - xc[i4 + 0],  y2 = xc[i2 + 1] - xc[i4 + 1], z2 = xc[i2 + 2] - xc[i4 + 2],
                 x3 = xc[i3 + 0] - xc[i4 + 0],  y3 = xc[i3 + 1] - xc[i4 + 1], z3 = xc[i3 + 2] - xc[i4 + 2];
                 
    const double jac = fabs(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
    //coefficients for derivatives
    const double a1 = (y2*z3-y3*z2)/jac, a2 = (y3*z1-y1*z3)/jac, a3 = (y1*z2-y2*z1)/jac, 
                 b1 = (x3*z2-x2*z3)/jac, b2 = (x3*y1-x1*y3)/jac, b3 = (x2*z1-x1*z2)/jac, 
                 c1 = (x2*y3-x3*y2)/jac, c2 = (x1*z3-x3*z1)/jac, c3 = (x1*y2-x2*y1)/jac;
     
     const double //r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), r3_n = r_old_n.at(ial[3]),
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]), r3_m = r_old_m.at(ial[3]);
                             
     const double //u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), u6_n = u_old_n.at(ial[6]), u7_n = u_old_n.at(ial[7]), u8_n = u_old_n.at(ial[8]), u9_n = u_old_n.at(ial[9]),
                  u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), u6_m = u_old_m.at(ial[6]), u7_m = u_old_m.at(ial[7]), u8_m = u_old_m.at(ial[8]), u9_m = u_old_m.at(ial[9]),
                  //v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), v6_n = v_old_n.at(ial[6]), v7_n = v_old_n.at(ial[7]), v8_n = v_old_n.at(ial[8]), v9_n = v_old_n.at(ial[9]),
                  v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), v6_m = v_old_m.at(ial[6]), v7_m = v_old_m.at(ial[7]), v8_m = v_old_m.at(ial[8]), v9_m = v_old_m.at(ial[9]),
                  //w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), w6_n = w_old_n.at(ial[6]), w7_n = w_old_n.at(ial[7]), w8_n = w_old_n.at(ial[8]), w9_n = w_old_n.at(ial[9]),
                  w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]), w6_m = w_old_m.at(ial[6]), w7_m = w_old_m.at(ial[7]), w8_m = w_old_m.at(ial[8]), w9_m = w_old_m.at(ial[9]);
     
        
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
     
     A33 function33[10][10] = {{A33_00,A33_01,A33_02,A33_03,A33_04,A33_05,A33_06,A33_07,A33_08,A33_09},
		                       {A33_10,A33_11,A33_12,A33_13,A33_14,A33_15,A33_16,A33_17,A33_18,A33_19},
		                       {A33_20,A33_21,A33_22,A33_23,A33_24,A33_25,A33_26,A33_27,A33_28,A33_29},
		                       {A33_30,A33_31,A33_32,A33_33,A33_34,A33_35,A33_36,A33_37,A33_38,A33_39},
		                       {A33_40,A33_41,A33_42,A33_43,A33_44,A33_45,A33_46,A33_47,A33_48,A33_49},
		                       {A33_50,A33_51,A33_52,A33_53,A33_54,A33_55,A33_56,A33_57,A33_58,A33_59},
		                       {A33_60,A33_61,A33_62,A33_63,A33_64,A33_65,A33_66,A33_67,A33_68,A33_69},
		                       {A33_70,A33_71,A33_72,A33_73,A33_74,A33_75,A33_76,A33_77,A33_78,A33_79},
		                       {A33_80,A33_81,A33_82,A33_83,A33_84,A33_85,A33_86,A33_87,A33_88,A33_89},
		                       {A33_90,A33_91,A33_92,A33_93,A33_94,A33_95,A33_96,A33_97,A33_98,A33_99}};
		                                
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function33[i][j](a,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function33[i][j](b,a,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function33[i][j](a,b,a,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function33[i][j](a,a,b,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)) 
			           +0.01224884*(function33[i][j](c,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function33[i][j](d,c,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function33[i][j](c,d,c,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function33[i][j](c,c,d,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m))
			           +0.007091003*(function33[i][j](f,e,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function33[i][j](e,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function33[i][j](e,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function33[i][j](e,f,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)
			           +function33[i][j](f,e,f,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m) 
			           +function33[i][j](f,f,e,a1,a2,a3,b1,b2,b3,c1,c2,c3,dt,mu,lambda,kp,r0_m,r1_m,r2_m,r3_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,u6_m,u7_m,u8_m,u9_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,v6_m,v7_m,v8_m,v9_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m,w6_m,w7_m,w8_m,w9_m)));
		  }
	   }

  
}
/*typedef double (*A00) (double x, double y, double z);
typedef double (*B0) (double x, double y, double z);
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
	
	A00 function00[4][10] = {{A00_00,A00_01,A00_02,A00_03},
		                     {A00_10,A00_11,A00_12,A00_13},
		                     {A00_20,A00_21,A00_22,A00_23},
		                     {A00_30,A00_31,A00_32,A00_33}};
     
       B0 function0[4] = {B0_0,B0_1,B0_2,B0_3};
             
      for(int i=0; i<=3; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function00[i][j](a,a,a) +function00[i][j](b,a,a) +function00[i][j](a,b,a) +function00[i][j](a,a,b)) +0.01224884*(function00[i][j](c,c,c) +function00[i][j](d,c,c) +function00[i][j](c,d,c) +function00[i][j](c,c,d))
			            +0.007091003*(function00[i][j](f,e,e) + function00[i][j](e,f,e) +function00[i][j](e,e,f) +function00[i][j](e,f,f)+ function00[i][j](f,e,f) +function00[i][j](f,f,e)));
		  }
	   }
	
	for(int i=0; i<=9; ++i)
      {
			  fe[i]=jac*(0.01878132*(function0[i](a,a,a) +function0[i](b,a,a) +function0[i](a,b,a) +function0[i](a,a,b)) +0.01224884*(function0[i](c,c,c) +function0[i](d,c,c) +function0[i](c,d,c) +function0[i](c,c,d))
			            +0.007091003*(function0[i](f,e,e) + function0[i](e,f,e) +function0[i](e,e,f) +function0[i](e,f,f)+ function0[i](f,e,f) +function0[i](f,f,e)));
	   }
	

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
     
      A01 function01[4][10] = {{A01_00,A01_01,A01_02,A01_03,A01_04,A01_05,A01_06,A01_07,A01_08,A01_09},
		                       {A01_10,A01_11,A01_12,A01_13,A01_14,A01_15,A01_16,A01_17,A01_18,A01_19},
		                       {A01_20,A01_21,A01_22,A01_23,A01_24,A01_25,A01_26,A01_27,A01_28,A01_29},
		                       {A01_30,A01_31,A01_32,A01_33,A01_34,A01_35,A01_36,A01_37,A01_38,A01_39}};
     
             
      for(int i=0; i<=3; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function01[i][j](a,a,a) +function01[i][j](b,a,a) +function01[i][j](a,b,a) +function01[i][j](a,a,b)) +0.01224884*(function01[i][j](c,c,c) +function01[i][j](d,c,c) +function01[i][j](c,d,c) +function01[i][j](c,c,d))
			            +0.007091003*(function01[i][j](f,e,e) + function01[i][j](e,f,e) +function01[i][j](e,e,f) +function01[i][j](e,f,f)+ function01[i][j](f,e,f) +function01[i][j](f,f,e)));
		  }
	   }

    ske[0][0] = 0.01878132*(A01_00(a,a,a) +A01_00(b,a,a) +A01_00(a,b,a) +A01_00(a,a,b)) +0.01224884*(A01_00(c,c,c) +A01_00(d,c,c) +A01_00(c,d,c) +A01_00(c,c,d))+0.007091003*(A01_00(f,e,e) + A01_00(e,f,e) +A01_00(e,e,f) +A01_00(e,f,f)+ A01_00(f,e,f) +A01_00(f,f,e));
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

    fe[0] = fe[1] = fe[2]= fe[3] =  0;
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
     
     A02 function01[4][10] = { {A01_00,A01_01,A01_02,A01_03,A01_04,A01_05,A01_06,A01_07,A01_08,A01_09},
		                       {A01_10,A01_11,A01_12,A01_13,A01_14,A01_15,A01_16,A01_17,A01_18,A01_19},
		                       {A01_20,A01_21,A01_22,A01_23,A01_24,A01_25,A01_26,A01_27,A01_28,A01_29},
		                       {A01_30,A01_31,A01_32,A01_33,A01_34,A01_35,A01_36,A01_37,A01_38,A01_39}};
     
             
      for(int i=0; i<=3; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function01[i][j](a,a,a) +function01[i][j](b,a,a) +function01[i][j](a,b,a) +function01[i][j](a,a,b)) +0.01224884*(function01[i][j](c,c,c) +function01[i][j](d,c,c) +function01[i][j](c,d,c) +function01[i][j](c,c,d))
			            +0.007091003*(function01[i][j](f,e,e) + function01[i][j](e,f,e) +function01[i][j](e,e,f) +function01[i][j](e,f,f)+ function01[i][j](f,e,f) +function01[i][j](f,f,e)));
		  }
	   }

    
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
    
    A03 function03[4][10] = {{A03_00,A03_01,A03_02,A03_03,A03_04,A03_05,A03_06,A03_07,A03_08,A03_09},
		                     {A03_10,A03_11,A03_12,A03_13,A03_14,A03_15,A03_16,A03_17,A03_18,A03_19},
		                     {A03_20,A03_21,A03_22,A03_23,A03_24,A03_25,A03_26,A03_27,A03_28,A03_29},
		                     {A03_30,A03_31,A03_32,A03_33,A03_34,A03_35,A03_36,A03_37,A03_38,A03_39}};
    
             
      for(int i=0; i<=3; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function03[i][j](a,a,a) +function03[i][j](b,a,a) +function03[i][j](a,b,a) +function03[i][j](a,a,b)) +0.01224884*(function03[i][j](c,c,c) +function03[i][j](d,c,c) +function03[i][j](c,d,c) +function03[i][j](c,c,d))
			            +0.007091003*(function03[i][j](f,e,e) + function03[i][j](e,f,e) +function03[i][j](e,e,f) +function03[i][j](e,f,f)+ function03[i][j](f,e,f) +function03[i][j](f,f,e)));
		  }
	   }

    
}
//A10
typedef double (*A10) (double x, double y, double z);
typedef double (*B1) (double x, double y, double z);
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
     
     A10 function10[10][4] = {{A10_00,A10_01,A10_02,A10_03},
		                      {A10_10,A10_11,A10_12,A10_13},
		                      {A10_20,A10_21,A10_22,A10_23},
		                      {A10_30,A10_31,A10_32,A10_33},
		                      {A10_40,A10_41,A10_42,A10_43},
		                      {A10_50,A10_51,A10_52,A10_53},
		                      {A10_60,A10_61,A10_62,A10_63},
		                      {A10_70,A10_71,A10_72,A10_73},
		                      {A10_80,A10_81,A10_82,A10_83},
		                      {A10_90,A10_91,A10_92,A10_93}};
		                      
	 B1 function1[10] = {B1_0,B1_1,B1_2,B1_3,B1_4,B1_5,B1_6,B1_7,B1_8,B1_9};
	 
	 for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=3; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function10[i][j](a,a,a) +function10[i][j](b,a,a) +function10[i][j](a,b,a) +function10[i][j](a,a,b)) +0.01224884*(function10[i][j](c,c,c) +function10[i][j](d,c,c) +function10[i][j](c,d,c) +function10[i][j](c,c,d))
			            +0.007091003*(function10[i][j](f,e,e) + function10[i][j](e,f,e) +function10[i][j](e,e,f) +function10[i][j](e,f,f)+ function10[i][j](f,e,f) +function10[i][j](f,f,e)));
		  }
	   }
    
             
      for(int i=0; i<=9; ++i)
      {
			  fe[i]=jac*(0.01878132*(function1[i](a,a,a) +function1[i](b,a,a) +function1[i](a,b,a) +function1[i](a,a,b)) +0.01224884*(function1[i](c,c,c) +function1[i](d,c,c) +function1[i](c,d,c) +function1[i](c,c,d))
			            +0.007091003*(function1[i](f,e,e) + function1[i](e,f,e) +function1[i](e,e,f) +function1[i](e,f,f)+ function1[i](f,e,f) +function1[i](f,f,e)));
	   }

    
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
    
    A11 function11[10][10] = {{A11_00,A11_01,A11_02,A11_03,A11_04,A11_05,A11_06,A11_07,A11_08,A11_09},
		                      {A11_10,A11_11,A11_12,A11_13,A11_14,A11_15,A11_16,A11_17,A11_18,A11_19},
		                      {A11_20,A11_21,A11_22,A11_23,A11_24,A11_25,A11_26,A11_27,A11_28,A11_29},
		                      {A11_30,A11_31,A11_32,A11_33,A11_34,A11_35,A11_36,A11_37,A11_38,A11_39},
		                      {A11_40,A11_41,A11_42,A11_43,A11_44,A11_45,A11_46,A11_47,A11_48,A11_49},
		                      {A11_50,A11_51,A11_52,A11_53,A11_54,A11_55,A11_56,A11_57,A11_58,A11_59},
		                      {A11_60,A11_61,A11_62,A11_63,A11_64,A11_65,A11_66,A11_67,A11_68,A11_69},
		                      {A11_70,A11_71,A11_72,A11_73,A11_74,A11_75,A11_76,A11_77,A11_78,A11_79},
		                      {A11_80,A11_81,A11_82,A11_83,A11_84,A11_85,A11_86,A11_87,A11_88,A11_89},
		                      {A11_90,A11_91,A11_92,A11_93,A11_94,A11_95,A11_96,A11_97,A11_98,A11_99}};
  
             
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function11[i][j](a,a,a) +function11[i][j](b,a,a) +function11[i][j](a,b,a) +function11[i][j](a,a,b)) +0.01224884*(function11[i][j](c,c,c) +function11[i][j](d,c,c) +function11[i][j](c,d,c) +function11[i][j](c,c,d))
			            +0.007091003*(function11[i][j](f,e,e) + function11[i][j](e,f,e) +function11[i][j](e,e,f) +function11[i][j](e,f,f)+ function11[i][j](f,e,f) +function11[i][j](f,f,e)));
		  }
	   }

  
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
     
     A12 function12[10][10] = {{A12_00,A12_01,A12_02,A12_03,A12_04,A12_05,A12_06,A12_07,A12_08,A12_09},
		                       {A12_10,A12_11,A12_12,A12_13,A12_14,A12_15,A12_16,A12_17,A12_18,A12_19},
		                       {A12_20,A12_21,A12_22,A12_23,A12_24,A12_25,A12_26,A12_27,A12_28,A12_29},
		                       {A12_30,A12_31,A12_32,A12_33,A12_34,A12_35,A12_36,A12_37,A12_38,A12_39},
		                       {A12_40,A12_41,A12_42,A12_43,A12_44,A12_45,A12_46,A12_47,A12_48,A12_49},
		                       {A12_50,A12_51,A12_52,A12_53,A12_54,A12_55,A12_56,A12_57,A12_58,A12_59},
		                       {A12_60,A12_61,A12_62,A12_63,A12_64,A12_65,A12_66,A12_67,A12_68,A12_69},
		                       {A12_70,A12_71,A12_72,A12_73,A12_74,A12_75,A12_76,A12_77,A12_78,A12_79},
		                       {A12_80,A12_81,A12_82,A12_83,A12_84,A12_85,A12_86,A12_87,A12_88,A12_89},
		                       {A12_90,A12_91,A12_92,A12_93,A12_94,A12_95,A12_96,A12_97,A12_98,A12_99}};
             
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function12[i][j](a,a,a) +function12[i][j](b,a,a) +function12[i][j](a,b,a) +function12[i][j](a,a,b)) +0.01224884*(function12[i][j](c,c,c) +function12[i][j](d,c,c) +function12[i][j](c,d,c) +function12[i][j](c,c,d))
			            +0.007091003*(function12[i][j](f,e,e) + function12[i][j](e,f,e) +function12[i][j](e,e,f) +function12[i][j](e,f,f)+ function12[i][j](f,e,f) +function12[i][j](f,f,e)));
		  }
	   }

    
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
      
      const double r0_n = r_old_n.at(ial[0]), r0_m = r_old_m.at(ial[0]),
                 r1_n = r_old_n.at(ial[1]), r1_m = r_old_m.at(ial[1]),
                 r2_n = r_old_n.at(ial[2]), r2_m = r_old_m.at(ial[2]),
                 r3_n = r_old_n.at(ial[3]), r3_m = r_old_m.at(ial[3]);
                            
     double a=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
     
     A13 function13[10][10] = {{A13_00,A13_01,A13_02,A13_03,A13_04,A13_05,A13_06,A13_07,A13_08,A13_09},
		                       {A13_10,A13_11,A13_12,A13_13,A13_14,A13_15,A13_16,A13_17,A13_18,A13_19},
		                       {A13_20,A13_21,A13_22,A13_23,A13_24,A13_25,A13_26,A13_27,A13_28,A13_29},
		                       {A13_30,A13_31,A13_32,A13_33,A13_34,A13_35,A13_36,A13_37,A13_38,A13_39},
		                       {A13_40,A13_41,A13_42,A13_43,A13_44,A13_45,A13_46,A13_47,A13_48,A13_49},
		                       {A13_50,A13_51,A13_52,A13_53,A13_54,A13_55,A13_56,A13_57,A13_58,A13_59},
		                       {A13_60,A13_61,A13_62,A13_63,A13_64,A13_65,A13_66,A13_67,A13_68,A13_69},
		                       {A13_70,A13_71,A13_72,A13_73,A13_74,A13_75,A13_76,A13_77,A13_78,A13_79},
		                       {A13_80,A13_81,A13_82,A13_83,A13_84,A13_85,A13_86,A13_87,A13_88,A13_89},
		                       {A13_90,A13_91,A13_92,A13_93,A13_94,A13_95,A13_96,A13_97,A13_98,A13_99}};
            
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function13[i][j](a,a,a) +function13[i][j](b,a,a) +function13[i][j](a,b,a) +function13[i][j](a,a,b)) +0.01224884*(function13[i][j](c,c,c) +function13[i][j](d,c,c) +function13[i][j](c,d,c) +function13[i][j](c,c,d))
			            +0.007091003*(function13[i][j](f,e,e) + function13[i][j](e,f,e) +function13[i][j](e,e,f) +function13[i][j](e,f,f)+ function13[i][j](f,e,f) +function13[i][j](f,f,e)));
		  }
	   }

    
}
//A20
typedef double (*A20) (double x, double y, double z);
typedef double (*B2) (double x, double y, double z);
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
     
     A20 function20[10][4] = {{A20_00,A20_01,A20_02,A20_03},
		                      {A20_10,A20_11,A20_12,A20_13},
		                      {A20_20,A20_21,A20_22,A20_23},
		                      {A20_30,A20_31,A20_32,A20_33},
		                      {A20_40,A20_41,A20_42,A20_43},
		                      {A20_50,A20_51,A20_52,A20_53},
		                      {A20_60,A20_61,A20_62,A20_63},
		                      {A20_70,A20_71,A20_72,A20_73},
		                      {A20_80,A20_81,A20_82,A20_83},
		                      {A20_90,A20_91,A20_92,A20_93}};
		                      
	 B2 function2[10] = {B2_0,B2_1,B2_2,B2_3,B2_4,B2_5,B2_6,B2_7,B2_8,B2_9};
	 
	 for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=3; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function20[i][j](a,a,a) +function20[i][j](b,a,a) +function20[i][j](a,b,a) +function20[i][j](a,a,b)) +0.01224884*(function20[i][j](c,c,c) +function20[i][j](d,c,c) +function20[i][j](c,d,c) +function20[i][j](c,c,d))
			            +0.007091003*(function20[i][j](f,e,e) + function20[i][j](e,f,e) +function20[i][j](e,e,f) +function20[i][j](e,f,f)+ function20[i][j](f,e,f) +function20[i][j](f,f,e)));
		  }
	   }
    
             
      for(int i=0; i<=9; ++i)
      {
			  fe[i]=jac*(0.01878132*(function2[i](a,a,a) +function2[i](b,a,a) +function2[i](a,b,a) +function2[i](a,a,b)) +0.01224884*(function2[i](c,c,c) +function2[i](d,c,c) +function2[i](c,d,c) +function2[i](c,c,d))
			            +0.007091003*(function2[i](f,e,e) + function2[i](e,f,e) +function2[i](e,e,f) +function2[i](e,f,f)+ function2[i](f,e,f) +function2[i](f,f,e)));
	   }

   
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
    
    A21 function21[10][10] = {{A21_00,A21_01,A21_02,A21_03,A21_04,A21_05,A21_06,A21_07,A21_08,A21_09},
		                      {A21_10,A21_11,A21_12,A21_13,A21_14,A21_15,A21_16,A21_17,A21_18,A21_19},
		                      {A21_20,A21_21,A21_22,A21_23,A21_24,A21_25,A21_26,A21_27,A21_28,A21_29},
		                      {A21_30,A21_31,A21_32,A21_33,A21_34,A21_35,A21_36,A21_37,A21_38,A21_39},
		                      {A21_40,A21_41,A21_42,A21_43,A21_44,A21_45,A21_46,A21_47,A21_48,A21_49},
		                      {A21_50,A21_51,A21_52,A21_53,A21_54,A21_55,A21_56,A21_57,A21_58,A21_59},
		                      {A21_60,A21_61,A21_62,A21_63,A21_64,A21_65,A21_66,A21_67,A21_68,A21_69},
		                      {A21_70,A21_71,A21_72,A21_73,A21_74,A21_75,A21_76,A21_77,A21_78,A21_79},
		                      {A21_80,A21_81,A21_82,A21_83,A21_84,A21_85,A21_86,A21_87,A21_88,A21_89},
		                      {A21_90,A21_91,A21_92,A21_93,A21_94,A21_95,A21_96,A21_97,A21_98,A21_99}};
		                                
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function21[i][j](a,a,a) +function21[i][j](b,a,a) +function21[i][j](a,b,a) +function21[i][j](a,a,b)) +0.01224884*(function21[i][j](c,c,c) +function21[i][j](d,c,c) +function21[i][j](c,d,c) +function21[i][j](c,c,d))
			            +0.007091003*(function21[i][j](f,e,e) + function21[i][j](e,f,e) +function21[i][j](e,e,f) +function21[i][j](e,f,f)+ function21[i][j](f,e,f) +function21[i][j](f,f,e)));
		  }
	   }

  
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
    
    A22 function22[10][10] = {{A22_00,A22_01,A22_02,A22_03,A22_04,A22_05,A22_06,A22_07,A22_08,A22_09},
		                      {A22_10,A22_11,A22_12,A22_13,A22_14,A22_15,A22_16,A22_17,A22_18,A22_19},
		                      {A22_20,A22_21,A22_22,A22_23,A22_24,A22_25,A22_26,A22_27,A22_28,A22_29},
		                      {A22_30,A22_31,A22_32,A22_33,A22_34,A22_35,A22_36,A22_37,A22_38,A22_39},
		                      {A22_40,A22_41,A22_42,A22_43,A22_44,A22_45,A22_46,A22_47,A22_48,A22_49},
		                      {A22_50,A22_51,A22_52,A22_53,A22_54,A22_55,A22_56,A22_57,A22_58,A22_59},
		                      {A22_60,A22_61,A22_62,A22_63,A22_64,A22_65,A22_66,A22_67,A22_68,A22_69},
		                      {A22_70,A22_71,A22_72,A22_73,A22_74,A22_75,A22_76,A22_77,A22_78,A22_79},
		                      {A22_80,A22_81,A22_82,A22_83,A22_84,A22_85,A22_86,A22_87,A22_88,A22_89},
		                      {A22_90,A22_91,A22_92,A22_93,A22_94,A22_95,A22_96,A22_97,A22_98,A22_99}};
		                               
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function22[i][j](a,a,a) +function22[i][j](b,a,a) +function22[i][j](a,b,a) +function22[i][j](a,a,b)) +0.01224884*(function22[i][j](c,c,c) +function22[i][j](d,c,c) +function22[i][j](c,d,c) +function22[i][j](c,c,d))
			            +0.007091003*(function22[i][j](f,e,e) + function22[i][j](e,f,e) +function22[i][j](e,e,f) +function22[i][j](e,f,f)+ function22[i][j](f,e,f) +function22[i][j](f,f,e)));
		  }
	   }

   
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
    
    A23 function23[10][10] = {{A23_00,A23_01,A23_02,A23_03,A23_04,A23_05,A23_06,A23_07,A23_08,A23_09},
		                      {A23_10,A23_11,A23_12,A23_13,A23_14,A23_15,A23_16,A23_17,A23_18,A23_19},
		                      {A23_20,A23_21,A23_22,A23_23,A23_24,A23_25,A23_26,A23_27,A23_28,A23_29},
		                      {A23_30,A23_31,A23_32,A23_33,A23_34,A23_35,A23_36,A23_37,A23_38,A23_39},
		                      {A23_40,A23_41,A23_42,A23_43,A23_44,A23_45,A23_46,A23_47,A23_48,A23_49},
		                      {A23_50,A23_51,A23_52,A23_53,A23_54,A23_55,A23_56,A23_57,A23_58,A23_59},
		                      {A23_60,A23_61,A23_62,A23_63,A23_64,A23_65,A23_66,A23_67,A23_68,A23_69},
		                      {A23_70,A23_71,A23_72,A23_73,A23_74,A23_75,A23_76,A23_77,A23_78,A23_79},
		                      {A23_80,A23_81,A23_82,A23_83,A23_84,A23_85,A23_86,A23_87,A23_88,A23_89},
		                      {A23_90,A23_91,A23_92,A23_93,A23_94,A23_95,A23_96,A23_97,A23_98,A23_99}};
		                              
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function23[i][j](a,a,a) +function23[i][j](b,a,a) +function23[i][j](a,b,a) +function23[i][j](a,a,b)) +0.01224884*(function23[i][j](c,c,c) +function23[i][j](d,c,c) +function23[i][j](c,d,c) +function23[i][j](c,c,d))
			            +0.007091003*(function23[i][j](f,e,e) + function23[i][j](e,f,e) +function23[i][j](e,e,f) +function23[i][j](e,f,f)+ function23[i][j](f,e,f) +function23[i][j](f,f,e)));
		  }
	   }

    

}
//A30
typedef double (*A30) (double x, double y, double z);
typedef double (*B3) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double mu, double lambda, double kp, double r0_n, double r1_n, double r2_n, double r3_n, double r0_m, double r1_m, double r2_m, double r3_m, 
double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, double u6_n, double u7_n, double u8_n, double u9_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, 
double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, double v6_n, double v7_n, double v8_n, double v9_n, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, 
double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, double w6_n, double w7_n, double w8_n, double w9_n, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m );
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
    
    A30 function30[10][4] = {{A30_00,A30_01,A30_02,A30_03},
		                     {A30_10,A30_11,A30_12,A30_13},
		                     {A30_20,A30_21,A30_22,A30_23},
		                     {A30_30,A30_31,A30_32,A30_33},
		                     {A30_40,A30_41,A30_42,A30_43},
		                     {A30_50,A30_51,A30_52,A30_53},
		                     {A30_60,A30_61,A20_62,A30_63},
		                     {A30_70,A30_71,A30_72,A30_73},
		                     {A30_80,A30_81,A30_82,A30_83},
		                     {A30_90,A30_91,A30_92,A30_93}};
		                      
	 B3 function3[10] = {B3_0,B3_1,B3_2,B3_3,B3_4,B3_5,B3_6,B3_7,B3_8,B3_9};
	 
	 for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=3; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function30[i][j](a,a,a) +function30[i][j](b,a,a) +function30[i][j](a,b,a) +function30[i][j](a,a,b)) +0.01224884*(function30[i][j](c,c,c) +function30[i][j](d,c,c) +function30[i][j](c,d,c) +function30[i][j](c,c,d))
			            +0.007091003*(function30[i][j](f,e,e) + function30[i][j](e,f,e) +function30[i][j](e,e,f) +function30[i][j](e,f,f)+ function30[i][j](f,e,f) +function30[i][j](f,f,e)));
		  }
	   }
    
             
      for(int i=0; i<=9; ++i)
      {
			  fe[i]=jac*(0.01878132*(function3[i](a,a,a, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			      +function3[i](b,a,a, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			      +function3[i](a,b,a, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			      +function3[i](a,a,b, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m)) 
			      +0.01224884*(function3[i](c,c,c, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			      +function3[i](d,c,c, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			      +function3[i](c,d,c, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			      +function3[i](c,c,d, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m))
			      +0.007091003*(function3[i](f,e,e, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			      + function3[i](e,f,e, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			      +function3[i](e,e,f, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			      +function3[i](e,f,f, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m)
			      + function3[i](f,e,f, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			      +function3[i](f,f,e, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m)));
	   }

  
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
     
     A31 function31[10][10] = {{A31_00,A31_01,A31_02,A31_03,A31_04,A31_05,A31_06,A31_07,A31_08,A31_09},
		                       {A31_10,A31_11,A31_12,A31_13,A31_14,A31_15,A31_16,A31_17,A31_18,A31_19},
		                       {A31_20,A31_21,A31_22,A31_23,A31_24,A31_25,A31_26,A31_27,A31_28,A31_29},
		                       {A31_30,A31_31,A31_32,A31_33,A31_34,A31_35,A31_36,A31_37,A31_38,A31_39},
		                       {A31_40,A31_41,A31_42,A31_43,A31_44,A31_45,A31_46,A31_47,A31_48,A31_49},
		                       {A31_50,A31_51,A31_52,A31_53,A31_54,A31_55,A31_56,A31_57,A31_58,A31_59},
		                       {A31_60,A31_61,A31_62,A31_63,A31_64,A31_65,A31_66,A31_67,A31_68,A31_69},
		                       {A31_70,A31_71,A31_72,A31_73,A31_74,A31_75,A31_76,A31_77,A31_78,A31_79},
		                       {A31_80,A31_81,A31_82,A31_83,A31_84,A31_85,A31_86,A31_87,A31_88,A31_89},
		                       {A31_90,A31_91,A31_92,A31_93,A31_94,A31_95,A31_96,A31_97,A31_98,A31_99}};
		                                
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function31[i][j](a,a,a) +function31[i][j](b,a,a) +function31[i][j](a,b,a) +function31[i][j](a,a,b)) +0.01224884*(function31[i][j](c,c,c) +function31[i][j](d,c,c) +function31[i][j](c,d,c) +function31[i][j](c,c,d))
			            +0.007091003*(function31[i][j](f,e,e) + function31[i][j](e,f,e) +function31[i][j](e,e,f) +function31[i][j](e,f,f)+ function31[i][j](f,e,f) +function31[i][j](f,f,e)));
		  }
	   }

    
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
    
    A32 function32[10][10] = { {A32_00,A32_01,A32_02,A32_03,A32_04,A32_05,A32_06,A32_07,A32_08,A32_09},
		                       {A32_10,A32_11,A32_12,A32_13,A32_14,A32_15,A32_16,A32_17,A32_18,A32_19},
		                       {A32_20,A32_21,A32_22,A32_23,A32_24,A32_25,A32_26,A32_27,A32_28,A32_29},
		                       {A32_30,A32_31,A32_32,A32_33,A32_34,A32_35,A32_36,A32_37,A32_38,A32_39},
		                       {A32_40,A32_41,A32_42,A32_43,A32_44,A32_45,A32_46,A32_47,A32_48,A32_49},
		                       {A32_50,A32_51,A32_52,A32_53,A32_54,A32_55,A32_56,A32_57,A32_58,A32_59},
		                       {A32_60,A32_61,A32_62,A32_63,A32_64,A32_65,A32_66,A32_67,A32_68,A32_69},
		                       {A32_70,A32_71,A32_72,A32_73,A32_74,A32_75,A32_76,A32_77,A32_78,A32_79},
		                       {A32_80,A32_81,A32_82,A32_83,A32_84,A32_85,A32_86,A32_87,A32_88,A32_89},
		                       {A32_90,A32_91,A32_92,A32_93,A32_94,A32_95,A32_96,A32_97,A32_98,A32_99}};
		                                
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function32[i][j](a,a,a) +function32[i][j](b,a,a) +function32[i][j](a,b,a) +function32[i][j](a,a,b)) +0.01224884*(function32[i][j](c,c,c) +function32[i][j](d,c,c) +function32[i][j](c,d,c) +function32[i][j](c,c,d))
			            +0.007091003*(function32[i][j](f,e,e) + function32[i][j](e,f,e) +function32[i][j](e,e,f) +function32[i][j](e,f,f)+ function32[i][j](f,e,f) +function32[i][j](f,f,e)));
		  }
	   }

   
}
//A33
typedef double (*A33) (double x, double y, double z, double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3, double dt, double mu, double lambda, double kp, double r0_n, double r1_n, double r2_n, double r3_n, double r0_m, double r1_m, double r2_m, double r3_m, 
double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, double u6_n, double u7_n, double u8_n, double u9_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double u6_m, double u7_m, double u8_m, double u9_m, 
double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, double v6_n, double v7_n, double v8_n, double v9_n, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double v6_m, double v7_m, double v8_m, double v9_m, 
double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, double w6_n, double w7_n, double w8_n, double w9_n, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m, double w6_m, double w7_m, double w8_m, double w9_m);
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
     
     A33 function33[10][10] = {{A33_00,A33_01,A33_02,A33_03,A33_04,A33_05,A33_06,A33_07,A33_08,A33_09},
		                       {A33_10,A33_11,A33_12,A33_13,A33_14,A33_15,A33_16,A33_17,A33_18,A33_19},
		                       {A33_20,A33_21,A33_22,A33_23,A33_24,A33_25,A33_26,A33_27,A33_28,A33_29},
		                       {A33_30,A33_31,A33_32,A33_33,A33_34,A33_35,A33_36,A33_37,A33_38,A33_39},
		                       {A33_40,A33_41,A33_42,A33_43,A33_44,A33_45,A33_46,A33_47,A33_48,A33_49},
		                       {A33_50,A33_51,A33_52,A33_53,A33_54,A33_55,A33_56,A33_57,A33_58,A33_59},
		                       {A33_60,A33_61,A33_62,A33_63,A33_64,A33_65,A33_66,A33_67,A33_68,A33_69},
		                       {A33_70,A33_71,A33_72,A33_73,A33_74,A33_75,A33_76,A33_77,A33_78,A33_79},
		                       {A33_80,A33_81,A33_82,A33_83,A33_84,A33_85,A33_86,A33_87,A33_88,A33_89},
		                       {A33_90,A33_91,A33_92,A33_93,A33_94,A33_95,A33_96,A33_97,A33_98,A33_99}};
		                                
      for(int i=0; i<=9; ++i)
      {
		  for(int j=0; j<=9; ++j)
		  {
			  ske[i][j]=jac*(0.01878132*(function33[i][j](a,a,a, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			           +function33[i][j](b,a,a, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			           +function33[i][j](a,b,a, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			           +function33[i][j](a,a,b, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m)) 
			           +0.01224884*(function33[i][j](c,c,c, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			           +function33[i][j](d,c,c, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			           +function33[i][j](c,d,c, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			           +function33[i][j](c,c,d, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m))
			           +0.007091003*(function33[i][j](f,e,e, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			           +function33[i][j](e,f,e, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			           +function33[i][j](e,e,f, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			           +function33[i][j](e,f,f, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m)
			           +function33[i][j](f,e,f, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m) 
			           +function33[i][j](f,f,e, a1, a2, a3, b1, b2, b3, c1, c2, c3, dt, mu, lambda, kp, r0_n, r1_n, r2_n, r3_n, r0_m, r1_m, r2_m, r3_m, u0_n, u1_n, u2_n, u3_n, u4_n, u5_n, u6_n, u7_n, u8_n, u9_n, u0_m, u1_m, u2_m, u3_m, u4_m, u5_m, u6_m, u7_m, u8_m, u9_m, v0_n, v1_n, v2_n, v3_n, v4_n, v5_n, v6_n, v7_n, v8_n, v9_n, v0_m, v1_m, v2_m, v3_m, v4_m, v5_m, v6_m, v7_m, v8_m, v9_m, 
w0_n, w1_n, w2_n, w3_n, w4_n, w5_n, w6_n, w7_n, w8_n, w9_n, w0_m, w1_m, w2_m, w3_m, w4_m, w5_m, w6_m, w7_m, w8_m, w9_m)));
		  }
	   }

  
}*/
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



//CRS_Matrix::CRS_Matrix(Mesh const & mesh, int ndof_v)
 //: _mesh(mesh), _nrows(0), _nnz(0), _id(0), _ik(0), _sk(0), _ndof_v(ndof_v)
//{
    //Derive_Matrix_Pattern();
    //Skalar2VectorMatrix(ndof_v);
    //return;
//}

//void CRS_Matrix::Derive_Matrix_Pattern()
//{
    //int const nelem(_mesh.Nelems());
    //int const ndof_e(_mesh.NdofsElement());
    //auto const &ia(_mesh.GetConnectivity());
////  Determine the number of matrix rows
    //_nrows = *max_element(ia.cbegin(), ia.cbegin() + ndof_e * nelem);
    //++_nrows;                                 // node numberng: 0 ... nnode-1
    //assert(*min_element(ia.cbegin(), ia.cbegin() + ndof_e * nelem) == 0); // numbering starts with 0 ?

////  Collect for each node those nodes it is connected to (multiple entries)
////  Detect the neighboring nodes
    //vector< list<int> > cc(_nrows);             //  cc[i] is the  list of nodes a node i is connected to
    //for (int i = 0; i < nelem; ++i)
    //{
        //int const idx = ndof_e * i;
        //for (int k = 0; k < ndof_e; ++k)
        //{
            //list<int> &cck = cc.at(ia[idx + k]);
            //cck.insert( cck.end(), ia.cbegin() + idx, ia.cbegin() + idx + ndof_e );
        //}
    //}
////  Delete the multiple entries
    //_nnz = 0;
    //for (auto &it : cc)
    //{
        //it.sort();
        //it.unique();
        //_nnz += it.size();
        //// cout << it.size() << " :: "; copy(it->begin(),it->end(), ostream_iterator<int,char>(cout,"  ")); cout << endl;
    //}

//// CSR data allocation
    //_id.resize(_nrows + 1);                  // Allocate memory for CSR row pointer
    //_ik.resize(_nnz);                        // Allocate memory for CSR column index vector

////  copy CSR data
    //_id[0] = 0;                                 // begin of first row
    //for (size_t i = 0; i < cc.size(); ++i)
    //{
        ////cout << i << "   " << nid.at(i) << endl;;
        //const list<int> &ci = cc.at(i);
        //const auto nci = static_cast<int>(ci.size());
        //_id[i + 1] = _id[i] + nci; // begin of next line
        //copy(ci.begin(), ci.end(), _ik.begin() + _id[i] );
    //}

    //assert(_nnz == _id[_nrows]);
    //_sk.resize(_nnz);                      // Allocate memory for CSR column index vector
    //return;
//}

//void CRS_Matrix::Skalar2VectorMatrix(int ndof_v)
//{
    //this->Debug();
    //cout << "\n########################\n";
    //if (1 == ndof_v) return;
    //assert(2 == ndof_v);

    //auto old_id = _id;
    //auto old_ik = _ik;

    //_sk.resize(ndof_v * ndof_v * _sk.size(), -1.0);
    //_id.resize(ndof_v * (_id.size() - 1) + 1);
    //_ik.resize(ndof_v * ndof_v * _ik.size(), -7);

    //_id[0] = 0;
    //for (int kold = 0; kold < Nrows(); ++kold) {
        //int nr = old_id[kold + 1] - old_id[kold];
        
        //for (int ii=1; ii<=ndof_v; ++ii){
			//_id[ndof_v * kold + ii] = _id[ndof_v * kold + ii-1] + ndof_v * nr;
		//}
  
    //}

    //for (int newrow = 0; newrow < ndof_v * Nrows(); ++newrow ) {
        //int oldrow = newrow / ndof_v;
        //int idx = _id[newrow];
        ////cout << " ("<< newrow<<")  " << idx;
        //for (int oid = old_id[oldrow]; oid < old_id[oldrow + 1]; ++oid) {
            //int oldcol = old_ik[oid];
            
            //for (int ii=0; ii < ndof_v; ++ii){
				//_ik[idx + ii] = ndof_v * oldcol+ii;
				//}
         
            //idx += ndof_v;
        //}
    //}

    //_nrows =  _id.size() - 1;
    //_nnz   =  _ik.size();

    //return;
//}

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

//void CRS_Matrix::CalculateLaplace(vector<double> &f)
//{
    //assert(_mesh.NdofsElement() == 4);               // only for triangular, linear elements
    ////cout << _nnz << " vs. " << _id[_nrows] << "  " << _nrows<< endl;
    //assert(_nnz == _id[_nrows]);

    //for (int k = 0; k < _nrows; ++k)
    //{
        //_sk[k] = 0.0;
    //}
    //for (int k = 0; k < _nrows; ++k)
    //{
        //f[k] = 0.0;
    //}

    ////double ske[3][3], fe[3];               // move inside loop (==> thread private)
    ////  Loop over all elements
    //auto const nelem = _mesh.Nelems();
    //auto const &ia   = _mesh.GetConnectivity();
    //auto const &xc   = _mesh.GetCoords();

//#pragma omp parallel for
    //for (int i = 0; i < nelem; ++i)
    //{
        //double ske[4][4], fe[4];             // OpenMP: Thread private
        //CalcElem(ia.data()+4 * i, xc.data(), ske, fe);
        ////AddElem(ia.data()+3 * i, ske, fe, _id.data(), _ik.data(), _sk.data(), f.data()); // GH: deprecated
        //AddElem_3(ia.data()+4 * i, ske, fe, f);
    //}

    ////Debug();

    //return;
//}

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


//void CRS_Matrix::AddElem_3(int const ial[4], double const ske[4][4], double const fe[4], vector<double> & f)
//{
    //for (int i = 0; i < 4; ++i)
    //{
        //const int ii  = ial[i];           // row ii (global index)
        //for (int j = 0; j < 4; ++j)       // no symmetry assumed
        //{
            //const int jj = ial[j];        // column jj (global index)
            //int ip = fetch(ii,jj);        // find column entry jj in row ii
//#ifndef NDEBUG                 // compiler option -DNDEBUG switches off the check
            //if (ip<0)          // no entry found !!
            //{
                //cout << "Error in AddElem: (" << ii << "," << jj << ") ["
                     //<< ial[0] << "," << ial[1] << "," << ial[2] << "," << ial[3] << "]\n";
                //assert(ip>=0);
            //}
//#endif
//#pragma omp atomic
            //_sk[ip] += ske[i][j];
        //}
//#pragma omp atomic
        //f[ii] += fe[i];
    //}
//}
// ############################################################################################################################################################################################################


Matrix::Matrix(int const nrows, int const ncols)
    : _nrows(nrows), _ncols(ncols), _dd(0)
{}


Matrix::~Matrix()
{}



CRS_Matrix::CRS_Matrix()
    : Matrix(0, 0), _nnz(0), _id(0), _ik(0), _sk(0)
{}

//CRS_Matrix::CRS_Matrix(const std::string &file) : Matrix(0, 0), _nnz(0), _id(0), _ik(0), _sk(0)
//{
    //readBinary(file);
    //_nrows = static_cast<int>(size(_id) - 1);
    //_ncols = _nrows;
//}


CRS_Matrix::~CRS_Matrix()
{}


void CRS_Matrix::JacobiSmoother(std::vector<double> const &f, std::vector<double> &u,
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

void CRS_Matrix::writeBinary(const std::string &file)
{
    vector<int> cnt(size(_id) - 1);
    for (size_t k = 0; k < size(cnt); ++k) {
        cnt[k] = _id[k + 1] - _id[k];
    }
    //adjacent_difference( cbegin(_id)+1, cend(_id), cnt );
    write_binMatrix(file, cnt, _ik, _sk);
}

void CRS_Matrix::readBinary(const std::string &file)
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
    : CRS_Matrix(), _mesh(mesh), _ndof_v(ndof_v)
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

bool CRS_Matrix::CheckSymmetry() const
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


bool CRS_Matrix::CheckRowSum() const
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

bool CRS_Matrix::CheckMproperty() const
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

bool CRS_Matrix::ForceMproperty()
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

bool CRS_Matrix::CheckMatrix() const
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

void CRS_Matrix::GetDiag_M(vector<double> &d) const
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

void DefectRestrict(CRS_Matrix const &SK1, BisectInterpolation const &P,
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

