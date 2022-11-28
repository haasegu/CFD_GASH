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

typedef double (*A00) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m,  
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m);

typedef double (*A01) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m);

typedef double (*A02) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m);

typedef double (*A10) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m, double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m);

typedef double (*A11) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m);

typedef double (*A12) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m);

typedef double (*A20) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double a, double gamma, double r0_m, double r1_m, double r2_m, double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m);

typedef double (*A21) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m);

typedef double (*A22) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m);

typedef double (*A30) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m,
double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m);

typedef double (*A31) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m);

typedef double (*A32) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_m, double r1_m, double r2_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m);

typedef double (*A33) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m);

typedef double (*B0) (double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double r0_n, double r1_n, double r2_n, double r0_m, double r1_m, double r2_m,
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m);

typedef double (*B1) (double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double u0_n, double u1_n, double u2_n, double u3_n, double u4_n, double u5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m);

typedef double (*B2) (double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double mu, double lambda, double a, double gamma, double r0_m, double r1_m, double r2_m,
double v0_n, double v1_n, double v2_n, double v3_n, double v4_n, double v5_n, double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m);

typedef double (*B3) (double t_ni, double x, double  y, double a1, double a2, double a3, double b1, double b2, double b3, double dt, double kt, double cp, double r0_m, double r1_m, double r2_m, double w0_n, double w1_n, double w2_n, double w3_n, double w4_n, double w5_n, 
double u0_m, double u1_m, double u2_m, double u3_m, double u4_m, double u5_m, double v0_m, double v1_m, double v2_m, double v3_m, double v4_m, double v5_m, double w0_m, double w1_m, double w2_m, double w3_m, double w4_m, double w5_m);


void CalcElem_Navier_Stokes(int const ial[], double const xc[], std::vector<std::vector<double>> &ske, std::vector<double> &fe, const std::vector<double> &r_old_n, const std::vector<double> &r_old_m, const std::vector<double> &u_old_n, const std::vector<double> &u_old_m, const std::vector<double> &v_old_n, const std::vector<double> &v_old_m, 
const std::vector<double> &w_old_n, const std::vector<double> &w_old_m, const double dt, const double mu, const double lambda, const double a, const double gamma, const double t_ni, const double kt, const double cp)
{
    const int    i1  = 2 * ial[0],   i2 = 2 * ial[1],   i3 = 2 * ial[2];
    const double a1 = xc[i1 + 0],  b1 = xc[i1 + 1],
                 a2 = xc[i2 + 0],  b2 = xc[i2 + 1],
                 a3 = xc[i3 + 0],  b3 = xc[i3 + 1];
             
     //cout<<"(x1,y1,z1) = ("<<xc[i1 + 0]<<","<<xc[i1 + 1]<<","<<xc[i1 + 2]<<")"<<"  (x2,y2,z2) = ("<<xc[i2 + 0]<<","<<xc[i2 + 1]<<","<<xc[i2 + 2]<<")"<<"   (x3,y3,z3) = ("<<xc[i3 + 0]<<","<<xc[i3 + 1]<<","<<xc[i3 + 2]<<")"<<"    (x4,y4,z4) = ("<<xc[i4 + 0]<<","<<xc[i4 + 1]<<","<<xc[i4 + 2]<<")"<<endl;
                 
    const double jac1 = a1*b2 - a2*b1 - a1*b3 + a3*b1 + a2*b3 - a3*b2;
    const double jac = std::abs(jac1);
    
    const double r0_n = r_old_n.at(ial[0]), r1_n = r_old_n.at(ial[1]), r2_n = r_old_n.at(ial[2]), 
                 r0_m = r_old_m.at(ial[0]), r1_m = r_old_m.at(ial[1]), r2_m = r_old_m.at(ial[2]);
                            
     const double u0_n = u_old_n.at(ial[0]), u1_n = u_old_n.at(ial[1]), u2_n = u_old_n.at(ial[2]), u3_n = u_old_n.at(ial[3]), u4_n = u_old_n.at(ial[4]), u5_n = u_old_n.at(ial[5]), 
                  u0_m = u_old_m.at(ial[0]), u1_m = u_old_m.at(ial[1]), u2_m = u_old_m.at(ial[2]), u3_m = u_old_m.at(ial[3]), u4_m = u_old_m.at(ial[4]), u5_m = u_old_m.at(ial[5]), 
                  v0_n = v_old_n.at(ial[0]), v1_n = v_old_n.at(ial[1]), v2_n = v_old_n.at(ial[2]), v3_n = v_old_n.at(ial[3]), v4_n = v_old_n.at(ial[4]), v5_n = v_old_n.at(ial[5]), 
                  v0_m = v_old_m.at(ial[0]), v1_m = v_old_m.at(ial[1]), v2_m = v_old_m.at(ial[2]), v3_m = v_old_m.at(ial[3]), v4_m = v_old_m.at(ial[4]), v5_m = v_old_m.at(ial[5]), 
                  w0_n = w_old_n.at(ial[0]), w1_n = w_old_n.at(ial[1]), w2_n = w_old_n.at(ial[2]), w3_n = w_old_n.at(ial[3]), w4_n = w_old_n.at(ial[4]), w5_n = w_old_n.at(ial[5]), 
                  w0_m = w_old_m.at(ial[0]), w1_m = w_old_m.at(ial[1]), w2_m = w_old_m.at(ial[2]), w3_m = w_old_m.at(ial[3]), w4_m = w_old_m.at(ial[4]), w5_m = w_old_m.at(ial[5]); 
                  
     
     
	//double aa=0.3108859, b=1-3*a, c=0.09273525, d=1-3*c, e=0.454463, f=0.5-e;
	
	//const double  x = (a1-a3)*\xi+(a2-a3)*\eta+a3, y = (b1-b3)*\xi+(b2-b3)*\eta+b3; 
	
	const double  x1 = (a1-a3)*(0.0)+(a2-a3)*(0.0)+a3, y1 = (b1-b3)*(0.0)+(b2-b3)*(0.0)+b3;
	const double  x2 = (a1-a3)*(1.0)+(a2-a3)*(0.0)+a3, y2 = (b1-b3)*(1.0)+(b2-b3)*(0.0)+b3;
	const double  x3 = (a1-a3)*(0.0)+(a2-a3)*(1.0)+a3, y3 = (b1-b3)*(0.0)+(b2-b3)*(1.0)+b3;
	const double  x4 = (a1-a3)*(0.0)+(a2-a3)*(0.5)+a3, y4 = (b1-b3)*(0.0)+(b2-b3)*(0.5)+b3;
	const double  x5 = (a1-a3)*(0.5)+(a2-a3)*(0.0)+a3, y5 = (b1-b3)*(0.5)+(b2-b3)*(0.0)+b3;
	const double  x6 = (a1-a3)*(0.5)+(a2-a3)*(0.5)+a3, y6 = (b1-b3)*(0.5)+(b2-b3)*(0.5)+b3;
	const double  x7 = (a1-a3)*(1.0/3)+(a2-a3)*(1.0/3)+a3, y7 = (b1-b3)*(1.0/3)+(b2-b3)*(1.0/3)+b3;
	
	
	
	
	
     A00 function00[3][3] = {{A00_00,A00_01,A00_02},
		                     {A00_10,A00_11,A00_12},
		                     {A00_20,A00_21,A00_22}};
		                     
       
      for(int i=0; i<=2; ++i)
      {
		  for(int j=0; j<=2; ++j)
		  {
			  ske[0 + (1-1)*6+i][0 + (1-1)*6+j]=jac*((3.0/120)*(function00[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function00[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function00[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(8.0/120)*(function00[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function00[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function00[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(27.0/120)*function00[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m));
		  }
	   }
	   
      A01 function01[3][6] =  {{A01_00,A01_01,A01_02,A01_03,A01_04,A01_05},
		                       {A01_10,A01_11,A01_12,A01_13,A01_14,A01_15},
		                       {A01_20,A01_21,A01_22,A01_23,A01_24,A01_25}};
     
              
      for(int i=0; i<=2; ++i)
      {
		  for(int j=0; j<=5; ++j)
		  {
			  ske[0 + (1-1)*6+i][3 + (1-1)*6+j]=jac*((3.0/120)*(function01[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m) 
			                                                   +function01[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m) 
			                                                   +function01[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m))
			                                        +(8.0/120)*(function01[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m) 
			                                                   +function01[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m) 
			                                                   +function01[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m))
			                                        +(27.0/120)*function01[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m));
		  }
	   }
  
     A02 function02[3][6] =  { {A02_00,A02_01,A02_02,A02_03,A02_04,A02_05},
		                       {A02_10,A02_11,A02_12,A02_13,A02_14,A02_15},
		                       {A02_20,A02_21,A02_22,A02_23,A02_24,A02_25}};
              
      for(int i=0; i<=2; ++i)
      {
		  for(int j=0; j<=5; ++j)
		  {
			  ske[0 + (1-1)*6+i][3 + (2-1)*6+j]=jac*((3.0/120)*(function02[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m) 
			                                                   +function02[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m) 
			                                                   +function02[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m))
			                                        +(8.0/120)*(function02[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m) 
			                                                   +function02[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m) 
			                                                   +function02[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m))
			                                        +(27.0/120)*function02[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m));
		  }
	   }
	   
	 for(int i=0; i<=2; ++i)
      {
		  for(int j=0; j<=5; ++j)
		  {
			  ske[0 + (1-1)*6+i][3 + (3-1)*6+j]=0;
		  }
	   }
	   
	 A10 function10[10][4] = {{A10_00,A10_01,A10_02},
		                      {A10_10,A10_11,A10_12},
		                      {A10_20,A10_21,A10_22},
		                      {A10_30,A10_31,A10_32},
		                      {A10_40,A10_41,A10_42},
		                      {A10_50,A10_51,A10_52}};
	 
	 for(int i=0; i<=5; ++i)
      {
		  for(int j=0; j<=2; ++j)
		  {
			  ske[3 + (1-1)*6+i][0 + (1-1)*6+j]=jac*((3.0/120)*(function10[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function10[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function10[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(8.0/120)*(function10[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function10[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function10[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(27.0/120)*function10[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m));
		  }
	   }
	   
	   
    A11 function11[6][6] =   {{A11_00,A11_01,A11_02,A11_03,A11_04,A11_05},
		                      {A11_10,A11_11,A11_12,A11_13,A11_14,A11_15},
		                      {A11_20,A11_21,A11_22,A11_23,A11_24,A11_25},
		                      {A11_30,A11_31,A11_32,A11_33,A11_34,A11_35},
		                      {A11_40,A11_41,A11_42,A11_43,A11_44,A11_45},
		                      {A11_50,A11_51,A11_52,A11_53,A11_54,A11_55}};
  
       
      for(int i=0; i<=5; ++i)
      {
		  for(int j=0; j<=5; ++j)
		  {
			  ske[3 + (1-1)*6+i][3 + (1-1)*6+j]=jac*((3.0/120)*(function11[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function11[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function11[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(8.0/120)*(function11[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function11[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function11[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(27.0/120)*function11[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m));
		  }
	   }
	   
	A12 function12[6][6] =   {{A12_00,A12_01,A12_02,A12_03,A12_04,A12_05},
		                      {A12_10,A12_11,A12_12,A12_13,A12_14,A12_15},
		                      {A12_20,A12_21,A12_22,A12_23,A12_24,A12_25},
		                      {A12_30,A12_31,A12_32,A12_33,A12_34,A12_35},
		                      {A12_40,A12_41,A12_42,A12_43,A12_44,A12_45},
		                      {A12_50,A12_51,A12_52,A12_53,A12_54,A12_55}};
  
      for(int i=0; i<=5; ++i)
      {
		  for(int j=0; j<=5; ++j)
		  {
			  ske[3 + (1-1)*6+i][3 + (2-1)*6+j]=jac*((3.0/120)*(function12[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m) 
			                                                   +function12[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m) 
			                                                   +function12[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m))
			                                        +(8.0/120)*(function12[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m) 
			                                                   +function12[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m) 
			                                                   +function12[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m))
			                                        +(27.0/120)*function12[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m));
		  }
	   }
	   
	 for(int i=0; i<=5; ++i)
      {
		  for(int j=0; j<=5; ++j)
		  {
			  ske[3 + (1-1)*6+i][3 + (2-1)*6+j]=0;
		  }
	   }
	   
	 A20 function20[10][4] = {{A20_00,A20_01,A20_02},
		                      {A20_10,A20_11,A20_12},
		                      {A20_20,A20_21,A20_22},
		                      {A20_30,A20_31,A20_32},
		                      {A20_40,A20_41,A20_42},
		                      {A20_50,A20_51,A20_52}};
	 
	 for(int i=0; i<=5; ++i)
      {
		  for(int j=0; j<=2; ++j)
		  {
			  ske[3 + (2-1)*6+i][0 + (1-1)*6+j]=jac*((3.0/120)*(function20[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function20[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function20[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(8.0/120)*(function20[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function20[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function20[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(27.0/120)*function20[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m));
		  }
	   }
	   
	A21 function21[6][6] =   {{A21_00,A21_01,A21_02,A21_03,A21_04,A21_05},
		                      {A21_10,A21_11,A21_12,A21_13,A21_14,A21_15},
		                      {A21_20,A21_21,A21_22,A21_23,A21_24,A21_25},
		                      {A21_30,A21_31,A21_32,A21_33,A21_34,A21_35},
		                      {A21_40,A21_41,A21_42,A21_43,A21_44,A21_45},
		                      {A21_50,A21_51,A21_52,A21_53,A21_54,A21_55}};
   
    
	                                 
      for(int i=0; i<=5; ++i)
      {
		  for(int j=0; j<=5; ++j)
		  {
			  ske[3 + (2-1)*6+i][3 + (1-1)*6+j]=jac*((3.0/120)*(function21[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function21[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function21[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(8.0/120)*(function21[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function21[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function21[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(27.0/120)*function21[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m));
		  }
	   }

   
    A22 function22[6][6] =   {{A22_00,A22_01,A22_02,A22_03,A22_04,A22_05},
		                      {A22_10,A22_11,A22_12,A22_13,A22_14,A22_15},
		                      {A22_20,A22_21,A22_22,A22_23,A22_24,A22_25},
		                      {A22_30,A22_31,A22_32,A22_33,A22_34,A22_35},
		                      {A22_40,A22_41,A22_42,A22_43,A22_44,A22_45},
		                      {A22_50,A22_51,A22_52,A22_53,A22_54,A22_55}};
	                            
      for(int i=0; i<=5; ++i)
      {
		  for(int j=0; j<=5; ++j)
		  {
			  ske[3 + (2-1)*6+i][3 + (2-1)*6+j]=jac*((3.0/120)*(function22[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function22[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function22[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(8.0/120)*(function22[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function22[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function22[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(27.0/120)*function22[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,mu,lambda,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m));
		  }
	   }
	   
	 for(int i=0; i<=5; ++i)
      {
		  for(int j=0; j<=5; ++j)
		  {
			  ske[3 + (2-1)*6+i][3 + (3-1)*6+j]=0;
		  }
	   }
	   
	   
	  A30 function30[10][4] ={{A30_00,A30_01,A30_02},
		                      {A30_10,A30_11,A30_12},
		                      {A30_20,A30_21,A30_22},
		                      {A30_30,A30_31,A30_32},
		                      {A30_40,A30_41,A30_42},
		                      {A30_50,A30_51,A30_52}};
	 
	 for(int i=0; i<=5; ++i)
      {
		  for(int j=0; j<=2; ++j)
		  {
			  ske[3 + (3-1)*6+i][0 + (1-1)*6+j]=jac*((3.0/120)*(function30[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                                   +function30[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                                   +function30[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m))
			                                        +(8.0/120)*(function30[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                                   +function30[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                                   +function30[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m))
			                                        +(27.0/120)*function30[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m));
		  }
	   }
	   
	A31 function31[6][6] =   {{A31_00,A31_01,A31_02,A31_03,A31_04,A31_05},
		                      {A31_10,A31_11,A31_12,A31_13,A31_14,A31_15},
		                      {A31_20,A31_21,A31_22,A31_23,A31_24,A31_25},
		                      {A31_30,A31_31,A31_32,A31_33,A31_34,A31_35},
		                      {A31_40,A31_41,A31_42,A31_43,A31_44,A31_45},
		                      {A31_50,A31_51,A31_52,A31_53,A31_54,A31_55}};
   
    
	                                 
      for(int i=0; i<=5; ++i)
      {
		  for(int j=0; j<=5; ++j)
		  {
			  ske[3 + (3-1)*6+i][3 + (1-1)*6+j]=jac*((3.0/120)*(function31[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                                   +function31[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                                   +function31[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m))
			                                        +(8.0/120)*(function31[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                                   +function31[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                                   +function31[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m))
			                                        +(27.0/120)*function31[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m));
		  }
	   }

   
    A32 function32[6][6] =   {{A32_00,A32_01,A32_02,A32_03,A32_04,A32_05},
		                      {A32_10,A32_11,A32_12,A32_13,A32_14,A32_15},
		                      {A32_20,A32_21,A32_22,A32_23,A32_24,A32_25},
		                      {A32_30,A32_31,A32_32,A32_33,A32_34,A32_35},
		                      {A32_40,A32_41,A32_42,A32_43,A32_44,A32_45},
		                      {A32_50,A32_51,A32_52,A32_53,A32_54,A32_55}};
	                            
      for(int i=0; i<=5; ++i)
      {
		  for(int j=0; j<=5; ++j)
		  {
			  ske[3 + (3-1)*6+i][3 + (2-1)*6+j]=jac*((3.0/120)*(function32[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                                   +function32[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                                   +function32[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m))
			                                        +(8.0/120)*(function32[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                                   +function32[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                                   +function32[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m))
			                                        +(27.0/120)*function32[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,r0_m,r1_m,r2_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m));
		  }
	   }
	   
	 A33 function33[6][6] =  {{A33_00,A33_01,A33_02,A33_03,A33_04,A33_05},
		                      {A33_10,A33_11,A33_12,A33_13,A33_14,A33_15},
		                      {A33_20,A33_21,A33_22,A33_23,A33_24,A33_25},
		                      {A33_30,A33_31,A33_32,A33_33,A33_34,A33_35},
		                      {A33_40,A33_41,A33_42,A33_43,A33_44,A33_45},
		                      {A33_50,A33_51,A33_52,A33_53,A33_54,A33_55}};
	   
	 for(int i=0; i<=5; ++i)
      {
		  for(int j=0; j<=5; ++j)
		  {
			  ske[3 + (3-1)*6+i][3 + (3-1)*6+j]=jac*((3.0/120)*(function33[i][j](x1,y1,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function33[i][j](x2,y2,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function33[i][j](x3,y3,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(8.0/120)*(function33[i][j](x4,y4,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function33[i][j](x5,y5,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                                   +function33[i][j](x6,y6,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                                        +(27.0/120)*function33[i][j](x7,y7,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m));
		  }
	   } 
	   
	   
	   
	   
	B0 function0[3] = {B0_0,B0_1,B0_2};
       
	for(int i=0; i<=2; ++i)
      {
			  fe[0 + (1-1)*6+i]=jac*((3.0/120)*(function0[i](x1,y1,a1,a2,a3,b1,b2,b3,dt,r0_n,r1_n,r2_n,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                   +function0[i](x2,y2,a1,a2,a3,b1,b2,b3,dt,r0_n,r1_n,r2_n,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                   +function0[i](x3,y3,a1,a2,a3,b1,b2,b3,dt,r0_n,r1_n,r2_n,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                        +(8.0/120)*(function0[i](x4,y4,a1,a2,a3,b1,b2,b3,dt,r0_n,r1_n,r2_n,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                   +function0[i](x5,y5,a1,a2,a3,b1,b2,b3,dt,r0_n,r1_n,r2_n,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                   +function0[i](x6,y6,a1,a2,a3,b1,b2,b3,dt,r0_n,r1_n,r2_n,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                        +(27.0/120)*function0[i](x7,y7,a1,a2,a3,b1,b2,b3,dt,r0_n,r1_n,r2_n,r0_m,r1_m,r2_m,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m));
	   } 
	   
	 B1 function1[6] = {B1_0,B1_1,B1_2,B1_3,B1_4,B1_5};  
	   
	 for(int i=0; i<=5; ++i)
      {
			  fe[3 + (1-1)*6+i]= jac*((3.0/120)*(function1[i](t_ni,x1,y1,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                    +function1[i](t_ni,x2,y2,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                    +function1[i](t_ni,x3,y3,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                         +(8.0/120)*(function1[i](t_ni,x4,y4,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                    +function1[i](t_ni,x5,y5,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                    +function1[i](t_ni,x6,y6,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                         +(27.0/120)*function1[i](t_ni,x7,y7,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,u0_n,u1_n,u2_n,u3_n,u4_n,u5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m));
	   } 
	   
	 B2 function2[6] = {B2_0,B2_1,B2_2,B2_3,B2_4,B2_5};
	   
	 for(int i=0; i<=5; ++i)
      {
			  fe[3 + (2-1)*6+i]=jac*((3.0/120)*(function2[i](t_ni,x1,y1,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                   +function2[i](t_ni,x2,y2,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                   +function2[i](t_ni,x3,y3,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                        +(8.0/120)*(function2[i](t_ni,x4,y4,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                   +function2[i](t_ni,x5,y5,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m) 
			                                   +function2[i](t_ni,x6,y6,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m))
			                        +(27.0/120)*function2[i](t_ni,x7,y7,a1,a2,a3,b1,b2,b3,dt,mu,lambda,a,gamma,r0_m,r1_m,r2_m,v0_n,v1_n,v2_n,v3_n,v4_n,v5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m));
	   }
	  
	 B3 function3[6] = {B3_0,B3_1,B3_2,B3_3,B3_4,B3_5};
	   
	 for(int i=0; i<=5; ++i)
      {
			  fe[3 + (3-1)*6+i]=jac*((3.0/120)*(function3[i](t_ni,x1,y1,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                   +function3[i](t_ni,x2,y2,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                   +function3[i](t_ni,x3,y3,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m))
			                        +(8.0/120)*(function3[i](t_ni,x4,y4,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                   +function3[i](t_ni,x5,y5,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m) 
			                                   +function3[i](t_ni,x6,y6,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m))
			                        +(27.0/120)*function3[i](t_ni,x7,y7,a1,a2,a3,b1,b2,b3,dt,kt,cp,r0_m,r1_m,r2_m,w0_n,w1_n,w2_n,w3_n,w4_n,w5_n,u0_m,u1_m,u2_m,u3_m,u4_m,u5_m,v0_m,v1_m,v2_m,v3_m,v4_m,v5_m,w0_m,w1_m,w2_m,w3_m,w4_m,w5_m));
	   }

  
}


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
    cout << "_nrows " << _nrows << "      d  " << d.size() << endl;
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

/*
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
}*/
/*
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
}*/
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

void FEM_Matrix::ApplyDirichletBC_Box1(std::vector<double> const &u, std::vector<double> &f,
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

