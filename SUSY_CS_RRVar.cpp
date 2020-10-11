//
//  RRvarCoherent.cpp
//  DCSB
//
//  Created by Mason Biamonte on 7/4/11.
//  Copyright 2011 University of Houston. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <complex>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <HeaderFiles/GaussianExpansion.h>

#define PI 3.14159265358979323846264338


using namespace std;
using namespace Eigen;

typedef complex<double> dcomp;

dcomp I (0.0,1.0);
double n = 0.7511255120;
//double D = 22.5*2.0*PI;                 //density of gridpoints in phase space.
double A = 1.0;                           //scaling parameter.
const int M_tilde = 14;
const int N1 = 1;
const int N2 = 1;
const int K1 = 2*N1+1;
const int K2 = 2*N2+1;
const int M = K1*K2;                      //number of basis functions (X-pattern).
Vector2d a[M];                           //array of vectors storing (q,k) values.

double coeff[M_tilde];
double m[M_tilde];


dcomp EI (double x) //Euler identity.
{
	dcomp z(cos(x),sin(x));
	return z;
}

//______________________________SUSY COHERENT STATES_____________________________________________________________________________//

double W_st(double x) //Superpotential used to construct coherent states.
{
	return x*x*x;
} 

double dW_st(double x) //Derivative of superpotential above.
{
	return 3.0*x*x;
}

double int_W_st(double x) //Integral of superpotential above.
{
	return x*x*x*x/4.0;
}

double V_st(double x) //Potential corresponding to superpotential above.
{
	return W_st(x)*W_st(x)-dW_st(x);
}

//_______________________________________________________________________________________________________________________________*/

double V_sys(double x) //Potential of system under investigation.
{
	return pow(x,6.0)-3.0*x*x;
}

//_______________________________________________________________________________________________________________________________//

double FWHMRoots(double x, void * params)
{
    struct norm_params * p = 
    (struct norm_params *) params;
    double B = A;
    B = p->C;
    
    return A*int_W_st(x) + log(0.5);
}

//-------------------------------------------------------------------------------//

double uncertainty() //Root-solver computes the FWHM separation distance. (*)
{
    int status;
    int iter = 0, max_iter = 500;
    
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    
    double r = 0.0;
    double x_lo = 0.0, x_hi = 10.0;
    double C = A;
    
    gsl_function F;
    norm_params pidgeon = {C};
    
    F.function = &FWHMRoots;
    F.params = &pidgeon;
    
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    
    do{
        
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0.0, 1e-06);
        
    } while(status == GSL_CONTINUE && iter < max_iter);
    
    gsl_root_fsolver_free(s);
    
    return 2.0*r;	
}

//-------------------------------------------------------------------------------//

double norm_func(double x, void * params) //Provides integrand to compute normalization constant.
{
    struct norm_params * r = 
    (struct norm_params *) params;
    A = r->C;
    
    return exp(-2.0*A*int_W_st(x));
}

//-------------------------------------------------------------------------------//

double normalization() //Computes normalization constant. (*)
{
    double result, error;
    //double A = 1.0; //scaling parameter
    norm_params r = {A};
    
    gsl_integration_workspace * w = 
    gsl_integration_workspace_alloc(10000);
    
    gsl_function N;
    N.function = &norm_func;
    N.params = &r;
    
    gsl_integration_qag 
    (&N, -100.0, 100.0, 1, 0, 10000, 6, w, &result, &error);
    
    return 1.0/(sqrt(result));
    
    gsl_integration_workspace_free(w);
}

//_______________________________________________________________________________________________________________________________//

double q(int j, double D) //expvals of position on the phase-space grid (*)
{
	double delta_x = uncertainty();
    return j*delta_x*sqrt(2.0*PI/D);
}

double k(int i, double D) //expvals of momentum on the phase-space grid (*)
{
	double delta_x = uncertainty();
    return (i/(delta_x))*sqrt(2.0*PI/D);
}

//_______________________________________________________________________________________________________________________________//

void populate_lattice(double D) //fills array a[] with appropriate expectation values of position and momentum
{
	int c = 0;
    
    //---------------------------------SQUARE GRID-----------------------------//
    
    for(int i = N1; i >= -N1; i--)
    for(int j = -N2; j <= N2; j++)
    {
            a[c](0) = q(j,D);
            a[c](1) = k(i,D);
            c++;
    }

    //-------------------------------------------------------------------------*/

}

void populate_params(double array1[M_tilde], double array2[M_tilde])
{
    params(array1,array2);
}

double sigma (int i, int j)
{
    populate_params(coeff,m);
    return 1/(m[i]+m[j]);
}

double delta_k (int i, int j, double D)
{
    populate_lattice(D);
    return a[i](1)-a[j](1);
}

double rho_delta (int i, int j, double D)
{
    populate_lattice(D);
    return a[i](1)*a[i](0)-a[j](1)*a[j](0);
}

double beta (int i, int j, int r, int s, double D)
{
    populate_lattice(D);
    populate_params(coeff,m);
    
    return m[r]*a[i](0)+m[s]*a[j](0);
}

double gamma (int i, int j, int r, int s, double D)
{
    populate_lattice(D);
    populate_params(coeff,m);
    
    return m[r]*a[i](0)*a[i](0)+m[s]*a[j](0)*a[j](0);
}

double h2 (double (*f)(int,int), double (*g)(int,int), int i, int j, int r, int s)
{
    return (*f)(i,j)*(*f)(i,j)*(*g)(r,s)*
    (*g)(r,s)-2.0*((*g)(r,s)-1.0);
}

double h4 (double (*f)(int,int), double (*g)(int,int), int i, int j, int r, int s)
{
    return (*g)(r,s)*(*g)(r,s)*( 
    pow((*f)(i,j),4.0)*pow((*g)(r,s),2.0) - 
    12.0*pow((*f)(i,j),2.0)*((*g)(r,s)-1.0) 
    + 12.0 ) - 24.0*(*g)(r,s) + 12.0;
}

double h6 (double (*f)(int,int), double (*g)(int,int), int i, int j, int r, int s)
{
    return pow((*f)(i,j)*(*g)(r,s),6.0) 
    - 120.0*pow((*g)(r,s)-1.0,3.0) + 
    180.0*pow((*f)(i,j)*(*g)(r,s)*
    ((*g)(r,s)-1.0),2.0) - 
    30.0*pow((*f)(i,j)*(*g)(i,j),4.0)
    *((*g)(r,s)-1.0);
}

double f1 (double (*f)(int,int), double (*g)(int,int), double (*h)(int,int,int,int), int i, int j, int r, int s)
{
    return (15.0/8.0)+(45.0/16.0)*( 
    pow(2.0*(*h)(i,j,r,s)*(*g)(r,s), 2.0)
    - h2(&f,&g,i,j,r,s) ) + (15.0/32.0)*
    ( pow(2.0*(*h)(i,j,r,s)*(*g)(r,s), 4.0) 
    - 6.0*pow(2.0*(*h)(i,j,r,s)*(*g)(r,s), 2.0)
     *h2(&f,&g,i,j,r,s) + h4(&f,&g,i,j,r,s) ) + 
    (1.0/64.0)*(pow(2.0*(*h)(i,j,r,s)*(*g)(r,s), 6.0)
    - 15.0*pow(2.0*(*h)(i,j,r,s)*(*g)(r,s), 4.0)*
    h2(&f,&g,i,j,r,s) + 15.0*pow(2.0*(*h)(i,j,r,s)*
    (*g)(r,s), 2.0)*h4(&f,&g,i,j,r,s) - 
    h6(&f,&g,i,j,r,s) );
}

double f2 (double (*f)(int,int), double (*g)(int,int), double (*h)(int,int,int,int), int i, int j, int r, int s)
{
    return 0.25 + (1.0/8.0)*( pow(2.0*(*h)(i,j,r,s)*(*g)(r,s),2.0) - h2(&f,&g,i,j,r,s) );
}

dcomp eta_k (int s, int j, double D)
{
    populate_lattice(D);
    populate_params(coeff,m);
    
    return 4.0*m[s]*a[j](0)*a[j](0)+
    4*I*m[s]*a[j](1)*a[j](0)-2.0*m[s]-a[j](1)*a[j](1);
}

dcomp overlap (int i, int j, double D)
{
    double N = normalization();
    double sum = 0.0;
   
    populate_params(coeff,m);
    populate_lattice(D);
    
    for(int r = 0; r < M_tilde; r++)
    for(int s = 0; s < M_tilde; s++)
    {
        sum += coeff[r]*coeff[s]*cos(rho_delta(i,j,D))*
        cos(delta_k(i,j,D)*beta(i,j,r,s,D)*sigma(r,s))*
        exp(sigma(r,s)*beta(i,j,r,s,D)-gamma(i,j,r,s,D))*
        sqrt(PI*sigma(r,s))*exp(-0.25*delta_k(i,j,D)*d
        elta_k(i,j,D)*sigma(r,s);
    }
                                                                                                                                                                     
    return sum;
                                                                                                                                                                     
}

dcomp hamiltonian (int i, int j, double D)
{
    double N = normalization();
    double sum = 0.0;
   
    populate_params(coeff,m);
    populate_lattice(D);
    
    for(int r = 0; r < M_tilde; r++)
    for(int s = 0; s < M_tilde; s++)
    {
        sum += coeff[r]*coeff[s]*cos(rho_delta(i,j,D))*
        cos(delta_k(i,j,D)*beta(i,j,r,s,D)*sigma(r,s))*
        exp(sigma(r,s)*beta(i,j,r,s,D)-gamma(i,j,r,s,D))*
        sqrt(PI*sigma(r,s))*exp( -0.25*delta_k(i,j,D)*
        delta_k(i,j,D)*sigma(r,s) )*( f1(&delta_k,&sigma,&beta,i,j,r,s) + (3.0-4.0*m[s])*f2(&delta_k,&sigma,&beta,i,j,r,s) - eta_k(s,j,D) );
    }
                                                                                                                                                                     
    return sum;

}
                                                                                                                                                                

//-------------------------------------------------------------------------------------------------------------------------//

int main(){
    
    clock_t start, end;
    
    start = clock();
    
    MatrixXcd H(M,M);
    MatrixXcd S(M,M);
    VectorXcd v[M];
    
    double D = 190.0*2.0*PI;
    populate_lattice(D);
    populate_params();
    
    for(int i = 0; i < M; i++)
    for(int j = 0; j < M; j++)
    {
          S(i,j) = overlap(i,j,D);
          H(i,j) = hamiltonian(i,j,D);
    }
    
    cout << H << endl;
          
    end = clock();
    
    cout << endl;
    
    cout << "Running Time: " << (double)(end-start)/CLOCKS_PER_SEC << " s" << endl;
    
   

}

