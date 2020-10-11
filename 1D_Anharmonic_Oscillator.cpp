#include <stdio.h>
#include <iostream>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <itpp/itbase.h>

#define PI 3.14159265358979323846264338

using namespace std;
using namespace itpp;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//

	double factorial (int y) //computes the factorial of an integer
	{	
		if (y>1) return (y*factorial(y-1));
		else return (1);
	}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//

	double herm (double x, int n) //Hermite polynomial generator; works properly
	{
		double Hn, Hn0, Hn1;

		Hn0 = 1.0;
		Hn1 = 2.0*x;
		Hn = 0.0;

		if (n==0) return Hn0;
		if (n==1) return Hn1;
		for (int i = 2; i<=n; i++)
			{
				Hn = 2*(x*Hn1 - (i-1)*Hn0);
			 	Hn0=Hn1;
				Hn1=Hn;
			}

			return Hn;
	}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//
	
	double phi (double x, int q) //H.O basis
	{
		return exp(-x*x/2.0)*herm(x,q)/sqrt(pow(2.0,q)*sqrt(PI)*factorial(q));
	}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//

	double D2 (double x, int n, double h) //Second derivative
	{
		return ((-phi(x+2.0*h, n))+(16.0*phi(x+h, n))-(30.0*phi(x, n))+(16.0*phi(x-h, n))-(phi(x-2.0*h, n)))/(12.0*(h*h));
	}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//
	
	double W (double x)
	{
		return pow(x,3.0);
	}
	
	double dW (double x)
	{
		return 3.0*x*x;
	}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//

	double V (double x) //Potential
	{
		return W(x)*W(x) - dW(x);
	}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//

	void Evals (const mat &A, vec &d, mat &V) //obtains the eigenvalues and eigenvectors of a given matrix
	{  
		eig_sym(A, d, V);	
	}	

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//

	struct InnerProduct_params //parameters
	{
		int m;
	  	int n; 
	};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//

	double InnerProduct (double x, void * params) //Function to integrate
	{  
		struct InnerProduct_params * p = (struct InnerProduct_params *) params; 
		int m = p->m;
		int n = p->n;
		double h = 0.001;
		return phi(x,m)*(-D2(x,n,h)+V(x)*phi(x,n));
	}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//


	double matrix_maker(int m, int n){

		double result1, error1;
		struct InnerProduct_params p = {m,n};

		gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);

		gsl_function F;
		F.function = &InnerProduct;
		F.params = &p;

		gsl_integration_qag (&F, -10, 10, 1e-3, 1e-5, 10000, 2, w, &result1, &error1);
		
		return (result1);
		
		gsl_integration_workspace_free(w);
	}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//

int main (){
	
	int N = 15;
	mat A1(N,N);
	mat V1(N,N);
	vec d1(N);
	vec d2(N);
	
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	{
		if (fabs(matrix_maker(i,j)) < 1e-6) A1(i,j) = 0.0;		
		else A1(i,j) = matrix_maker(i,j);
	}

	Evals(A1,d1,V1);

	for (int k = 0; k < 5; k++){

		cout << d1[k];
		cout << endl;
	}
	
	

	return 0;



}








