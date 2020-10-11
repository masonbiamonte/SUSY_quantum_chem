#include <iostream>
#include <itpp/itbase.h>

using namespace itpp;
using namespace std;

	double kron (int i, int j) 
	{
			if (i==j) return i*i + 0.5;
			else if ((i==j+1) || (i==j-1)) return -0.5;
			else if ((i==j+2) || (i==j-2)) return -0.25;
			else return 0.0;
			return 0;	
	}

	void Evals (const mat &H, vec &d, mat &V)
	{
		eig_sym (H, d, V);
	}	

int main ()
{

	int i;
	int j;
	int I;
	int J;
	int n_start = 10;
	int n = n_start+1;
	int M = (2*(n_start))+1;
	long double diff;
	long double eps = 1e-8;
	long double swap_eigenvalue;
	long double new_eigenvalue;

	mat H2(2*n+1, 2*n+1);	
	vec d2(2*n+1);
	mat H1(M,M);
	mat V(M,M);
	vec d1(M);
	
	
	for (i = -n_start; i <= n_start; i++)
	for (j = -n_start; j <= n_start; j++)
	{
		J = n_start-j;
		I = n_start-i;
		H1(I,J) = kron(i,j);
	}

	for (n=n_start+1; n<10000; n++)
	{
		for (i = -n; i <= n; i++)
		for (j = -n; j <= n; j++)
		{
			J = n-j;
			I = n-i;	
			H2(I,J) = kron(i,j);
		}

		swap_eigenvalue = d1(0);
		new_eigenvalue = d2(0);
		diff = fabs(new_eigenvalue - swap_eigenvalue);
		//cout << "Difference: " << diff << endl;

		if (diff < eps) {
			//cout << "More accurate eigenvalue: " << new_eigenvalue << endl;
			break;
			
		}

		if (diff >= eps) swap_eigenvalue = new_eigenvalue;
	}

	

	Evals (H1, d1, V); 
	Evals  (H2, d2, V);
	
	//cout << "[ "<< d1(0) << ", " << d1(2) << ", " <<  d1(4) << ", " << d1(6) << ", " << d1(8) << " ]" << endl;
	cout << "[ "<< d2(0) << ", " << d2(2) << ", " <<  d2(4) << ", " << d2(6) << ", " << d2(8) << " ]" << endl;
	//cout << V << endl;
	//cout << H2 << endl;
	//cout << H1 << endl;
	return 0;
}

