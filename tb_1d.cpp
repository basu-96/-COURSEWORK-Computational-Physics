#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <iterator>
#include <ctime>
// #include "functions.h"
#include <gsl/gsl_eigen.h>

using namespace std;

double delta_func(int, int);
double lorenzian(double, double);

int main()
{
	//-----number of sites
	int N = 50;
	//hamiltonian: H(m,n) = -t(delta(m,n+1) + delta(m,n-1));
	double t = 1.0;
	vector<double> eigen_array(N,0.0);
	// vector<vector<double> > hamiltonian_matrix(N, vector<double>(N, 0.0));
	gsl_matrix * m = gsl_matrix_alloc(N,N);

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			gsl_matrix_set(m, i, j, -t*delta_func(i,(j+1)%N)-t*delta_func((i+1)%N,j));
		}
	} 	

	// for(int i = 0; i < N; i++)
	// {
	// 	for(int j = 0; j < N; j++)
	// 	{
	// 		cout << gsl_matrix_get(m,i,j) << '\t';
	// 	}
	// 	cout << endl;
	// }
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(N);
	gsl_vector *eval = gsl_vector_alloc(N);
	gsl_matrix *evec = gsl_matrix_alloc(N,N);
	gsl_eigen_symmv(m, eval, evec, w);	
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
	ofstream outfile;
	outfile.open("tight_binding_1d.dat");
	for(int i = 0; i < N; i++)
	{
		double eval_i = gsl_vector_get(eval, i);
		eigen_array[i] = eval_i;
		outfile << eval_i << endl;
	}
	outfile.close();
	gsl_eigen_symmv_free (w);
	gsl_vector_free (eval);
	gsl_matrix_free (m);

	double max_eigen = *max_element(eigen_array.begin(), eigen_array.end());
	double min_eigen = *min_element(eigen_array.begin(), eigen_array.end());
	int n_points = 5000;// number of data points for DOS calculation
	double gamma = (max_eigen - min_eigen)/(double)N;
	cout << gamma << endl;
	gamma = 0.5;
	outfile.open("tb_1d_dos.dat");
	for(int i = 0; i <= n_points; i++)
	{
		// int i = 0;
		double x = min_eigen + (max_eigen - min_eigen)*(double)i/(double)n_points;
		double sum = 0;
		for(int j = 0; j < N; j++)
		{
			double eval_i = gsl_vector_get(eval, j);
			sum = sum + lorenzian(gamma, x - eval_i);
		}
		outfile << x << '\t' << sum << endl;
	}
	outfile.close();

	return 0;
}
double delta_func(int m, int n)
{
	if (m == n) return 1.0;
	else return 0.0;
}
double lorenzian(double gamma, double x)
{
	return (1/(2*M_PI))*gamma/(x*x + 0.25*gamma*gamma);
}