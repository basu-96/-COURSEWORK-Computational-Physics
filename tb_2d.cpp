// ---------------------------------
// To compile:
// $ g++ -std=c++11 tb_2d.cpp -lgsl -lgslcblas -lm
// ---------------------------------
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <iterator>
#include <ctime>
#include <gsl/gsl_eigen.h>

using namespace std;

vector<int> nn_index(int );
double delta_func(int, int);
double lorenzian(double, double);
double ran_num();

int main()
{
	int N;
	N = 4;
	double temp = 0.0001
	double t = 1.0;
	double U = 6.0;
	double mu = U/2;
	int matrix_size = 2*N*N;
	vector<double> eigen_array(matrix_size,0.0);

	gsl_matrix * m1 = gsl_matrix_alloc(N*N,N*N);
	gsl_matrix * m = gsl_matrix_alloc(matrix_size, matrix_size);
	for(int i = 0; i < N*N; i++)
	{
		for (int j = 0; j < N*N; j++)
		{
			gsl_matrix_set(m1, i, j, 0.0);
			gsl_matrix_set(m1, i, j, 0.0);
			gsl_matrix_set(m, i, j, 0.0);
			gsl_matrix_set(m, i, j, 0.0);
			gsl_matrix_set(m, i+N*N, j+N*N, 0.0);
		}
	}
	
	cout << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			gsl_matrix_set(m1, i*N+j, i*N+(j+1)%N, -t);
			gsl_matrix_set(m1, i*N+j, ((i+1)%N)*N+j, -t);
			// cout << i*N+j << '\t' << i*N+(j+1)%N << '\t' <<  ((i+1)%N)*N+j << endl;
			gsl_matrix_set(m1, i*N+(j+1)%N, i*N+j, -t);
			gsl_matrix_set(m1, ((i+1)%N)*N+j, i*N+j, -t);
		}
	}
	cout << endl;

	for(int i = 0; i < N*N; i++)
	{
		for(int j = 0; j < N*N; j++)
		{
			gsl_matrix_set(m, i, j,  gsl_matrix_get(m1, i, j));
			gsl_matrix_set(m, i+N*N, j+N*N,  gsl_matrix_get(m1, i, j));
		}
	}
	
	// generating random values of n_nup and n_down at each site;
	vector<double> n_density(matrix_size, 0.0);
	for(int i = 0; i < matrix_size; i++)
	{
		double random_num = ran_num();
		n_density[i] = random_num;
		gsl_matrix_set(m,i,i, random_num);
	}
	ofstream outfile1;
	outfile1.open("hamiltonian_matrix_N_4.dat");

	for(int i = 0; i < matrix_size; i++)
	{
		for(int j = 0; j < matrix_size; j++)
		{
			outfile1 << gsl_matrix_get(m, i, j) << '\t';
		}
		outfile1 << endl;
	}
	outfile1.close();

	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(matrix_size);
	gsl_vector *eval = gsl_vector_alloc(matrix_size);
	gsl_matrix *evec = gsl_matrix_alloc(matrix_size,matrix_size);
	gsl_eigen_symmv(m, eval, evec, w);
	// gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
	ofstream outfile;
	outfile.open("mean_field.dat");
	for(int i = 0; i < matrix_size; i++)
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
	// cout << "min :" << min_eigen << endl;
	// cout << "max :" << max_eigen << endl;
	cout << "BW :" << max_eigen - min_eigen << endl;
	
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
double ran_num()
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);
	return dis(gen);
}