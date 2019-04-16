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
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace std;

int main()
{
	int N_SITES = 2;
	double e1 = 0.0;
	double e2 = 2.0;
	// double t = 1.0;
	double omega = -5.0;
	int num_z_points = 1000;
	double eta = 0.01;
	double pi = acos(-1);
	gsl_complex z, g11, g22, G11, G22, t;
	GSL_SET_COMPLEX(&t, 1, 0);
	ofstream outfile;
	outfile.open("greens_func_2_site.dat");
	for(int i = 0; i <= num_z_points; i++)
	{
		GSL_SET_COMPLEX(&z, omega, eta);
		g11 = gsl_complex_inverse(gsl_complex_sub_real(z,e1));
		g22 = gsl_complex_inverse(gsl_complex_sub_real(z,e2));

		G11 = gsl_complex_inverse(gsl_complex_sub(gsl_complex_inverse(g11), gsl_complex_mul(gsl_complex_mul(gsl_complex_conjugate(t), g22), t)));
		G22 = gsl_complex_inverse(gsl_complex_sub(gsl_complex_inverse(g22), gsl_complex_mul(gsl_complex_mul(gsl_complex_conjugate(t), g11), t)));

		outfile << omega << '\t' << -GSL_IMAG(G11)/pi << '\t' << -GSL_IMAG(G22)/pi << endl;


		// cout << GSL_REAL(z) << '\t' << GSL_IMAG(z) << endl;
		omega = omega + 10.0/(double)num_z_points;
	}

	outfile.close();
	

	return 0;
}