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
#include <complex>

using namespace std;

int main()
{
	int N_SITES = 2;
	double e1 = 0.0;
	double e2 = 5.0;
	double t = 1.0;
	vector<vector<double> > ham_mat(N_SITES,vector<double>(N_SITES,0.0));
	ham_mat[0][0] = e1;
	ham_mat[1][1] = e2;
	ham_mat[0][1] = ham_mat[1][0] = -t;
	complex<double> z (1,2);

	cout << 1/z << endl;
	

	return 0;
}