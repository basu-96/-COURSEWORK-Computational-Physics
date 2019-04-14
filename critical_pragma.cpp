#include<iostream>
#include<fstream>
#include<iomanip>
#include<math.h>
#include<cstdlib>
#include<sstream>
#include<vector>
#include <algorithm>
#include <iterator>
#include <omp.h>

using namespace std;


int square(int, int);

int main()
{
	// system("rm critical_pragma_test.dat");
	fstream outfile;
	outfile.open("critical_pragma_test.dat", ios::app);
	omp_set_dynamic(0);
	#pragma omp parallel for schedule(guided) num_threads(3)
	for (int i = 0; i < 100; i++)
	{
		#pragma omp critical
		outfile << i << endl;
	}
	outfile.close();
	return 0;
}

int square(int x, int y)
{
	return pow(x,2) + pow(y,2);
}