#include<iostream>
#include<fstream>
#include<iomanip>
#include<math.h>
#include<cstdlib>
#include<sstream>
#include<vector>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <omp.h>

#include "gnuplot.h"

using namespace::std;
using std::ios;

// #define THREADS omp_get_max_threads()

double ran_num();
double energy(vector<vector<double> > &, int);
double magnetization(vector<vector<double> > &, int);
vector<vector<double> > config_gen(int); // updates the spin_grid.dat file

int main()
{
	// omp_set_dynamic(0);
	int grid_size;
	int system_sweeps;
	double T = 1.5;
	system_sweeps = 500;
	grid_size = 5;
	vector<vector<double> > spin_grid(grid_size,vector<double>(grid_size,0));
	spin_grid = config_gen(grid_size);

	
	
	return 0;
}

double ran_num()
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);
	return dis(gen);
}
double energy(vector<vector<double> > &spin_grid, int grid_size)
{
	double result = 0;
	for(int i = 0; i < grid_size; i++)
	{
		for(int j = 0; j < grid_size; j++)
		{
			result = result + (-spin_grid[i][j]) * (spin_grid[(i-1 + grid_size)%grid_size][j] + spin_grid[(i+1 + grid_size)%grid_size][j] + spin_grid[i][(j-1 + grid_size)%grid_size] + spin_grid[i][(j + 1 + grid_size)%grid_size]);
		}
	}
	result = result / 2;
	return result;
}

double magnetization(vector<vector<double> > &spin_grid, int grid_size)
{
	double result = 0;
	for( int i = 0; i < grid_size; i++)
	{
		for(int j = 0; j < grid_size; j++)
		{
			result = result + spin_grid[i][j];
		}
	}
	return abs(result);
}
vector<vector<double>> config_gen(int grid_size)
{
	vector<vector<double> > spin_grid(grid_size,vector<double>(grid_size,0));
	for(int i = 0; i < grid_size; i++)
	{
		for(int j = 0; j < grid_size; j++)
		{
			double random_number = ran_num();
			// cout << i << '\t' << j << '\t' << random_number << endl;
			if (random_number < 0.5)
			{
				spin_grid[i][j] = 1;
			}
			else
			{
				spin_grid[i][j] = -1;
			}
		}
	}
	return spin_grid;
}