// ---------------------------------
// To compile:
// $ g++ -std=c++11 single_spin_flip_MC.cpp
// ---------------------------------

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <ctime>
#include <omp.h>

using namespace std;

// #define THREADS omp_get_max_threads()

double ran_num();
double energy(vector<vector<double> > &, int);
double magnetization(vector<vector<double> > &, int);

int main()
{
	// omp_set_dynamic(0);
	int grid_size;
	int system_sweeps;
	double T = 1.5;
	system_sweeps = 500;
	grid_size = 50;
	vector<vector<double> > spin_grid(grid_size,vector<double>(grid_size,0));
	// generating the nitial random configuration
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
	// display the grid
	// for(int i = 0; i < grid_size; i++)
	// {
	// 	for(int j = 0; j < grid_size; j++)
	// 	{
	// 		cout << spin_grid[i][j] << '\t';
	// 	}
	// 	cout << '\n';
	// }
	// cout << endl;
	// cout << magnetization(spin_grid,grid_size) << endl;
	
	cout << energy(spin_grid, grid_size) << endl;
	// new config and old config
	cout << "Number of system-sweeps to complete :" << system_sweeps << endl;
	
	vector<vector<double> > old_config;
	vector<vector<double> > new_config;
	old_config = spin_grid;
	ofstream outfile1, outfile2, outfile3;
	outfile1.open("energy_vs_system_sweep.dat");
	outfile2.open("average_energy_vs_temperature.dat");
	outfile3.open("average_mag_vs_temperature.dat");

	// #pragma omp parallel for firstprivate(new_config, old_config) schedule(guided) num_threads(3)
	for(int temp = 100; temp > 0; temp = temp - 1)
	{
		// #pragma omp critical
		double T = (double)temp / (double)10;
		vector<double> tail_energy; 
		vector<double>  tail_magnetization;
		for(int sys_swp = 0; sys_swp < system_sweeps; sys_swp++)
	{
		for(int i = 0; i < grid_size; i++)
		{
			for(int j = 0; j < grid_size; j++)
			{
				new_config = old_config;
				new_config[i][j] = -new_config[i][j];
				double del_energy = energy(new_config,grid_size) - energy(old_config,grid_size);
				double r = ran_num();
				if (r < exp(-del_energy/T))
				{
					old_config = new_config;
				}
			}
		}	
		outfile1 << sys_swp << '\t' << energy(new_config, grid_size)/pow(grid_size,2) << endl;
		if(sys_swp > 250)
		{
			tail_energy.push_back(energy(new_config, grid_size)/pow(grid_size,2));
			tail_magnetization.push_back(magnetization(new_config, grid_size)/pow(grid_size,2));
		}
		// cout << sys_swp << endl;
	}		
		outfile2 << T << '\t' << (double) accumulate(tail_energy.begin(), tail_energy.end(), 0.0) / tail_energy.size() << endl;
		outfile3 << T << '\t' << (double) accumulate(tail_magnetization.begin(), tail_magnetization.end(), 0.0) / tail_magnetization.size() << endl;
		cout << T << '\t' << accumulate(tail_energy.begin(), tail_energy.end(), 0.0) / tail_energy.size(); 
		cout << '\t' << accumulate(tail_magnetization.begin(), tail_magnetization.end(), 0.0) / tail_magnetization.size() << endl;
	}

	outfile1.close();
	outfile2.close();
	outfile3.close();
	// cout << "Execuation time for " << system_sweeps << "sweeps : " << difftime(tend, tstart) << "seconds" << endl;
	// cout << "Size of energy tail : " << tail_energy.size() << endl;
	
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