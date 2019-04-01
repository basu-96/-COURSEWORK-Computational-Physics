#ifndef FUNCTIONS
#define FUNCTIONS

#include<vector>

using namespace std;

vector<double> cross(vector<double> &array_1, vector<double> &array_2)
{
	vector<double> result(3);
	result[0] = array_1[1]*array_2[2] - array_1[2]*array_2[1];
	result[1] = array_1[2]*array_2[0] - array_1[0]*array_2[2];
	result[2] = array_1[0]*array_2[1] - array_1[1]*array_2[0]; 
	return result;
}

double dot(vector<double> &array_1, vector<double> &array_2)
{
	double result = 0.0;
	result = array_1[0]*array_2[0] + array_1[1]*array_2[1] +
				array_1[2]*array_2[2];

	return result;
}

double random_num(double x, double y)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(x, y);
	return dis(gen);
}

void display_vector(vector<double> &array)
{
	int array_size = array.size();
	for(int i = 0; i < array_size; i++)
	{
		cout << array[i] << '\t'; 
	}
	cout << endl;
}

double norm(vector<double> &array)
{
	double norm = 0.0;
	int array_size = array.size();
	for(int i = 0; i < array_size; i++)
	{
		norm = norm + pow(array[i], 2);
	}
	norm = sqrt(norm);
}

vector<double> sum(vector<double> &array_1, vector<double> &array_2)
{
	int array_1_size = array_1.size();
	int array_2_size = array_2.size();
	if(array_1_size != array_2_size)
	{
		cout << "Warning! : Vectors being added are not of the same size." << endl;
		terminate();
	}
	vector<double> sum_array(array_1_size,0.0);

	for(int i = 0; i < array_1_size; i++)
	{
		sum_array[i] = array_1[i] + array_2[i];
	}
	return sum_array;
}

int check_inside_zone(vector<double> &current_q_vector, double dist_diff_tolerance, int number_BZ_surfaces, vector<vector<double> > &normal_unit_vec_BZ_surfaces, vector<double> &BZ_surface_normal_distance)
{
	int inside_zone = 1;
	for (int i = 0; i < number_BZ_surfaces; i++)
	{
		double current_q_vector_projection = abs(dot(current_q_vector, normal_unit_vec_BZ_surfaces[i]));
		double dist_diff = current_q_vector_projection - BZ_surface_normal_distance[i];
		if (dist_diff > dist_diff_tolerance)
		{
			inside_zone = 0;
			break;
		}
	}
	return inside_zone;
}

#endif