#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <iterator>
#include <ctime>
#include "functions.h"
#include <gsl/gsl_eigen.h>

using namespace std;


int main()
{
	// Primitive lattice vectors
	vector<double> a1(3,0.0);
	vector<double> a2(3,0.0);
	vector<double> a3(3,0.0);
	a1[0] = 2.0;
	a1[1] = 0.0;
	a1[2] = 0.0;

	a2[0] = 0.0;
	a2[1] = 2.0;
	a2[2] = 0.0;

	a3[0] = 0.0;
	a3[1] = 0.0;
	a3[2] = 1.0;
	// Reciprocal lattice vectors
	vector<double> b1 = cross(a2,a3);
	vector<double> b2 = cross(a3,a1);
	// Unit vectors of b1 and b2
	vector<double> k1(3,0.0);
	vector<double> k2(3,0.0);

	vector<double> temp1 = cross(a2,a3);
	double reciprocal_factor = 2*M_PI/dot(a1,temp1);
	for(int i = 0; i < b1.size(); i++)
	{
		b1[i] = reciprocal_factor*b1[i];
		b2[i] = reciprocal_factor*b2[i];
	}
	for(int i = 0; i < b1.size(); i++)
	{
		k1[i] = b1[i]/norm(b1);
		k2[i] = b2[i]/norm(b2);
	}

	// Generating the k points
	int N1 = 100; // Number of lattice points in the direction of a1
	int N2 = 100; // Number of lattice points in the direction of a2
	vector<vector<double> > current_q_vectors (N1*N2, vector<double>(3, 0.0));
	int q_vectors_count = 0;

	for(int i = 0; i < N1; i++)
	{
		for(int j = 0; j < N2; j++)
		{
			q_vectors_count = q_vectors_count + 1;
			vector<double> temp (3,0.0);
			for(int t = 0; t < 3; t++)
			{
				temp[t] = (double) i / N1 * b1[t] + (double) j / N2 * b2[t];
			} 
			current_q_vectors[q_vectors_count] = temp;
		}
	}
	ofstream outfile;
	outfile.open("raw_q_vectors.dat");
	for (int i = 0; i < q_vectors_count; i++)
	{
		outfile << current_q_vectors[i][0] << '\t' << current_q_vectors[i][1] << endl;
	}
	outfile.close();
	//------------------Generating translation vectors -----------------------------
	// int translation_limit_q1 = 3;
	// int translation_limit_q2 = 3;

	// int num_translation_vectors = (2*translation_limit_q1 + 1) * (2*translation_limit_q1 + 1);

	// int translation_vectors_count = 0;

	// vector<vector<double> > translation_vectors (num_translation_vectors, vector<double> (3,0.0));

	// for (int i = -translation_limit_q1; i <= translation_limit_q1; i++)
	// {
	// 	for(int j = -translation_limit_q2; j <= translation_limit_q2; j++)
	// 	{
	// 		vector<double> temp1 (3,0.0);
	// 		vector<double> temp2 (3,0.0);
	// 		for(int t = 0; t < 3; t++)
	// 		{
	// 			temp1[t] = i*b1[t];
	// 			temp2[t] = j*b2[t];
	// 		}
	// 		translation_vectors[translation_vectors_count] = sum(temp1, temp2);
	// 		translation_vectors_count = translation_vectors_count + 1;
	// 	} 
	// }
	
	// //------------------Defining first BZ surfaces ------------------------------
	// int number_BZ_surfaces = 6;
	
	// vector<vector<double> > normal_unit_vec_BZ_surfaces (number_BZ_surfaces, vector<double> (3,0.0));
	// vector<double> temp (3,0.0);

	// normal_unit_vec_BZ_surfaces[0] = k1;

	// normal_unit_vec_BZ_surfaces[1] = k2;

	// normal_unit_vec_BZ_surfaces[2] = sum(k1,k2);
	
	// for(int i = 0; i < 3; i++)
	// {
	// 	temp[i] = -k1[i];
	// }
	// normal_unit_vec_BZ_surfaces[3] = temp;
	// for(int i = 0; i < 3; i++)
	// {
	// 	temp[i] = -k2[i];
	// }
	// normal_unit_vec_BZ_surfaces[4] = temp;
	// for(int i = 0; i < 3; i++)
	// {
	// 	temp[i] = -k1[i]-k2[i];
	// }
	// normal_unit_vec_BZ_surfaces[5] = temp;

	// //-----------------------Defining normal distances of surfaces from (0,0) 
	// vector<double> BZ_surface_normal_distance (number_BZ_surfaces,0.0);

	// BZ_surface_normal_distance[0] = norm(b1)/2;
	// BZ_surface_normal_distance[1] = norm(b2)/2;
	// temp = sum(b1, b2);
	// BZ_surface_normal_distance[2] = norm(temp)/2;
	// BZ_surface_normal_distance[3] = norm(b1)/2;
	// BZ_surface_normal_distance[4] = norm(b2)/2;
	// BZ_surface_normal_distance[5] = norm(temp)/2;


	// // ---------------------------- Translating the vectors ----------------------
	// vector<vector<double> > first_BZ_q_vectors (N1*N2, vector<double>(3, 0.0));

	// double dist_diff_tolerance = 1e-12;

	// for (int i = 0; i < q_vectors_count; i++)
	// {
	// 	vector<double> current_q_vector = current_q_vectors[i];
	// 	vector<double> current_q_vector_translated (3,0.0);

	// 	for(int j = 0; j < num_translation_vectors; j++)
	// 	{
	// 		current_q_vector_translated = sum(current_q_vector, translation_vectors[j]);

	// 		int inside_zone =  check_inside_zone(current_q_vector_translated, dist_diff_tolerance, number_BZ_surfaces, normal_unit_vec_BZ_surfaces, BZ_surface_normal_distance);
	// 		if (inside_zone == 1)
	// 		{
	// 			break;
	// 		}			
	// 	}
	// 	first_BZ_q_vectors[i] = current_q_vector_translated;
	// }
	
	// outfile.open("first_BZ_q_vectors.dat");

	// for(int i = 0; i < q_vectors_count; i++)
	// {
	// 	outfile << first_BZ_q_vectors[i][0] << '\t' << first_BZ_q_vectors[i][1] << endl;
	// }
	// outfile.close();


	//-----------------------GSL Test Routine -----------------
	outfile.open("eigen_values.dat");
	for (int t = 0; t < q_vectors_count; t++)
	{

		vector<double> q = current_q_vectors[t];

	  // double data[] = { 1.0  , 1/2.0, 1/3.0, 1/4.0,
   //                  1/2.0, 1/3.0, 1/4.0, 1/5.0,
   //                  1/3.0, 1/4.0, 1/5.0, 1/6.0,
   //                  1/4.0, 1/5.0, 1/6.0, 1/7.0 };
		double data[] = {0 , -2 * cos(q[0]) , -2 * cos(q[1]),
						-2 * cos(q[0]) , 0 ,  0,
						-2 * cos(q[1]) , 0 , 0};

  gsl_matrix_view m
    = gsl_matrix_view_array (data, 3, 3);

  gsl_vector *eval = gsl_vector_alloc (3);
  gsl_matrix *evec = gsl_matrix_alloc (3, 3);

  gsl_eigen_symmv_workspace * w =
    gsl_eigen_symmv_alloc (3);

  gsl_eigen_symmv (&m.matrix, eval, evec, w);

  gsl_eigen_symmv_free (w);

  gsl_eigen_symmv_sort (eval, evec,
                        GSL_EIGEN_SORT_ABS_ASC);

  outfile << q[0] << '\t' << q[1] << '\t';

  {
    int i;

    for (i = 0; i < 3; i++)
      {
        double eval_i
           = gsl_vector_get (eval, i);
        gsl_vector_view evec_i
           = gsl_matrix_column (evec, i);

        // printf ("eigenvalue = %g\n", eval_i);
           cout << "eigenvalue = " << eval_i << endl;
           outfile << eval_i << '\t';
        // printf ("eigenvector = \n");
        // gsl_vector_fprintf (stdout,
        //                     &evec_i.vector, "%g");
      }

      outfile << endl;
  }

  gsl_vector_free (eval);
  gsl_matrix_free (evec);
}
outfile.close();

	return 0;
}