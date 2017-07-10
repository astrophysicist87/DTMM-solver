#ifndef LIB_H
#define LIB_H

#include <iostream>
#include <complex>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
//#include <eigen3/Eigen/Dense>
//#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <boost/math/special_functions/bernoulli.hpp>

using namespace std;
using Eigen::Matrix2cd;

/*USAGE: debugger(__LINE__, __FILE__);*/
void inline debugger(int cln, const char* cfn)
{
	cerr << "You made it to " << cfn << ":" << cln << "!" << endl;
	return;
}

//miscellaneous constants
extern complex<double> i;

//initializations
extern double threshold;
extern int n_base_pts;
extern int n_max;
extern double xi;
extern double xf;

extern vector<double> bernoulli_numbers;
extern vector<double> x_grid;
extern vector<double> base_pts;
extern vector<double> base_wts;

extern vector<Matrix2cd> Umatrix;					//n_base_pts (x) 2x2 matrices
extern vector<vector<Matrix2cd> > Omega_vector;		//n_max (x) n_base_pts (x) 2x2 matrices
extern vector<Matrix2cd> Omega_total;				//n_base_pts (x) 2x2 matrices
extern vector<vector<vector<Matrix2cd> > > Smatrix;	//n_max (x) n_max-1 (x) n_base_pts (x) 2x2 matrices

//function to define the problem
double k(double x);		//function value
double kp(double x);	//function derivative

long factorial(long n);
Matrix2cd comm(Matrix2cd A, Matrix2cd B);
Matrix2cd commN(Matrix2cd A, Matrix2cd B, const int N);

void initialize_everything();
void set_Umatrix();
Eigen::Vector2cd get_ICs(complex<double> u_xi, complex<double> up_xi);
void compute_Omega_at_x(int ix, int n_max);
Matrix2cd compute_Omega_term_at_x(int ix, int n);
Matrix2cd compute_S_j_n_at_x(int ix, int n, int j);
vector<Matrix2cd> compute_S_j_n(int n, int j);

vector<Matrix2cd> compute_Q(int n_max, int & n_final);

//End of file
#endif
