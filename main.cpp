#include <iostream>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
//#include <eigen3/Eigen/Dense>
//#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <boost/math/special_functions/bernoulli.hpp>

using namespace std;
using Eigen::Matrix2cd;

#include "gauss_quadrature.h"
#include "lib.h"
#include "matrixInterpolation.h"

//miscellaneous constants
complex<double> i(0, 1);

//initializations
const int n_iterations = 40;
double threshold = 1.e-4;
int n_base_pts = 501;
int n_max = 10;
double xi = 1.0, xf = 1.1;
double u0 = 1.0, up0 = 1.0;


//define coefficient function in differential equation u''(x) + k^2(x) u(x) == 0 to be solved
double k(double x)
{
	return (x);	//simple toy example
}

double kp(double x)
{
	return (1.0);	//simple toy example
}

int main(int argc, char *argv[])
{
	double delta_x = xf - xi;
	for (int iter = 0; iter < n_iterations; ++iter)
	{
		//set up anything needed for calculation
		initialize_everything();

		int n_final = 0;

		//calculate time-ordered exponential
		vector<Matrix2cd> results = compute_Q(n_max, n_final);
		cout << "********************************************" << endl
				<< "result(" << n_final << ") = " << results[results.size()-1] << endl
				<< "********************************************" << endl;

		Eigen::Vector2cd ICs = get_ICs(u0, up0);

		vector<Eigen::Vector2cd> results_at_x;

		for (int ix = 0; ix < results.size(); ++ix)
		{
			double x = x_grid[ix];
			Eigen::Vector2cd result_at_x = results[ix] * ICs;
			results_at_x.push_back( result_at_x );
			complex<double> final_result_at_x = result_at_x(0) * exp(-i*x*k(x)) + result_at_x(1) * exp(i*x*k(x));
			cerr << "Output: " << x << "   " << final_result_at_x.real() << "   " << final_result_at_x.imag() << endl;
		}

		Eigen::Vector2cd result_at_xf = interpolateVectors_at_xf(xf, x_grid, results_at_x);
		Eigen::Vector2cd result_p_at_xf = interpolateMatrices_at_xf(xf, x_grid, Umatrix) * result_at_xf;

		double uf = ( result_at_xf(0) * exp(-i*xf*k(xf)) + result_at_xf(1) * exp(i*xf*k(xf)) ).real();	//assumes up is real
		double factor = k(xf) + xf*kp(xf);
		complex<double> factor1 = -i * result_at_xf(0) * factor + result_p_at_xf(0);
		complex<double> factor2 = i * result_at_xf(1) * factor + result_p_at_xf(1);
		double upf = ( factor1 * exp(-i*xf*k(xf)) + factor2 * exp(i*xf*k(xf)) ).real();	//assumes up0 is real

		cout << xi << "   " << u0 << "   " << up0 << ";   " << xf << "   " << uf << "   " << upf << endl;

		xi += delta_x;
		xf += delta_x;
		u0 = uf;
		up0 = upf;
	}

	return (0);
}
