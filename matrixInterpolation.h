#ifndef MATRIX_INTERPOLATION_H
#define MATRIX_INTERPOLATION_H

#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/special_functions/bernoulli.hpp>

#include "gauss_quadrature.h"
#include "lib.h"

using namespace std;
using Eigen::Matrix2cd;

/*class matrixInterpolator
{
	private:
		//range from -1.0 to 1.0
		vector<double> base_points, base_weights;
	public:
		matrixInterpolator(int n_points_in) {n_points = n_points_in;}
		void Start();
		void Stop();
		void Reset();
		double printTime();
};*/

extern int n_base_pts;
extern vector<double> base_pts, base_wts;

//**********************************************************************
long binarySearch(vector<double> A, double value, bool skip_out_of_range = true, bool verbose = false);
vector<Matrix2cd> interpolateMatrices(vector<double> x_grid, vector<Matrix2cd> integrand, double a, double b, int kind);
Eigen::Vector2cd interpolateVectors_at_xf(double xf, vector<double> x_grid, vector<Eigen::Vector2cd> integrand);
Matrix2cd interpolateMatrices_at_xf(double xf, vector<double> x_grid, vector<Matrix2cd> integrand);
Matrix2cd interpCubicNonDirect(vector<double> x, vector<Matrix2cd> y, double xi, bool returnflag = false, double default_return_value = 0);

#endif
