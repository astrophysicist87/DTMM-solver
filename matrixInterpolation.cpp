#include <iostream>
#include <iomanip>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/special_functions/bernoulli.hpp>

#include "matrixInterpolation.h"

using namespace std;
using Eigen::Matrix2cd;

//**********************************************************************
long binarySearch(vector<double> A, double value, bool skip_out_of_range /*== true*/, bool verbose /*== false*/)
// Return the index of the largest number less than value in the list A
// using binary search. Index starts with 0.
// If skip_out_of_range is set to true, then it will return -1 for those
// samples that are out of the table range (default is true).
{
	int idx_i, idx_f, idx;
	int length = A.size();
	   
	idx_i = 0;
	idx_f = length-1;

	if(value > A[idx_f])
	{
		if (verbose) cerr << "binarySearch: desired value is too large, exceeding the end of the table: value = " << value << " and A[idx_f] = " << A[idx_f] << endl;
		if (skip_out_of_range) return -1;
		exit(1);
	}
   if(value < A[idx_i])
	{
		if (verbose) cerr << "binarySearch: desired value is too small, exceeding the beginning of table: value = " << value << " and A[idx_i] = " << A[idx_i] << endl;
		if (skip_out_of_range) return -1;
		exit(1);
	}
	idx = (int) floor((idx_f+idx_i)/2.);
	while((idx_f-idx_i) > 1)
	{
		if(A[idx] < value)
			idx_i = idx;
		else
			idx_f = idx;
		idx = (int) floor((idx_f+idx_i)/2.);
	}
	return(idx_i);
}

//**********************************************************************
vector<Matrix2cd> interpolateMatrices(vector<double> x_grid, vector<Matrix2cd> integrand, double a, double b, int kind)
{
	vector<Matrix2cd> interpolatedIntegrand;
	double cen = 0.5 * (a + b), hw = 0.5 * (b - a);

	if (kind == 1)	//cubic
	{
		for (int ib = 0; ib < n_base_pts; ++ib)
		{
			double x_point = cen + hw * base_pts[ib];
			interpolatedIntegrand.push_back( interpCubicNonDirect(x_grid, integrand, x_point, true) );	//true means use linear extrapolation at the end
		}
	}
	else			//linear is default
	{
		for (int ib = 0; ib < n_base_pts; ++ib)
		{
			double x_point = cen + hw * base_pts[ib];
			long idx = binarySearch(x_grid, x_point, true);	//true means abort when point is outside of range (shouldn't happen)
			interpolatedIntegrand.push_back( integrand[idx] + (integrand[idx+1]-integrand[idx])/(x_grid[idx+1]-x_grid[idx])*(x_point-x_grid[idx]) );
		}
	}

	return ( interpolatedIntegrand );
}

//**********************************************************************
Eigen::Vector2cd interpolateVectors_at_xf(double xf, vector<double> x_grid, vector<Eigen::Vector2cd> integrand)
{
	Eigen::Vector2cd interpolatedIntegrand;
	int len = x_grid.size();
	interpolatedIntegrand = integrand[len-2] + (integrand[len-1]-integrand[len-2])/(x_grid[len-1]-x_grid[len-2])*(xf-x_grid[len-2]);

	return ( interpolatedIntegrand );
}

//**********************************************************************
Matrix2cd interpolateMatrices_at_xf(double xf, vector<double> x_grid, vector<Matrix2cd> integrand)
{
	Matrix2cd interpolatedIntegrand;
	int len = x_grid.size();
	interpolatedIntegrand = integrand[len-2] + (integrand[len-1]-integrand[len-2])/(x_grid[len-1]-x_grid[len-2])*(xf-x_grid[len-2]);

	return ( interpolatedIntegrand );
}

//**********************************************************************
Matrix2cd interpCubicNonDirect(vector<double> x, vector<Matrix2cd> y, double xi, bool returnflag /*= false*/, double default_return_value /* = 0*/)
// Returns the interpreted value of y=y(x) at x=x0 using cubic polynomial interpolation method.
// -- x,y: the independent and dependent double tables; x is *NOT* assumed to be equal spaced but it has to be increasing
// -- xi: where the interpolation should be performed
{
	int size = x.size();
	if (size==1) {cout<<"interpCubicNondirect warning: table size = 1"<<endl; return y[0];}

	// if close to left end:
	if (abs(xi-x[0])<(x[1]-x[0])*1e-30) return y[0];

	// find x's integer index
	long idx = binarySearch(x, xi, true);

	if (idx < 0 || idx >= size-1)
	{
		if (!returnflag)	//i.e., if returnflag is false, exit
		{
			cerr << "interpCubicNonDirect(): index out of range!  Aborting!" << endl
				<< "interpCubicNonDirect(): size = " << size << ", x0 = " << xi << ", " << "idx=" << idx << endl;
			exit(1);
		}
		else
		{
			cerr << "interpCubicNonDirect(): Warning - index out of range!" << endl
				<< "interpCubicNonDirect(): size = " << size << ", x0 = " << xi << ", " << "idx=" << idx << endl
				<< "interpCubicNonDirect(): proceeding with linear extrapolation!" << endl;
			idx = (idx<0) ? 0 : size-2;	//uses linear extrapolation
		}
		//else return default_return_value;
	}

	if (idx==0)
	{
		// use linear interpolation at the left end
		return y[0] + (y[1]-y[0])/(x[1]-x[0])*(xi-x[0]);
	}
	else if (idx==size-2)
	{
		// use linear interpolation at the right end
		return y[size-2] + (y[size-1]-y[size-2] )/(x[size-1]-x[size-2] )*(xi-x[size-2]);
	}
	else
	{
		// use cubic interpolation
		Matrix2cd y0 = y[idx-1], y1 = y[idx], y2 = y[idx+1], y3 = y[idx+2];
		Matrix2cd y01=y0-y1, y02=y0-y2, y03=y0-y3, y12=y1-y2, y13=y1-y3, y23=y2-y3;
		long double x0 = x[idx-1], x1 = x[idx], x2 = x[idx+1], x3 = x[idx+2];
		long double x01=x0-x1, x02=x0-x2, x03=x0-x3, x12=x1-x2, x13=x1-x3, x23=x2-x3;
		long double x0s=x0*x0, x1s=x1*x1, x2s=x2*x2, x3s=x3*x3;
		long double denominator = x01*x02*x12*x03*x13*x23;
		Matrix2cd C0 = (x0*x02*x2*x03*x23*x3*y1
						+ x1*x1s*(x0*x03*x3*y2+x2s*(-x3*y0+x0*y3)+x2*(x3s*y0-x0s*y3))
						+ x1*(x0s*x03*x3s*y2+x2*x2s*(-x3s*y0+x0s*y3)+x2s*(x3*x3s*y0-x0*x0s*y3))
						+ x1s*(x0*x3*(-x0s+x3s)*y2+x2*x2s*(x3*y0-x0*y3)+x2*(-x3*x3s*y0+x0*x0s*y3))
						)/denominator;
		Matrix2cd C1 = (x0s*x03*x3s*y12
						+ x2*x2s*(x3s*y01+x0s*y13)
						+ x1s*(x3*x3s*y02+x0*x0s*y23-x2*x2s*y03)
						+ x2s*(-x3*x3s*y01-x0*x0s*y13)
						+ x1*x1s*(-x3s*y02+x2s*y03-x0s*y23)
						)/denominator;
		Matrix2cd C2 = (-x0*x3*(x0s-x3s)*y12
						+ x2*(x3*x3s*y01+x0*x0s*y13)
						+ x1*x1s*(x3*y02+x0*y23-x2*y03)
						+ x2*x2s*(-x3*y01-x0*y13)
						+ x1*(-x3*x3s*y02+x2*x2s*y03-x0*x0s*y23)
						)/denominator;
		Matrix2cd C3 = (x0*x03*x3*y12
						+ x2s*(x3*y01+x0*y13)
						+ x1*(x3s*y02+x0s*y23-x2s*y03)
						+ x2*(-x3s*y01-x0s*y13)
						+ x1s*(-x3*y02+x2*y03-x0*y23)
						)/denominator;
		return C0 + C1*xi + C2*xi*xi + C3*xi*xi*xi;
	}
}

//End of file
