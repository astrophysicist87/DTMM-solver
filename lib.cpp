#include <iostream>
#include <iomanip>

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

vector<double> bernoulli_numbers;
vector<double> x_grid;
vector<double> base_pts;
vector<double> base_wts;

vector<Matrix2cd> Umatrix;					//n_base_pts (x) 2x2 matrices
vector<vector<Matrix2cd> > Omega_vector;		//n_max (x) n_base_pts (x) 2x2 matrices
vector<Matrix2cd> Omega_total;				//n_base_pts (x) 2x2 matrices
vector<vector<vector<Matrix2cd> > > Smatrix;	//n_max (x) n_max-1 (x) n_base_pts (x) 2x2 matrices

long factorial(long n)
{
	return ( ( n==0 ) ? 1 : n*factorial(n-1) );
}

Matrix2cd comm(Matrix2cd A, Matrix2cd B)
{
	return ( A*B - B*A );

}

Matrix2cd commN(Matrix2cd A, Matrix2cd B, const int N)
{
	return ( (N==0) ? B : comm(A, commN(A, B, N-1)) );
}


Matrix2cd integrate(vector<Matrix2cd> integrand, double a, double b)
{
	double cen = 0.5 * (a + b), hw = 0.5 * (b - a);
	Matrix2cd result = Matrix2cd::Zero();

	//set matrix values interpolated to each point in integration
	//integrand ranges from xi to xf --> need to interpolate to get range from a to b
	int kind = 0;	// 0 - linear, 1 - cubic
	vector<Matrix2cd> integrand_at_x_points = interpolateMatrices(x_grid, integrand, a, b, kind);

	for (int ib = 0; ib < n_base_pts; ++ib)
	{
		double x_point = cen + hw * base_pts[ib];
		result += hw * base_wts[ib] * integrand_at_x_points[ib];
	}

	return (result);
}

void initialize_everything()
{
	//erase everything from previous iteration
	base_pts.clear();
	base_wts.clear();
	bernoulli_numbers.clear();
	x_grid.clear();
	Umatrix.clear();
	Omega_total.clear();
	Omega_vector.clear();
	Smatrix.clear();

	//initialize base integration points
	gauss_quadrature(n_base_pts, 1, 0.0, 0.0, -1.0, 1.0, &base_pts, &base_wts);

	//set needed bernoulli numbers
	for (int in = 0; in <= n_max; ++in)
		bernoulli_numbers.push_back( (in % 2 == 1) ? 0 : boost::math::bernoulli_b2n<double>(in / 2) );
	bernoulli_numbers[1] = -0.5;	//the only non-zero odd Bernoulli number

	//set up grid at which to calculate time-ordered exponential
	//for (int ix = 0; ix < n_grid_pts; ++ix)
	//	x_grid.push_back( xi + (xf - xi) / (double)(n_grid_pts - 1) );
	double cen = 0.5 * (xi + xf), hw = 0.5 * (xf - xi);
	for (int ix = 0; ix < (int)base_pts.size(); ++ix)
		x_grid.push_back( cen + hw * base_pts[ix] );

	//set U(x) needed for calculation
	set_Umatrix();

	//set the final Omega running sum (function of x, so a vector of length n_x_pts)
	Omega_total.resize(x_grid.size());
	Omega_vector.resize(n_max);
	Smatrix.resize(n_max);

	return;
}

void set_Umatrix()
{
	//Umatrix.reserve(n_base_pts);
	for (int ix = 0; ix < n_base_pts; ++ix)
	{
		Matrix2cd U_at_x(2,2);
		double x_local = x_grid[ix];
		double k_x_local = k(x_local);
		double kp_x_local = kp(x_local);

		U_at_x(0,0) = -1.0 + 2.0*i*x_local * k_x_local;
		U_at_x(0,1) = exp(2.0*i*x_local * k_x_local);
		U_at_x(1,0) = exp(-2.0*i*x_local * k_x_local);
		U_at_x(1,1) = -1.0 - 2.0*i*x_local * k_x_local;

		U_at_x *= (0.5 * kp_x_local / k_x_local);
		//cout << "********************************************" << endl
		//		<< setprecision(20) << setw(24) << x_local << "   " << U_at_x << endl
		//		<< "********************************************" << endl;
		Umatrix.push_back(U_at_x);
	}

	return;
}

//returns u1(xi) and u2(xi) corresponding to given solution ICs u(xi) and u'(xi)
Eigen::Vector2cd get_ICs(complex<double> u_xi, complex<double> up_xi)
{
	complex<double> k_xi = k(xi);
	complex<double> kp_xi = kp(xi);

	complex<double> e1 = exp(i*xi*k_xi);
	complex<double> e2 = e1*e1;
	complex<double> e3 = e1*e1*e1;
	complex<double> e4 = e1*e1*e1*e1;

	complex<double> prefactor = 0.5 * kp_xi / k_xi;
	complex<double> U00 = prefactor * ( -1.0 + 2.0*i*xi * k_xi );
	complex<double> U01 = prefactor * exp(2.0*i*xi * k_xi);
	complex<double> U10 = prefactor * exp(-2.0*i*xi * k_xi);
	complex<double> U11 = prefactor * ( -1.0 - 2.0*i*xi * k_xi );

	complex<double> denominator = U01 - e4*U10 - e2*( U00 - U11 - 2.0*i*k_xi - 2.0*i*xi*kp_xi );

	complex<double> u10 = e1*u_xi*U01 + e3*( u_xi*U11 - up_xi + i*u_xi*( k_xi + xi*kp_xi ) );
	complex<double> u20 = -e1 * ( u_xi * ( U00 + e2*U10 ) - up_xi - i * u_xi * ( k_xi + xi*kp_xi ) );

	Eigen::Vector2cd result;
	result(0) = u10 / denominator;
	result(1) = u20 / denominator;
	
	return ( result );
}

bool test_convergence(complex<double> previous, complex<double> current)
{
	return (abs( (current-previous) / previous ) <= threshold);
}

bool test_convergence(Matrix2cd previous, Matrix2cd current)
{
	return ( not ( test_convergence( previous(0,0), current(0,0) )
					and test_convergence( previous(0,1), current(0,1) )
					and test_convergence( previous(1,0), current(1,0) )
					and test_convergence( previous(1,1), current(1,1) )
			 ) );
}

vector<Matrix2cd> compute_Q(int n_max, int & n_final)
{
	int n = 1;
	bool NOT_CONVERGED = true;
	Matrix2cd Q_at_xf_current = Matrix2cd::Identity();	//initial result
	Matrix2cd Q_at_xf_previous = Matrix2cd::Identity();	//store previous result to check for convergence

	while (NOT_CONVERGED and n <= n_max)
	{
		cout << "Computing Q for n = " << n << endl;

		cout << "********************************************" << endl
				<< "result(" << n-1 << ") = " << Q_at_xf_current << endl
				<< "********************************************" << endl;


		//calculate full vector of Omega(x) at level n
		for (int ix = 0; ix < (int)x_grid.size(); ++ix)
		{
			compute_Omega_at_x(ix, n);
		}

		//test for convergence by looking at point closest to xf,
		//since this is the last one to converge
		Q_at_xf_previous = Q_at_xf_current;
		Q_at_xf_current = ( Omega_total[(int)x_grid.size() - 1] ).exp();	//computes exp(M) of matrix M at x==xf

		NOT_CONVERGED = test_convergence(Q_at_xf_previous, Q_at_xf_current);
		if (not NOT_CONVERGED) cout << "We've converged!" << endl;
		//cout << "n = " << n << ": " << Q_at_xf_previous << endl;
		//cout << "n = " << n << ": " << Q_at_xf_current << endl;
		n++;
	}

	//set Q for all x
	vector<Matrix2cd> Q_at_x;
	for (int ix = 0; ix < (int)x_grid.size(); ++ix)
		Q_at_x.push_back( ( Omega_total[ix] ).exp() );

	//record how many iterations were required for convergence
	n_final = n;
	
	//return (Q_at_xf_current);
	return (Q_at_x);
}

void compute_Omega_at_x(int ix, int n)
{
	Omega_vector[n-1].push_back( compute_Omega_term_at_x(ix, n) );	//note zero-indexing
	Omega_total[ix] += Omega_vector[n-1][ix];
	//if (n == n_max)
	//	cout << "********************************************" << endl
	//			<< "Omega_" << n << "(" << setprecision(20) << x_grid[ix] << ") = " << Omega_vector[n-1][ix] << endl
	//			<< "********************************************" << endl;

	return;
}

Matrix2cd compute_Omega_term_at_x(int ix, int n)
{
	double x = x_grid[ix];
	Matrix2cd result = Matrix2cd::Zero();
	if (n == 1)
	{
		result = integrate(Umatrix, xi, x);
	}
	else
	{
		vector<vector<Matrix2cd> > S_n;
		for (int j = 1; j <= n - 1; ++j)
		{
			double Bj = bernoulli_numbers[j];
			double jf = (double)factorial(j);
			vector<Matrix2cd> tmp;
			if (ix==0)
			{
				tmp = compute_S_j_n(n, j);
				S_n.push_back( tmp );
			}
			else
			{
				tmp = Smatrix[n-2][j-1];
			}
			result += (Bj/jf) * integrate(tmp, xi, x);			//note zero-indexing
		}
		if (ix==0)
			Smatrix[n-2] = S_n;				//Smatrix[0] <--> n==2
		//cout << "********************************************" << endl
		//		<< setprecision(20) << "len = " << Smatrix[0][0].size() << ": S^(1)_2(" << x << ") = " << Smatrix[0][0][ix] << endl
		//		<< "********************************************" << endl;
	}

	return ( result );
}

vector<Matrix2cd> compute_S_j_n(int n, int j)
{
	vector<Matrix2cd> local_Sjn;
	for (int ix = 0; ix < n_base_pts; ++ix)
	{
		local_Sjn.push_back( compute_S_j_n_at_x(ix, n, j) );
	}

	return ( local_Sjn );
}

Matrix2cd compute_S_j_n_at_x(int ix, int n, int j)
{
	Matrix2cd sum = Matrix2cd::Zero();
	if (j == 1)
	{
		sum = comm(Omega_vector[n-1 - 1][ix], Umatrix[ix]);							//note zero-indexing
	}
	else if (j == n - 1)
	{
		sum = commN(Omega_vector[1-1][ix], Umatrix[ix], n-1);						//note zero-indexing
	}
	else
	{
		for (int m = 1; m <= n - j; ++m)
		{
			int im = m-1;
			int ij = j-1;
			//cout << Smatrix.size() << endl;
			//cout << Smatrix[0].size() << endl;
			//cout << Smatrix[0][0].size() << endl;
			sum += comm( Omega_vector[im][ix], Smatrix[n-m-2][ij-1][ix] );		//note zero-indexing
		}
	}
	//if (ix==(n_base_pts - 1)/2)
	//	cout << "********************************************" << endl
	//			<< setprecision(20) << n << ", " << j << ": S^(1)_2(" << x_grid[ix] << ") = " << sum << endl
	//			<< "********************************************" << endl;


	return (sum);
}


//End of file
