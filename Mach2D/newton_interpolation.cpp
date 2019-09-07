#include "stdafx.h"

double* get_newton_interpolation_coefficients(int n, double* xn, double* fn, double* c)
{
	// Calculates the coefficients of the newtonian interpolation
	// n  -  Degree of the interpolatory polynomial
	// xn - x nodes
	// fn - f(x)
	// c  - Coefficients of the interpolation

	// Auxiliary variables
	int i, j;
	
	double** M;

	// Initializating
	M = new double*[n + 1];
	for (i = 0; i < n + 1; i++)
	{
		M[i] = new double[n + 1];
	}

	for (i = 0; i <= n; i += 1)
	{
		M[i][0] = fn[i];;
	}

	for (j = 1; j <= n; j += 1)					// already fixed loop - verify
	{
		for (i = j; i <= n; i += 1)				// already fixed loop - verify
		{
			M[i][j] = (M[i][j - 1] - M[i - 1][j - 1]) / (xn[i] - xn[i - j]);
		}
	}
	
	for (j = 0; j <= n; j += 1)					// already fixed loop - verify
	{
		c[j] = M[j][j];
	}

	for (i = 0; i < n + 1; i++)
	{
		delete[] M[i];
	}
	delete[] M;
	M = nullptr;

	return c;
}

double get_newton_interpolation_f0(int n, double x, double* xn, double* c)
{
	// Calculates the value of the polynomial at a specified point x.
	// n  -  Degree of the interpolatory polynomial
	// x  -  Value of x where the interpolation will be carried out
	// xn - x nodes
	// c  -  Coefficients of the interpolation
	double get_newton_interpolation_f0;

	// Auxiliary variables
	int i, j;

	double aux1, aux2;

	aux1 = 0.0;

	for (i = 0; i <= n; i += 1)					// already fixed loop - verify
	{
		aux2 = 1.0;
		for (j = 0; j <= i - 1; j += 1)			// already fixed loop - verify
		{
			aux2 = aux2 * (x - xn[j]);
		}
		aux1 = aux1 + c[i] * aux2;
	}
	get_newton_interpolation_f0 = aux1;
	return get_newton_interpolation_f0;
}

double get_newton_interpolation_f1(int n, double x, double* xn, double* c)
{
	// Calculates the value of the polynomial first derivative at a specified point x.
	// n  -  Degree of the interpolatory polynomial
	// x  -  Value of x where the interpolation will be carried out
	// xn - x nodes
	// c  -  Coefficients of the interpolation
	double get_newton_interpolation_f1;

	// Auxiliary variables
	int i, j, k;

	double aux1, aux2, aux3;

	aux1 = 0.0;

	for (i = 0; i <= n; i += 1)					// already fixed loop - verify
	{
		aux2 = 0.0;
		for (j = 0; j <= i - 1; j += 1)			// already fixed loop - verify
		{
			aux3 = 1.0;
			for (k = 0; k <= i - 1; k += 1)		// already fixed loop - verify
			{
				if (k != j)
				{
					aux3 = aux3 * (x - xn[k]);
				}
			}
			aux2 = aux2 + aux3;
		}
		aux1 = aux1 + c[i] * aux2;
	}
	get_newton_interpolation_f1 = aux1;
	return get_newton_interpolation_f1;
}

interpolation get_ni_f0_vec_d1(int n, double* xn, double* fn, double xc)
{
	// Searches in the vector xn(1:n) the value of x (xc) where the interpolation
	// will be carried out. If xc is found, then a interpolation of f(xc) with a FIRST degree
	// polynomial is performed based on xn and fn. Depending on the results, the variable res
	// returns: res = 1, meaning "no problems found", res = -1, meaning "out of range", and
	// res = -2, meaning "error".
	// n  -  Number of nodes in the vector xn
	// xn - Vector of the nodes
	// fn - Vector of f(xn)
	// xc - Value of x where the interpolation must be done
	double fc;	// Interpolation of f(xc)
	int res;	// Status of the result

	// Auxiliary variables
	int kc;		// Index for xc
	double* xl;	// Local values of xn
	double* fl;	// Local values of fn
	double* c;	// Coefficients of the newtonian interpolation

	// Initializating
	xl = new double[1 + 1];
	fl = new double[1 + 1];
	c = new double[1 + 1];

	res = 1;

	fc = 0.0;

	kc = searcher(n, xn, xc);

	if (kc == -1)
	{
		// Out of range
		res = -1;
		interpolation result;
		result.fc = fc;
		result.status = res;
		return result;
	}
	else if (kc == -2)
	{
		// It was found an error
		res = -2;
		interpolation result;
		result.fc = fc;
		result.status = res;
		return result;
	}

	memcpy(xl, xn + kc, (kc+1) * sizeof(double));							// VERIFY
	memcpy(fl, fn + kc, (kc + 1) * sizeof(double));

	c = get_newton_interpolation_coefficients(1, xl, fl, c);

	fc = get_newton_interpolation_f0(1, xc, xl, c);

	interpolation result;
	result.fc = fc;
	result.status = res;

	delete [] xl;
	xl = nullptr;
	delete[] fl;
	fl = nullptr;
	delete[] c;
	c = nullptr;

	return result;
}

interpolation get_ni_f0_vec_d2(int n, double* xn, double* fn, double xc)
{
	// Searches in the vector xn(1:n) the value of x (xc) where the interpolation
	// will be carried out. If xc is found, then a interpolation of f(xc) with a SECOND degree
	// polynomial is performed based on xn and fn. Depending on the results, the variable res
	// returns: res = 1, meaning "no problems found", res = -1, meaning "out of range", and
	// res = -2, meaning "error".
	// n  -  Number of nodes in the vector xn
	// xn - Vector of the nodes
	// fn - Vector of f(xn)
	// xc - Value of x where the interpolation must be done
	double fc;	// Interpolation of f(xc)
	int res;	// Status of the result

	// Auxiliary variables
	int kc;		 // Index for xc

	double* xl; // Local values of xn
	double* fl; // Local values of fn
	double* c;	// Coefficients of the newtonian interpolation

	// Initializating
	xl = new double[2 + 1];
	fl = new double[2 + 1];
	c = new double[2 + 1];

	res = 1;

	fc = 0.0;

	kc = searcher(n, xn, xc);

	if (kc == -1)
	{
		// Out of range
		res = -1;
		interpolation result;
		result.fc = fc;
		result.status = res;
		return result;
	}
	else if (kc == -2)
	{
		// It was found an error
		res = -2;
		interpolation result;
		result.fc = fc;
		result.status = res;
		return result;
	}

	if (kc == 0)
	{
		memcpy(xl, xn + kc, (kc + 2) * sizeof(double));
		memcpy(fl, fn + kc, (kc + 2) * sizeof(double));

		c = get_newton_interpolation_coefficients(2, xl, fl, c);

		fc = get_newton_interpolation_f0(2, xc, xl, c);
	}
	else if (kc == n - 1 - 1)
	{
		memcpy(xl, xn + kc - 1, (kc + 1) * sizeof(double));
		memcpy(fl, fn + kc - 1, (kc + 1) * sizeof(double));

		c = get_newton_interpolation_coefficients(2, xl, fl, c);

		fc = get_newton_interpolation_f0(2, xc, xl, c);
	}
	else
	{
		if (std::fabs(xc - xn[kc]) < std::fabs(xc - xn[kc + 1]))
		{
			memcpy(xl, xn + kc - 1, (kc + 1) * sizeof(double));
			memcpy(fl, fn + kc - 1, (kc + 1) * sizeof(double));

			c = get_newton_interpolation_coefficients(2, xl, fl, c);

			fc = get_newton_interpolation_f0(2, xc, xl, c);
		}
		else
		{
			memcpy(xl, xn + kc, (kc + 2) * sizeof(double));
			memcpy(fl, fn + kc, (kc + 2) * sizeof(double));

			c = get_newton_interpolation_coefficients(2, xl, fl, c);

			fc = get_newton_interpolation_f0(2, xc, xl, c);
		}
	}
	interpolation result;
	result.fc = fc;
	result.status = res;

	delete[] xl;
	xl = nullptr;
	delete[] fl;
	fl = nullptr;
	delete[] c;
	c = nullptr;

	return result;
}

int searcher(int n, double* x, double xc)
{
	// Searches for the xc in the vector x(1:n). If xc was found, then
	// kc is the index in the vector x such that x(kc) <= xc <= x(kc+1), else
	// kc = -1 if xc is out of [x(1),x(n)] and
	// kc = -2 if an error was found.
	// n  -  Number of nodes
	// x  -  Vector of independent variable
	// xc - Value of x to be found
	int kc;	// Index in the vector such that x(kc) <= xc <= x(kc+1)

	// Auxiliary variables
	int i, ki, km, kf;

	ki = 1 - 1;
	kf = n - 1;

	if ((x[ki] - xc) * (x[kf] - xc) > 0.0)
	{
		// xc do not belongs to the interval [x[ki], x[kf]]
		kc = -1;
		return kc;
	}
	else
	{
		for (i = 1 - 1; i <= n - 1; i += 1)		// already fixed loop - verify
		{
			km = (ki + kf) / 2;
			if ((x[ki] - xc) * (x[km] - xc) > 0.0)
			{
				ki = km;
			}
			else
			{
				kf = km;
			}
			if (kf - ki < 2)
			{
				break;
			}
		}

		kc = ki;

		if (x[kc] > xc || x[kc + 1] < xc)
		{
			// It was found an error
			kc = -2;
		}
	}
	return kc;
}