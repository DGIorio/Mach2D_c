#include "stdafx.h"

// Contains modules and subroutines related to the solution of the linear systems.


double norm_l1_5d(int nx, int ny, double* var, double* b, double** a, double norm)
{
	// Auxiliary variables
	int i, j, np, nps, npn, npw, npe;
	double aux1;

	// Norm is calculated taking into account only real volumes
	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			nps = np - nx;
			npn = np + nx;
			npw = np - 1;
			npe = np + 1;

			norm = norm + std::fabs(
				  a[np][0] * var[nps]
				+ a[np][1] * var[npw]
				+ a[np][2] * var[np]
				+ a[np][3] * var[npe]
				+ a[np][4] * var[npn]
				- b[np]);
		}
	}
	return norm;
}

double norm_l1_9d(int nx, int ny, double* var, double* b, double** a, double norm)
{
	// Auxiliary variables
	int i, j, np, nps, npn, npw, npe, npsw, npse, npnw, npne;
	double aux1;

	// Norm is calculated taking into account only real volumes
	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			nps = np - nx;
			npn = np + nx;
			npw = np - 1;
			npe = np + 1;
			npsw = nps - 1;
			npse = nps + 1;
			npnw = npn - 1;
			npne = npn + 1;

			norm = norm + std::fabs(
				  a[np][0] * var[npsw]
				+ a[np][1] * var[nps]
				+ a[np][2] * var[npse]
				+ a[np][3] * var[npw]
				+ a[np][4] * var[np]
				+ a[np][5] * var[npe]
				+ a[np][6] * var[npnw]
				+ a[np][7] * var[npn]
				+ a[np][8] * var[npne]
				- b[np]);

		}
	}
	return norm;
}

double norm_l1_b(int nx, int ny, double* b)
{
	// Auxiliary variables
	int i, j, np;
	double norm;

	// Norm is calculated taking into account only real volumes
	norm = 0.0;

	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;

			norm = norm + std::fabs(b[np]);
		}
	}
	return norm;
}