#include "stdafx.h"

//Contains user - DEPENDENT subroutines and functions
//such as initial and boundary conditions as well as fluid thermophysical properties.

coordinates get_grid_boundary_power_law(int nx, int ny, double akn, double aks, double la, double lb, double lr, double rb, double* x, double* y)
{
	int i, j, np;
	double aux;
	double const pi = 3.141592653589793238463;	//requires less calculation time
	//double* x;
	//double* y;
	coordinates result;

	//x = new double[nx*ny];
	//y = new double[nx*ny];

	// North boundary
	j = ny - 1;
	for (i = 0; i <= nx - 1 - 1; i += 1) 		// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		x[np] = la * std::pow((double)(i) / (double)(nx - 2), akn) - (la - lr);		// already fixed i
		aux = (x[np] - lr) / la;
		y[np] = lb * sqrt(1.0 - aux*aux);
	}

	// South boundary
	j = 1;
	for (i = 0; i <= nx - 1 - 1; i += 1)		// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		x[np] = lr * std::pow((double)(i) / (double)(nx - 2), aks);					// already fixed i
		y[np] = f_body_shape(x[np], lr, rb);
	}

	result.x = x;
	result.y = y;
	
	return result;
}

double f_body_shape(double x, double lr, double rb)
{
	// Defines the function g(x) of the body shape.
	//return pow(rb*(x / lr), 1.0);
	return rb * (x / lr);
}

double* set_cp(int nx, int ny, double* cp)
{
	// Calculates cp at the center of each real volume.
	//double* cp;
	//cp = new double[nx*ny];
	std::fill_n(cp, nx*ny, 1003.8);	// or for (int i = 0; i < nx*ny; i++)
	return cp;
}

double* set_gamma(int nx, int ny, double Rg, double* cp, double* gcp)
{
	// Calculates gamma = Cp / Cv at the center of each real volume based on Cp and Rg.
	int i;
	//double* gcp;
	//gcp = new double[nx*ny];
	//std::transform(gcp.begin(), gcp.end(), gcp.begin(), bind2nd(std::minus<float>(), Rg));
	for (i = 0; i <= nx*ny - 1; i += 1)
	{
		gcp[i] = cp[i] / (cp[i] - Rg);
	}
	return gcp;
}

coeffs_source set_bcu(int nx, int ny, double UF, double* xk, double* yk, double* alphae, double* betae, double* u, double* v, double* Uce, double* Vcw, double** a, double* b)
{
	// Calculates the coefficients and source of the linear system of u for the fictitious volumes based on the boundary conditions.
	//double** a(nx*ny, std::vector<double>(9));
	//double* b(nx*ny);

	// Auxiliary variables
	int i;
	int j;
	int k;
	int np;
	int npw;
	int npe;
	int nps;
	int npn;
	int npse;
	int npsw;
	int npne;
	int npnw;

	int npnnw;
	int npssw;
	int npnne;
	int npsse;
	int npss;
	int npnn;

	double aux;
	double u_eta;

	// West boundary
	i = 0;			// already fixed loop - verify
	j = 2;

	np = nx * (j - 1) + i;
	npn = np + nx;
	npe = np + 1;
	npne = npn + 1;
	npnn = npn + nx;
	npnne = npnn + 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}
	a[np][4] =  1.0;
	a[np][5] = -1.0;

	u_eta = 0.5*(-u[npnne] + 4.0*u[npne] - 3.0*u[npe]);

	b[np] = -betae[np] / alphae[np] * u_eta;

	for (j = 3; j <= ny-2; j += 1)		// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npn = np + nx;
		npse = nps + 1;
		npne = npn + 1;

		for (k = 0; k <= 9 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][4] = 1.0;
		a[np][5] = -1.0;

		u_eta = 0.5 * (u[npne] - u[npse]);

		b[np] = -betae[np] / alphae[np] * u_eta;
	}

	j = ny - 1;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npe = np + 1;
	npse = nps + 1;
	npss = nps - nx;
	npsse = npss + 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}
	a[np][4] = 1.0;
	a[np][5] = -1.0;

	u_eta = 0.5 * (u[npsse] - 4.0 * u[npse] + 3.0 * u[npe]);

	b[np] = -betae[np] / alphae[np] * u_eta;


	// East boundary
	i = nx - 1;			// already fixed loop - verify

	j = 2;

	np = nx * (j - 1) + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;
	npnnw = npnw + nx;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}
	a[np][3] = -1.0;
	a[np][4] = 1.0;

	u_eta = 0.5 * (-u[npnnw] + 4.0 * u[npnw] - 3.0 * u[npw]);

	b[np] = -Vcw[j] / Uce[npw] * u_eta;


	for (j = 3; j <= ny - 2; j += 1)		// already fixed loop - verify
	{

		np = nx * (j - 1) + i;
		nps = np - nx;
		npn = np + nx;
		npw = np - 1;
		npsw = nps - 1;
		npnw = npn - 1;

		for (k = 0; k <= 9 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][3] = -1.0;
		a[np][4] = 1.0;

		u_eta = 0.50 * (u[npnw] - u[npsw]);

		b[np] = -Vcw[j] / Uce[npw] * u_eta;
	}

	j = ny - 1;

	np = nx * (j - 1) + i;
	nps = np - nx;
	npw = np - 1;
	npsw = nps - 1;
	npssw = npsw - nx;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}
	a[np][3] = -1.0;
	a[np][4] = 1.0;

	u_eta = 0.5 * (3.0 * u[npw] - 4.0 * u[npsw] + u[npssw]);

	b[np] = -Vcw[j] / Uce[npw] * u_eta;


	// South boundary
	j = 1;
	for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		npn = np + nx;

		for (k = 0; k <= 9 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][7] = 1.0;
		a[np][4] = 1.0;

		aux = sqrt((u[npn] * u[npn] + v[npn] * v[npn]) / (xk[np] * xk[np] + yk[np] * yk[np]))
			* 1.0*sign(u[npn] * xk[np] + v[npn] * yk[np]);

			b[np] = 2.0 * aux * xk[np];
	}


	// North boundary
	j = ny;
	for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		for (k = 0; k <= 9 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][4] = 1.0;
		a[np][1] = 1.0;

		b[np] = 2.0*UF;
	}

	
	// Corners (extrapolation)
	// SW
	i = 0;			// already fixed loop - verify
	j = 1;
	np = nx * (j - 1) + i;
	npn = np + nx;
	npe = np + 1;
	npne = npn + 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][4] = 1.0;
		
	b[np] = (u[npn] + u[npe] + u[npne]) / 3.0;


	// SE
	i = nx-1;		// already fixed loop - verify
	j = 1;
	np = nx * (j - 1) + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;
		
	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][4] = 1.0;
		
	b[np] = (u[npn] + u[npw] + u[npnw]) / 3.0;

	// NW
	i = 0;			// already fixed loop - verify
	j = ny;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npe = np + 1;
	npse = nps + 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][4] = 1.0;

	b[np] = (u[nps] + u[npe] + u[npse]) / 3.0;

	// NE
	i = nx - 1;			// already fixed loop - verify
	j = ny;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npw = np - 1;
	npsw = nps - 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][4] = 1.0;

	b[np] = (u[nps] + u[npw] + u[npsw]) / 3.0;

	coeffs_source result;
	result.a = a;
	result.b = b;
	return result;
}

coeffs_source set_bcv(int nx, int ny, double* xk, double* yk, double* u, double* v, double* Uce, double* Vcw, double** a, double* b)
{
	// Calculates the coefficients and source of the linear system of v for the fictitious volumes based on the boundary conditions.
	//double** a(nx*ny, std::vector<double>(9));
	//double* b(nx*ny);

	// Auxiliary variables
	int i;
	int j;
	int k;
	int np;
	int npw;
	int npe;
	int nps;
	int npn;
	int npse;
	int npsw;
	int npne;
	int npnw;

	int npnnw;
	int npssw;

	double aux;
	double v_eta;

	// West boundary
	i = 0;			// already fixed loop - verify

	for (j = 2; j <= ny - 1; j += 1)			// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		for (k = 0; k <= 9 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][4] = 1.0;
		a[np][5] = 1.0;
		b[np] = 0.0;
	}

	// East boundary
	i = nx-1;			// already fixed loop - verify

	j = 2;

	np = nx * (j - 1) + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;
	npnnw = npnw + nx;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}
	a[np][3] = -1.0;
	a[np][4] = 1.0;

	v_eta = 0.5 * (-v[npnnw] + 4.0 * v[npnw] - 3.0 * v[npw]);

	b[np] = -Vcw[j] / Uce[npw] * v_eta;

	for (j = 3; j <= ny - 2; j += 1)			// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npn = np + nx;
		npw = np - 1;
		npsw = nps - 1;
		npnw = npn - 1;
		
		for (k = 0; k <= 9 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][3] = -1.0;
		a[np][4] = 1.0;

		v_eta = 0.5 * (v[npnw] - v[npsw]);

		b[np] = -Vcw[j] / Uce[npw] * v_eta;

	}

	j = ny - 1;

	np = nx * (j - 1) + i;
	nps = np - nx;
	npw = np - 1;
	npsw = nps - 1;
	npssw = npsw - nx;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}
	a[np][3] = -1.0;
	a[np][4] = 1.0;

	v_eta = 0.5 * (3.0 * v[npw] - 4.0 * v[npsw] + v[npssw]);

	b[np] = -Vcw[j] / Uce[npw] * v_eta;

	// South boundary
	j = 1;
	for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		npn = np + nx;

		for (k = 0; k <= 9 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][4] = 1.0;
		a[np][7] = 1.0;

		aux = sqrt((u[npn] * u[npn] + v[npn] * v[npn]) / (xk[np] * xk[np] + yk[np] * yk[np]))
			* 1.0*sign(u[npn] * xk[np] + v[npn] * yk[np]);

		b[np] = 2.0 * aux * yk[np];
	}

	// North boundary
	j = ny;
	for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
	{
		np = nx * (j - 1) + i;

		for (k = 0; k <= 9 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][4] = 1.0;
		a[np][1] = 1.0;

		b[np] = 0.0;
	}
	
	// Corners (extrapolation)
	// SW
	i = 0;				// already fixed loop - verify
	j = 1;
	np = nx * (j - 1) + i;
	npn = np + nx;
	npe = np + 1;
	npne = npn + 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][4] = 1.0;

	b[np] = (v[npn] + v[npe] + v[npne]) / 3.0;


	// SE
	i = nx - 1;				// already fixed loop - verify
	j = 1;
	np = nx * (j - 1) + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][4] = 1.0;

	b[np] = (v[npn] + v[npw] + v[npnw]) / 3.0;

	// NW
	i = 0;				// already fixed loop - verify
	j = ny;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npe = np + 1;
	npse = nps + 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][4] = 1.0;

	b[np] = (v[nps] + v[npe] + v[npse]) / 3.0;

	// NE
	i = nx - 1;			// already fixed loop - verify
	j = ny;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npw = np - 1;
	npsw = nps - 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][4] = 1.0;

	b[np] = (v[nps] + v[npw] + v[npsw]) / 3.0;

	coeffs_source result;
	result.a = a;
	result.b = b;
	return result;
}

coeffs_source set_bcT(int nx, int ny, double TF, double* Uce, double* Vcw, double* T, double* alphae, double* betae, double* betan, double* gamman, double** a, double* b)
{
	// Calculates the coefficients  and source of the linear system of T for the fictitious volumes based on the boundary conditions.
	//double** a(nx*ny, std::vector<double>(9));
	//double* b(nx*ny);

	// Auxiliary variables
	int i;
	int j;
	int k;
	int np;
	int npe;
	int npw;
	int nps;
	int npn;
	int npse;
	int npsw;
	int npne;
	int npnw;
	int npnn;
	int npnne;
	int npss;
	int npsse;
	int npnnw;
	int npssw;
	double Te;
	double Tk;

	// West boundary
	i = 0;				// already fixed loop - verify
	j = 2;
	np = nx * (j - 1) + i;
	npn = np + nx;
	npe = np + 1;
	npne = npn + 1;
	npnn = npn + nx;
	npnne = npnn + 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}
	a[np][4] = 1.0;
	a[np][5] = -1.0;

	Te = 0.5 * (-T[npnne] + 4.0 * T[npne] - 3.0 * T[npe]);

	b[np] = -betae[np] / alphae[np] * Te;

	for (j = 3; j <= ny - 2; j += 1)			// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npn = np + nx;
		npse = nps + 1;
		npne = npn + 1;

		for (k = 0; k <= 9 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][4] = 1.0;
		a[np][5] = -1.0;

		Te = 0.5 * (T[npne] - T[npse]);

		b[np] = -betae[np] / alphae[np] * Te;
	}

	j = ny - 1;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npe = np + 1;
	npse = nps + 1;
	npss = nps - nx;
	npsse = npss + 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}
	a[np][4] = 1.0;
	a[np][5] = -1.0;

	Te = 0.5 * (T[npsse] - 4.0 * T[npse] + 3.0 * T[npe]);

	b[np] = -betae[np] / alphae[np] * Te;

	// East boundary
	i = nx - 1;				// already fixed loop - verify
	j = 2;
	np = nx * (j - 1) + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;
	npnnw = npnw + nx;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}
	a[np][3] = -1.0;
	a[np][4] = 1.0;

	Te = 0.5 * (-T[npnnw] + 4.0 * T[npnw] - 3.0 * T[npw]);

	b[np] = -Vcw[j] / Uce[npw] * Te;

	for (j = 3; j <= ny - 2; j += 1)				// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npn = np + nx;
		npw = np - 1;
		npsw = nps - 1;
		npnw = npn - 1;

		for (k = 0; k <= 9 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][3] = -1.0;
		a[np][4] = 1.0;

		Te = 0.5 * (T[npnw] - T[npsw]);

		b[np] = -Vcw[j] / Uce[npw] * Te;
	}

	j = ny - 1;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npw = np - 1;
	npsw = nps - 1;
	npssw = npsw - nx;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}
	a[np][3] = -1.0;
	a[np][4] = 1.0;

	Te = 0.5 * (3.0 * T[npw] - 4.0 * T[npsw] + T[npssw]);

	b[np] = -Vcw[j] / Uce[npw] * Te;


	// South boundary
	j = 1;

	for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		npn = np + nx;
		npw = np - 1;
		npe = np + 1;
		npnw = npn - 1;
		npne = npn + 1;

		for (k = 0; k <= 9 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][4] = 1.0;
		a[np][7] = -1.0;

		Tk = 0.25 * (T[npe] + T[npne] - T[npw] - T[npnw]);

		b[np] = -betan[np] / gamman[np] * Tk;
	}

	// North boundary
	j = ny;

	for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		for (k = 0; k <= 9 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][4] = 1.0;
		a[np][1] = 1.0;
		b[np] = 2.0 * TF;
	}

	// Corners (extrapolation)
	// SW
	i = 0;				// already fixed loop - verify
	j = 1;
	np = nx * (j - 1) + i;
	npn = np + nx;
	npe = np + 1;
	npne = npn + 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][4] = 1.0;

	b[np] = (T[npn] + T[npe] + T[npne]) / 3.0;


	// SE
	i = nx - 1;				// already fixed loop - verify
	j = 1;
	np = nx * (j - 1) + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][4] = 1.0;

	b[np] = (T[npn] + T[npw] + T[npnw]) / 3.0;

	// NW
	i = 0;			// already fixed loop - verify
	j = ny;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npe = np + 1;
	npse = nps + 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][4] = 1.0;

	b[np] = (T[nps] + T[npe] + T[npse]) / 3.0;

	// NE
	i = nx - 1;			// already fixed loop - verify
	j = ny;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npw = np - 1;
	npsw = nps - 1;

	for (k = 0; k <= 9 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][4] = 1.0;

	b[np] = (T[nps] + T[npw] + T[npsw]) / 3.0;

	coeffs_source result;
	result.a = a;
	result.b = b;
	return result;
}

coeffs_source set_bcp(int nx, int ny, double* alphae, double* betae, double* betan, double* gamman, double* Uce, double* Vcw, double* p, double** a, double* b)
{
	// Calculates the coefficients and source of the linear system of the pressure deviation for the fictitious volumes based on the boundary conditions
	//double** a(nx*ny, std::vector<double>(5));
	//double* b(nx*ny);

	// Auxiliary variables
	int i;
	int j;
	int k;
	int np;
	int npe;
	int npw;
	int nps;
	int npn;
	int npse;
	int npsw;
	int npne;
	int npnw;
	int npnee;
	int npnww;
	int npnne;
	int npsse;
	int npnnw;
	int npssw;
	double pk;
	double pe;

	// West boundary
	i = 0;				// already fixed loop - verify
	j = 2;

	np = nx * (j - 1) + i;
	npn = np + nx;
	npe = np + 1;
	npne = npn + 1;
	npnne = npne + nx;

	for (k = 0; k <= 5 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][2] = 1.0;

	a[np][3] = -1.0;

	pe = (-p[npnne] + 4.0 * p[npne] - 3.0 * p[npe]) / 2.0;

	b[np] = p[npe] - p[np] - betae[np] / alphae[np] * pe;

	for (j = 3; j <= ny - 2; j += 1)			// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npn = np + nx;
		npe = np + 1;
		npse = nps + 1;
		npne = npn + 1;

		for (k = 0; k <= 5 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}

		a[np][2] = 1.0;

		a[np][3] = -1.0;

		pe = (p[npne] - p[npse]) / 2.0;

		b[np] = p[npe] - p[np] - betae[np] / alphae[np] * pe;
	}

	j = ny - 1;

	np = nx * (j - 1) + i;
	nps = np - nx;
	npe = np + 1;
	npse = nps + 1;
	npsse = npse - nx;

	for (k = 0; k <= 5 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][2] = 1.0;

	a[np][3] = -1.0;

	pe = (3.0 * p[npe] - 4.0 * p[npse] + p[npsse]) / 2.0;

	b[np] = p[npe] - p[np] - betae[np] / alphae[np] * pe;


	// East boundary
	i = nx - 1;			// already fixed loop - verify
	j = 2;

	np = nx * (j - 1) + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;
	npnnw = npnw + nx;

	for (k = 0; k <= 5 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}
	a[np][1] = -1.0;
	a[np][2] = 1.0;

	pe = (-p[npnnw] + 4.0 * p[npnw] - 3.0 * p[npw]) / 2.0;

	b[np] = p[npw] - p[np] - Vcw[j] / Uce[npw] * pe;

	for (j = 3; j <= ny - 2; j += 1)
	{

		np = nx * (j - 1) + i;
		nps = np - nx;
		npn = np + nx;
		npw = np - 1;
		npsw = nps - 1;
		npnw = npn - 1;

		for (k = 0; k <= 5 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}
		a[np][1] = -1.0;
		a[np][2] = 1.0;

		pe = (p[npnw] - p[npsw]) / 2.0;

		b[np] = p[npw] - p[np] - Vcw[j] / Uce[npw] * pe;
	}

	j = ny - 1;

	np = nx * (j - 1) + i;
	nps = np - nx;
	npw = np - 1;
	npsw = nps - 1;
	npssw = npsw - nx;

	for (k = 0; k <= 5 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}
	a[np][1] = -1.0;
	a[np][2] = 1.0;

	pe = (3.0 * p[npw] - 4.0 * p[npsw] + p[npssw]) / 2.0;

	b[np] = p[npw] - p[np] - Vcw[j] / Uce[npw] * pe;


	// South boundary
	j = 1;
	i = 1;				// already fixed loop - verify

	np = nx * (j - 1) + i;
	npn = np + nx;
	npne = npn + 1;
	npnee = npne + 1;

	for (k = 0; k <= 5 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][2] = 1.0;

	a[np][4] = -1.0;

	pk = (-p[npnee] + 4.0 * p[npne] - 3.0 * p[npn]) / 2.0;

	b[np] = p[npn] - p[np] - betan[np] / gamman[np] * pk;

	for (i = 2; i <= nx - 2 - 1; i += 1)			// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		npn = np + nx;
		npnw = npn - 1;
		npne = npn + 1;

		for (k = 0; k <= 5 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}

		a[np][2] = 1.0;

		a[np][4] = -1.0;

		pk = (p[npne] - p[npnw]) / 2.0;

		b[np] = p[npn] - p[np] - betan[np] / gamman[np] * pk;
	}

	i = nx - 1 - 1;				// already fixed loop - verify

	np = nx * (j - 1) + i;
	npn = np + nx;
	npnw = npn - 1;
	npnww = npnw - 1;

	for (k = 0; k <= 5 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][2] = 1.0;

	a[np][4] = -1.0;

	pk = (3.0 * p[npn] - 4.0 * p[npnw] + p[npnww]) / 2.0;

	b[np] = p[npn] - p[np] - betan[np] / gamman[np] * pk;


	// North boundary
	j = ny;

	for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
	{

		np = nx * (j - 1) + i;

		for (k = 0; k <= 5 - 1; k += 1)
		{
			a[np][k] = 0.0;
		}

		a[np][2] = 1.0;

		a[np][0] = 1.0;

		b[np] = 0.0;
	}


	// Corners

	// SW
	i = 0;				// already fixed loop - verify
	j = 1;
	np = nx * (j - 1) + i;
	npn = np + nx;
	npe = np + 1;
	npne = npn + 1;

	for (k = 0; k <= 5 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][2] = 1.0;

	b[np] = 0.0;


	// SE
	i = nx - 1;				// already fixed loop - verify
	j = 1;
	np = nx * (j - 1) + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;

	for (k = 0; k <= 5 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][2] = 1.0;

	b[np] = 0.0;


	// NW
	i = 0;				// already fixed loop - verify
	j = ny;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npe = np + 1;
	npse = nps + 1;

	for (k = 0; k <= 5 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][2] = 1.0;

	b[np] = 0.0;


	// NE
	i = nx - 1;			// already fixed loop - verify
	j = ny;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npw = np - 1;
	npsw = nps - 1;

	for (k = 0; k <= 5 - 1; k += 1)
	{
		a[np][k] = 0.0;
	}

	a[np][2] = 1.0;

	b[np] = 0.0;

	coeffs_source result;
	result.a = a;
	result.b = b;
	return result;
}

double* get_Vcw(int nx, int ny, double* xe, double* xk, double* ue, double* ve, double* Vcw)
{
	// Calculates the contravariant velocity V on the west face of the fictitious volumes of the east boundary

	//Inner variables
	int i;
	int j;
	int np;

	i = nx - 1 - 1;				// already fixed loop - verify
	for (j = 2; j <= ny - 1; j += 1)
	{
		np = nx * (j - 1) + i;
		Vcw[j] = ve[np] * xk[np] - ue[np] * xe[np];
	}

	return Vcw;
}

velocity get_u_v_extrapolation_to_fictitious(int nx, int ny, double UF, double* xk, double* yk, double* alphae, double* betae, double* Uce, double* Vcw, double* u, double* v)
{
	// Extrapolates the corrected velocities of the real volumes to the fictitious volumes based on the boundary conditions.
	
	// Auxiliary variables
	int i, j, np, npe, npw, nps, npn;
	int npssw, npnnw, npsw, npnw, npsse, npnne, npss, npnn, npse, npne;
	double aux, u_eta, v_eta;

	// West boundary ( only for u )
	i = 0;				// already fixed loop - verify
	j = 2;

	np = nx * (j - 1) + i;
	npn = np + nx;
	npe = np + 1;
	npne = npn + 1;
	npnn = npn + nx;
	npnne = npnn + 1;

	u_eta = 0.5 * (-u[npnne] + 4.0 * u[npne] - 3.0 * u[npe]);

	u[np] = u[npe] - betae[np] / alphae[np] * u_eta;

	for (j = 3; j <= ny - 2; j += 1)
	{
		np = nx * (j - 1) + i;
		npe = np + 1;
		nps = np - nx;
		npn = np + nx;
		npse = nps + 1;
		npne = npn + 1;

		u_eta = 0.5 * (u[npne] - u[npse]);

		u[np] = u[npe] - betae[np] / alphae[np] * u_eta;
	}

	j = ny - 1;

	np = nx * (j - 1) + i;
	nps = np - nx;
	npe = np + 1;
	npse = nps + 1;
	npss = nps - nx;
	npsse = npss + 1;

	u_eta = 0.5 * (u[npsse] - 4.0 * u[npse] + 3.0 * u[npe]);

	u[np] = u[npe] - betae[np] / alphae[np] * u_eta;


	// West boundary ( only for v )
	i = 0;				// already fixed loop - verify

	for (j = 2; j <= ny - 1; j += 1)
	{
		np = nx * (j - 1) + i;
		npe = np + 1;

		v[np] = -v[npe];
	}


	// East boundary
	i = nx - 1;				// already fixed loop - verify
	j = 2;

	np = nx * (j - 1) + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;
	npnnw = npnw + nx;

	u_eta = 0.5 * (-u[npnnw] + 4.0 * u[npnw] - 3.0 * u[npw]);
	v_eta = 0.5 * (-v[npnnw] + 4.0 * v[npnw] - 3.0 * v[npw]);

	u[np] = u[npw] - Vcw[j] / Uce[npw] * u_eta;
	v[np] = v[npw] - Vcw[j] / Uce[npw] * v_eta;

	for (j = 3; j <= ny - 2; j += 1)
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npn = np + nx;
		npw = np - 1;
		npsw = nps - 1;
		npnw = npn - 1;

		u_eta = 0.5 * (u[npnw] - u[npsw]);
		v_eta = 0.5 * (v[npnw] - v[npsw]);

		u[np] = u[npw] - Vcw[j] / Uce[npw] * u_eta;
		v[np] = v[npw] - Vcw[j] / Uce[npw] * v_eta;
	}

	j = ny - 1;

	np = nx * (j - 1) + i;
	nps = np - nx;
	npw = np - 1;
	npsw = nps - 1;
	npssw = npsw - nx;

	u_eta = 0.5 * (3.0 * u[npw] - 4.0 * u[npsw] + u[npssw]);
	v_eta = 0.5 * (3.0 * v[npw] - 4.0 * v[npsw] + v[npssw]);

	u[np] = u[npw] - Vcw[j] / Uce[npw] * u_eta;
	v[np] = v[npw] - Vcw[j] / Uce[npw] * v_eta;


	// South boundary
	j = 1;

	for (i = 1; i <= nx - 1 - 1; i += 1)				// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		npn = np + nx;

		aux = sqrt((u[npn] * u[npn] + v[npn] * v[npn]) / (xk[np] * xk[np] + yk[np] * yk[np]))
				*1.0*sign(u[npn] * xk[np] + v[npn] * yk[np]);

		u[np] = 2.0 * aux * xk[np] - u[npn];

		v[np] = 2.0 * aux * yk[np] - v[npn];
	}

	// North boundary (far field)
	j = ny;

	for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		nps = np - nx;

		u[np] = 2.0 * UF - u[nps];

		v[np] = -v[nps];
	}

	u = get_extrapolation_to_corners(nx, ny, u);

	v = get_extrapolation_to_corners(nx, ny, v);

	velocity result;
	result.u = u;
	result.v = v;
	return result;
}

double* get_T_extrapolation_to_fictitious(int nx, int ny, double TF, double* alphae, double* betae, double* betan, double* gamman, double* Uce, double* Vcw, double* T)
{
	// Extrapolates the temperature of the real volumes to the fictitious volumes based on the boundary conditions.

	// Auxiliary variables
	int i, j, np, npe, npw, nps, npn, npnw, npne, npse;
	int npnn, npnne, npss, npsse, npsw, npssw, npnnw;
	double Te, Tk;

	// West boundary (frontal symmetry line)
	i = 0;					// already fixed loop - verify
	j = 2;

	np = nx * (j - 1) + i;
	npn = np + nx;
	npe = np + 1;
	npne = npn + 1;
	npnn = npn + nx;
	npnne = npnn + 1;

	Te = 0.5 * (-T[npnne] + 4.0 * T[npne] - 3.0 * T[npe]);

	T[np] = T[npe] - betae[np] / alphae[np] * Te;

	for (j = 3; j <= ny - 2; j += 1)
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npn = np + nx;
		npe = np + 1;
		npse = nps + 1;
		npne = npn + 1;

		Te = 0.5 * (T[npne] - T[npse]);

		T[np] = T[npe] - betae[np] / alphae[np] * Te;
	}

	j = ny - 1;

	np = nx * (j - 1) + i;
	nps = np - nx;
	npe = np + 1;
	npse = nps + 1;
	npss = nps - nx;
	npsse = npss + 1;

	Te = 0.5 * (T[npsse] - 4.0 * T[npse] + 3.0 * T[npe]);

	T[np] = T[npe] - betae[np] / alphae[np] * Te;


	// East face
	i = nx - 1;				// already fixed loop - verify
	j = 2;

	np = nx * (j - 1) + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;
	npnnw = npnw + nx;

	Te = 0.5 * (-T[npnnw] + 4.0 * T[npnw] - 3.0 * T[npw]);

	T[np] = T[npw] - Vcw[j] / Uce[npw] * Te;

	for (j = 3; j <= ny - 2; j += 1)
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npn = np + nx;
		npw = np - 1;
		npsw = nps - 1;
		npnw = npn - 1;

		Te = 0.5 * (T[npnw] - T[npsw]);

		T[np] = T[npw] - Vcw[j] / Uce[npw] * Te;
	}

	j = ny - 1;

	np = nx * (j - 1) + i;
	nps = np - nx;
	npn = np + nx;
	npw = np - 1;
	npsw = nps - 1;
	npssw = npsw - nx;

	Te = 0.5 * (3.0 * T[npw] - 4.0 * T[npsw] + T[npssw]);

	T[np] = T[npw] - Vcw[j] / Uce[npw] * Te;


	// South boundary
	j = 1;

	for (i = 1; i <= nx - 1 - 1; i += 1)				// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		npn = np + nx;
		npw = np - 1;
		npe = np + 1;
		npnw = npn - 1;
		npne = npn + 1;

		Tk = 0.25 * (T[npe] + T[npne] - T[npw] - T[npnw]);

		T[np] = T[npn] - betan[np] / gamman[np] * Tk;
	}

	// North boundary (far field);
	j = ny;

	for (i = 1; i <= nx - 1 - 1; i += 1)				// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		nps = np - nx;

		T[np] = 2.0 * TF - T[nps];
	}

	T = get_extrapolation_to_corners(nx, ny, T);

	return T;
}

double* get_p_extrapolation_to_fictitious(int nx, int ny, double PF, double* alphae, double* betae, double* betan, double* gamman, double* Uce, double* Vcw, double* p)
{
	// Extrapolates the pressure of the real volumes to the fictitious volumes based on the boundary conditions.

	// Auxiliary variables
	int i, j, np, npe, npw, nps, npn;
	int npne, npnw, npse, npsw, npnww, npnee, npnne, npsse, npnnw, npssw;
	double pk, pe;

	// West boundary
	i = 0;				// already fixed loop - verify
	j = 2;

	np = nx * (j - 1) + i;
	npn = np + nx;
	npe = np + 1;
	npne = npn + 1;
	npnne = npne + nx;

	pe = (-p[npnne] + 4.0 * p[npne] - 3.0 * p[npe]) / 2.0;

	p[np] = p[npe] - betae[np] / alphae[np] * pe;

	for (j = 3; j <= ny - 2; j += 1)
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npn = np + nx;
		npe = np + 1;
		npse = nps + 1;
		npne = npn + 1;

		pe = (p[npne] - p[npse]) / 2.0;

		p[np] = p[npe] - betae[np] / alphae[np] * pe;
	}

	j = ny - 1;

	np = nx * (j - 1) + i;
	nps = np - nx;
	npe = np + 1;
	npse = nps + 1;
	npsse = npse - nx;

	pe = (3.0 * p[npe] - 4.0 * p[npse] + p[npsse]) / 2.0;

	p[np] = p[npe] - betae[np] / alphae[np] * pe;


	// East boundary
	i = nx - 1;				// already fixed loop - verify
	j = 2;

	np = nx * (j - 1) + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;
	npnnw = npnw + nx;

	pe = (-p[npnnw] + 4.0 * p[npnw] - 3.0 * p[npw]) / 2.0;

	p[np] = p[npw] - Vcw[j] / Uce[npw] * pe;

	for (j = 3; j <= ny - 2; j += 1)
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npn = np + nx;
		npw = np - 1;
		npsw = nps - 1;
		npnw = npn - 1;

		pe = (p[npnw] - p[npsw]) / 2.0;

		p[np] = p[npw] - Vcw[j] / Uce[npw] * pe;
	}

	j = ny - 1;

	np = nx * (j - 1) + i;
	nps = np - nx;
	npw = np - 1;
	npsw = nps - 1;
	npssw = npsw - nx;

	pe = (3.0 * p[npw] - 4.0 * p[npsw] + p[npssw]) / 2.0;

	p[np] = p[npw] - Vcw[j] / Uce[npw] * pe;


	// South boundary
	j = 1;
	i = 1;				// already fixed loop - verify

	np = nx * (j - 1) + i;
	npn = np + nx;
	npne = npn + 1;
	npnee = npne + 1;

	pk = (-p[npnee] + 4.0 * p[npne] - 3.0 * p[npn]) / 2.0;

	p[np] = p[npn] - betan[np] / gamman[np] * pk;

	for (i = 2; i <= nx - 2 - 1; i += 1)				// already fixed loop - verify
	{

		np = nx * (j - 1) + i;
		npn = np + nx;
		npnw = npn - 1;
		npne = npn + 1;

		pk = (p[npne] - p[npnw]) / 2.0;

		p[np] = p[npn] - betan[np] / gamman[np] * pk;
	}

	i = nx - 1 - 1;					// already fixed loop - verify

	np = nx * (j - 1) + i;
	npn = np + nx;
	npnw = npn - 1;
	npnww = npnw - 1;

	pk = (3.0 * p[npn] - 4.0 * p[npnw] + p[npnww]) / 2.0;

	p[np] = p[npn] - betan[np] / gamman[np] * pk;


	// North boundary (far field)
	j = ny;

	for (i = 1; i <= nx - 1 - 1; i += 1)				// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		nps = np - nx;

		p[np] = 2.0 * PF - p[nps];
	}

	p = get_extrapolation_to_corners(nx, ny, p);

	return p;
}

double* get_extrapolation_to_corners(int nx, int ny, double* f)
{
	// Extrapolation of a field to the fictitious volumes of the corners.

	// Auxiliary variables
	int i, j, np, npe, npw, npn, nps, npne, npnw, npse, npsw;

	// Corners (extrapolation)
	// SW
	i = 0;					// already fixed loop - verify
	j = 1;
	np = nx * (j - 1) + i;
	npn = np + nx;
	npe = np + 1;
	npne = npn + 1;

	f[np] = (f[npn] + f[npe] + f[npne]) / 3.0;

	// SE
	i = nx - 1;					// already fixed loop - verify		
	j = 1;
	np = nx * (j - 1) + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;

	f[np] = (f[npn] + f[npw] + f[npnw]) / 3.0;

	// NW
	i = 0;							// already fixed loop - verify
	j = ny;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npe = np + 1;
	npse = nps + 1;

	f[np] = (f[nps] + f[npe] + f[npse]) / 3.0;

	// NE
	i = nx - 1;						// already fixed loop - verify
	j = ny;
	np = nx * (j - 1) + i;
	nps = np - nx;
	npw = np - 1;
	npsw = nps - 1;

	f[np] = (f[nps] + f[npw] + f[npsw]) / 3.0;

	return f;
}

simplec_coeffs get_boundary_simplec_coefficients(int nx, int ny, double* xe, double* ye, double* due, double* dve, double* dun, double* dvn, double* de, double* dn)
{
	// Calculates the SIMPLEC coefficients at the boundary faces.

	// Auxiliary variables
	int i, j, np, npe, npw, npn;

	// West boundary
	i = 0;						// already fixed loop - verify

	for (j = 2; j <= ny - 1; j += 1)
	{
		np = nx * (j - 1) + i;
		npe = np + 1;

		due[np] = due[npe];

		dve[np] = 0.0;

		de[np] = 0.0;
	}

	// East boundary
	i = nx - 1 - 1;						// already fixed loop - verify

	for (j = 2; j <= ny - 1; j += 1)
	{
		np = nx * (j - 1) + i;
		npw = np - 1;

		due[np] = due[npw];

		dve[np] = dve[npw];

		de[np] = ye[np] * due[np] + xe[np] * dve[np];
	}

	// South boundary
	j = 1;

	for (i = 1; i <= nx - 1 - 1; i += 1)					// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		npn = np + nx;

		dun[np] = dun[npn];

		dvn[np] = dvn[npn];

		dn[np] = 0.0;
	}

	// North boundary
	j = ny - 1;

	for (i = 1; i <= nx - 1 - 1; i += 1)						// already fixed loop - verify
	{
		np = nx * (j - 1) + i;

		dun[np] = 0.0;

		dvn[np] = 0.0;

		dn[np] = 0.0;
	}

	simplec_coeffs result;
	result.due = due;
	result.dve = dve;
	result.dun = dun;
	result.dvn = dvn;
	result.de = de;
	result.dn = dn;
	return result;
}

velocity_face get_velocities_at_boundary_faces(int nx, int ny, double UF, double* xe, double* ye, double* xk, double* yk, double* u, double* v,
												double* ue, double* ve, double* un, double* vn, double* Uce, double* Vcn)
{
	// Calculates the velocities at boundary faces based on the boundary conditions.

	// Variables
	//double* ue(nx*ny);
	//double* ve(nx*ny);
	//double* un(nx*ny);
	//double* vn(nx*ny);
	//double* Uce(nx*ny);
	//double* Vcnb(nx*ny);

	// Axuliary variables
	int i, j, np, npe, npn;
	double aux;

	// West boundary ( frontal symmetry line )
	i = 0;							// already fixed loop - verify

	for (j = 2; j <= ny - 1; j += 1)
	{
		np = nx * (j - 1) + i;
		npe = np + 1;

		ue[np] = (u[np] + u[npe]) / 2.0;

		ve[np] = 0.0;

		Uce[np] = 0.0; // There is no mass flow through the symmetry line
	}

	// East boundary
	i = nx - 1 - 1;						// already fixed loop - verify

	for (j = 2; j <= ny - 1; j += 1)
	{
		np = nx * (j - 1) + i;
		npe = np + 1;

		ue[np] = (u[np] + u[npe]) / 2.0;

		ve[np] = (v[np] + v[npe]) / 2.0;

		Uce[np] = ue[np] * ye[np] - ve[np] * xe[np];
	}

	// South boundary (body surface)
	j = 1;

	for (i = 1; i <= nx - 1 - 1; i += 1)						// already fixed loop - verify
	{
		np = nx * (j - 1) + i;
		npn = np + nx;

		aux = sqrt((u[npn] * u[npn] + v[npn]*v[npn]) / (xk[np]*xk[np] + yk[np]*yk[np]))
			* 1.0*sign(u[npn] * xk[np] + v[npn] * yk[np]);

		un[np] = aux * xk[np];

		vn[np] = aux * yk[np];

		Vcn[np] = 0.0; // There is no mass flow through the body surface
	}

	// North boundary (far field surface)
	j = ny - 1;

	for (i = 1; i <= nx - 1 - 1; i += 1)						// already fixed loop - verify
	{
		np = nx * (j - 1) + i;

		un[np] = UF;

		vn[np] = 0.0;

		Vcn[np] = -un[np] * yk[np];
	}

	velocity_face result;
	result.ue = ue;
	result.ve = ve;
	result.un = un;
	result.vn = vn;
	result.Uce = Uce;
	result.Vcn = Vcn;
	return result;
}

initial_conditions get_initial_conditions(int nx, int ny, double PF, double TF, double Rg, double GF, double MF, double UF, double* xe, double* ye, double* xk, double* yk, double* alphae, double* betae, double* betan, double* gamman,
	double* p, double* pl, double* T, double* ro, double* roe, double* ron, double* u, double* v, double* ue, double* ve, double* un, double* vn, double* Uce, double* Vcn, double* Vcw)
{
	// Prescribes an initial conditions for some variables.

	// Auxiliary variables
	int i, j, np;
	double ROF = PF / (Rg * TF);			// Free stream density

	// Variables
	//double UF = MF * sqrt(GF * Rg * TF);	// Free-stream speed
	//double* p;
	//double* pl;
	//double* T;
	//double* ro;
	//double* roe;
	//double* ron;
	//double* u;
	//double* v;
	//double* ue;
	//double* ve;
	//double* un;
	//double* vn;
	//double* Uce;
	//double* Vcn;
	//double* Vcw;

	//p = new double[nx*ny];
	//pl = new double[nx*ny]();
	//T = new double[nx*ny];
	//ro = new double[nx*ny];
	//roe = new double[nx*ny];
	//ron = new double[nx*ny];
	//u = new double[nx*ny];
	//v = new double[nx*ny]();
	//ue = new double[nx*ny];
	//ve = new double[nx*ny]();
	//un = new double[nx*ny];
	//vn = new double[nx*ny]();
	//Uce = new double[nx*ny];
	//Vcn = new double[nx*ny];
	//Vcw = new double[ny]();

	UF = MF * sqrt(GF * Rg * TF);

	//std::fill_n(p, nx*ny, init);
	for (i = 0; i < nx*ny; i++)
	{
		p[i] = PF;
		T[i] = TF;
		ro[i] = ROF;
		roe[i] = ROF;
		ron[i] = ROF;
		u[i] = UF;
		ue[i] = UF;
		un[i] = UF;
		Uce[i] = UF;
	}
		
	velocity velocities = get_u_v_extrapolation_to_fictitious(nx, ny, UF, xk, yk, alphae, betae, Uce, Vcw, u, v); // Last two are inout
	u = velocities.u;
	v = velocities.v;

	// U and V on the interface between real volumes
	for (i = 1; i <= nx - 2 - 1; i += 1)					// already fixed loop - verify
	{
		for (j = 2; j <= ny - 1; j += 1)
		{
			np = nx * (j - 1) + i;

			Uce[np] = ue[np] * ye[np] - ve[np] * xe[np];
		}
	}

	for (i = 1; i <= nx - 1 - 1; i += 1)					// already fixed loop - verify
	{
		for (j = 2; j <= ny - 2; j += 1)
		{
			np = nx * (j - 1) + i;

			Vcn[np] = vn[np] * xk[np] - un[np] * yk[np];
		}
	}

	// u, v, U and V  on the boundary real faces
	velocity_face velocities_boundary = get_velocities_at_boundary_faces(nx, ny, UF, xe, ye, xk, yk, u, v, ue, ve, un, vn, Uce, Vcn);
	ue = velocities_boundary.ue;
	ve = velocities_boundary.ve;
	un = velocities_boundary.un;
	vn = velocities_boundary.vn;
	Uce = velocities_boundary.Uce;
	Vcn = velocities_boundary.Vcn;

	// Calculates the contravariant velocity V on the west face
	// of the fictitious volumes of the east boundary
	Vcw = get_Vcw(nx, ny, xe, xk, ue, ve, Vcw);


	// Values of the variables at fictitious nodes

	// Extrapolation of T according the boundary conditions
	T = get_T_extrapolation_to_fictitious(nx, ny, TF, alphae, betae, betan, gamman, Uce, Vcw, T); // Last is inout


	// Extrapolation of p according the boundary conditions
	p = get_p_extrapolation_to_fictitious(nx, ny, PF, alphae, betae, betan, gamman, Uce, Vcw, p); // Last is inout


	// Extrapolation of ro according to the state equation
	for (i = 0; i <= nx * ny - 1; i += 1)
	{
		ro[i] = p[i] / (Rg * T[i]);
	}

	initial_conditions result;
	result.UF = UF;
	result.p = p;
	result.pl = pl;
	result.T = T;
	result.ro = ro;
	result.roe = roe;
	result.ron = ron;
	result.velocity.u = u;
	result.velocity.v = v;
	result.velocity_face.ue = ue;
	result.velocity_face.ve = ve;
	result.velocity_face.un = un;
	result.velocity_face.vn = vn;
	result.velocity_face.Uce = Uce;
	result.velocity_face.Vcn = Vcn;
	result.Vcw = Vcw;
	return result;
}