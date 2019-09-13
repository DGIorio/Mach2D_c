#include "stdafx.h"

// Contains subroutines related to the grid generation and metrics calculation.


coordinates set_grid(int kg, int nx, int ny, double a1, double* x, double* y)
{
	// Calculates the coordinates (x,y) of grid using the coordinates of the boundary.
	// There are three options of spreading of the eta lines, depending on the values of kg:
	// 1 = uniform,
	// 2 = geometric progression,
	// 3 = power law.

	coordinates result;

	switch (kg)
	{
	case 1:
		result = get_uniform_grid(nx, ny, x, y);
		x = result.x;
		y = result.y;
		break;
	case 2:
		result = get_pg_grid(a1, nx, ny, x, y);
		x = result.x;
		y = result.y;
		break;
	case 3:
		result = get_power_grid(a1, nx, ny, x, y);
		x = result.x;
		y = result.y;
		break;
	default:
		std::cout << "Grid type set to uniform grid." << std::endl;
		result = get_uniform_grid(nx, ny, x, y);
		x = result.x;
		y = result.y;
		break;
	}

	return result;
}

coordinates get_uniform_grid(int nx, int ny, double* x, double* y)
{
	// Calculates the coordinates (x,y) of the grid with a uniform distribution of the eta lines.

	int i, j, np, nps;
	double xi, xf, yi, yf, dx, dy;

	for (i = 0; i <= nx - 1 - 1; i += 1)
	{
		j = 1;
		np = nx * (j - 1) + i;
		xi = x[np];
		yi = y[np];

		j = ny - 1;
		np = nx * (j - 1) + i;
		xf = x[np];
		yf = y[np];

		dx = (xf - xi) / (ny - 2);
		dy = (yf - yi) / (ny - 2);

		for (j = 2; j <= ny - 2; j += 1)
		{
			np = nx * (j - 1) + i;
			nps = np - nx;
			x[np] = x[nps] + dx;
			y[np] = y[nps] + dy;
		}
	}

	coordinates result;
	result.x = x;
	result.y = y;
	return result;
}

coordinates get_power_grid(double a1, int nx, int ny, double* x, double* y)
{
	// Calculates the coordinates (x,y) of the grid with a power law distribution of the eta lines.

	int i, j, np;
	double a, r, xi, xf, yi, yf;

	for (i = 0; i <= nx - 1 - 1; i += 1)
	{
		j = 1;
		np = nx * (j - 1) + i;
		xi = x[np];
		yi = y[np];

		j = ny - 1;
		np = nx * (j - 1) + i;
		xf = x[np];
		yf = y[np];

		r = sqrt( (xf - xi)*(xf - xi) + (yf - yi)*(yf - yi)) / a1;

		if (r > 1.0)
		{
			a = log(r) / log((double)(ny - 2));
		}
		else
		{
			a = 1.0;
		}

		for (j = 2; j <= ny - 2; j += 1)
		{
			np = nx * (j - 1) + i;

			x[np] = (xf - xi) * std::pow((double)(j - 1) / (double)(ny - 2), a) + xi;

			y[np] = (yf - yi) * std::pow((double)(j - 1) / (double)(ny - 2), a) + yi;
		}
	}

	coordinates result;
	result.x = x;
	result.y = y;
	return result;
}

coordinates get_pg_grid(double a1, int nx, int ny, double* x, double* y)
{
	// Calculates the coordinates (x,y) of the grid with a geometric progression distribution of the eta lines.

	int i, j, np, npn;
	double q, r, t, xi, xf, yi, yf;

	for (i = 0; i <= nx - 1 - 1; i += 1)
	{
		j = 1;
		np = nx * (j - 1) + i;
		xi = x[np];
		yi = y[np];

		j = ny - 1;
		np = nx * (j - 1) + i;
		xf = x[np];
		yf = y[np];

		r = sqrt((xf - xi)*(xf - xi) + (yf - yi)*(yf - yi)) / a1;

		q = get_gp_ratio(ny - 2, r);
		t = 0.0;

		for (j = 1; j <= ny - 3; j += 1)
		{
			np = nx * (j - 1) + i;
			npn = np + nx;

			t = t + std::pow(q,(j - 1)) / r; // t of j+1

			x[npn] = xi + (xf - xi) * t;

			y[npn] = yi + (yf - yi) * t;
		}
	}

	coordinates result;
	result.x = x;
	result.y = y;
	return result;
}

double get_gp_ratio(int n, double r)
{
	// Calculates the geometric progression exponent q.
	double q;

	int it;
	double qo, fo, fo1;

	qo = 2.0;

	if ((double)(n) > r) qo = 0.5;

	it = 0;

	for (;;)
	{
		it = it + 1;

		fo = std::pow(qo,n) + r * (1.0 - qo) - 1.0;

		fo1 = (double)(n)* std::pow(qo,(n - 1)) - r;

		q = qo - fo / fo1;

		if (std::fabs(q - qo) < 1.0e-15)
		{
			break;
		}

		if (it > 10000)
		{
			std::cout << "get_gp_ratio: Maximum number of iteractions was exceeded." << std::endl;
			exit(0);
		}

		qo = q;
	}

	return q;
}

coordinates get_real_centroids_xy(int opt, int nx, int ny, double* x, double* y, double* xp, double* yp)
{
	// Calculates the coordinates (xp,yp) of the centroids of the real volumes.
	// There are two options: (opt = 1) a simple mean and (opt = 2) a weighted mean.
	
	// Centroids of real volumes
	//xp;
	//yp;

	int i, j, np, nps, npw, npsw;
	double A1, A2;

	if (opt == 1) // Simple mean
	{
		// Centroids are given by the mean of the corners
		for (j = 2; j <= ny - 1; j += 1)
		{
			for (i = 1; i <= nx - 1 - 1; i += 1)
			{
				np = nx * (j - 1) + i;
				nps = np - nx;
				npw = np - 1;
				npsw = nps - 1;

				xp[np] = (x[np] + x[npw] + x[npsw] + x[nps]) * 0.25;
				yp[np] = (y[np] + y[npw] + y[npsw] + y[nps]) * 0.25;
			}
		}
	}
	else // weighted mean
	{
		// The rectangle is divided in two triangles of area A1 and A2.
		// The centroids of each triangle is calculated by the mean of the corners.
		// The centroid of the rectangle is given by the weighted mean of the centroids of each triangle.
		for (j = 2; j <= ny - 1; j += 1)
		{
			for (i = 1; i <= nx - 1 - 1; i += 1)
			{

				np = nx * (j - 1) + i;
				nps = np - nx;
				npw = np - 1;
				npsw = nps - 1;

				A1 = 0.5 * (x[npsw] * (y[nps] - y[np])
					+ x[nps] * (y[np] - y[npsw])
					+ x[np] * (y[npsw] - y[nps]));

				A2 = 0.5 * (x[npsw] * (y[np] - y[npw])
					+ x[np] * (y[npw] - y[npsw])
					+ x[npw] * (y[npsw] - y[np]));

				xp[np] = ((x[npsw] + x[nps] + x[np]) * A1 + (x[npsw] + x[np] + x[npw]) * A2)
					/ (3.0 * (A1 + A2));

				yp[np] = ((y[npsw] + y[nps] + y[np]) * A1 + (y[npsw] + y[np] + y[npw]) * A2)
					/ (3.0 * (A1 + A2));
			}
		}
	}

	coordinates result;
	result.x = xp;
	result.y = yp;
	return result;
}

metrics get_metrics(int nx, int ny, double* x, double* y, double* xp, double* yp,
	double* xe, double* ye, double* xen, double* yen, double* xk, double* yk, double* xke, double* yke,
	double* Jp, double* Je, double* Jn, double* alphae, double* gamman, double* betae, double* betan)
{
	// Calculates the metrics, jacobians and components of the metric tensor.
	//xe - x_eta at the center of east face of volume P (m)
	//ye - y_eta at the center of east face of volume P (m)
	//xen - x_eta at the center of north face of volume P (m)
	//yen - y_eta at the center of north face of volume P (m)
	//xk - x_csi at the center of north face  of volume P (m)
	//yk - y_csi at the center of north face of volume P (m)
	//xke - x_csi at the center of east face of volume P (m)
	//yke - y_csi at the center of east face of volume P (m)
	//Jp - Jacobian at the center of volume P (1/m2)
	//Je - Jacobian at the center of east face of volume P (1/m2)
	//Jn - Jacobian at the center of north face of volume P (1/m2)
	//alphae - (metric tensor component) alpha at the center of east face of volume P (m2)
	//gamman -(metric tensor component) gamma at the center of north face of volume P (m2)
	//betae - (metric tensor component) beta  at the center of east face of volume P (m2)
	//betan -(metric tensor component) beta  at the center of north face of volume P (m2)

	int i, j, npsw, nps, npw, np, npe, npn;
	double fw, fe, fs, fn, xkp, xep, ykp, yep;

	// Derivatives relatively to eta at the center of the east face
	for (j = 2; j <= ny - 1; j += 1)
	{
		for (i = 0; i <= nx - 1 - 1; i += 1)
		{
			np = nx * (j - 1) + i;
			nps = np - nx;

			xe[np] = x[np] - x[nps];
			ye[np] = y[np] - y[nps];

			alphae[np] = xe[np]*xe[np] + ye[np]*ye[np];
		}
	}

	// Derivatives relatively to csi at the center of the north face
	for (j = 1; j <= ny - 1; j += 1)
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)
		{
			np = nx * (j - 1) + i;
			npw = np - 1;

			xk[np] = x[np] - x[npw];
			yk[np] = y[np] - y[npw];

			gamman[np] = xk[np]*xk[np] + yk[np]*yk[np];
		}
	}

	// Derivatives relatively to csi at the center of the east face (only inner real faces)
	for (j = 2; j <= ny - 1; j += 1)
	{
		for (i = 1; i <= nx - 2 - 1; i += 1)
		{
			np = nx * (j - 1) + i;
			npe = np + 1;

			xke[np] = xp[npe] - xp[np];
			yke[np] = yp[npe] - yp[np];
		}
	}

	// Derivatives relatively to csi at the center of the east face (only west boundary of the domain)
	i = 1;
	for (j = 2; j <= ny - 1; j += 1)
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npw = np - 1;
		npsw = nps - 1;

		fw = (x[npw] + x[npsw]) * 0.5;
		fe = (x[np] + x[nps]) * 0.5;

		xke[npw] = -3.0 * fw + 4.0 * xp[np] - fe;

		fw = (y[npw] + y[npsw]) * 0.5;
		fe = (y[np] + y[nps]) * 0.5;

		yke[npw] = -3.0 * fw + 4.0 * yp[np] - fe;
	}

	// Derivatives relatively to csi at the center of the east face (only east boundary of the domain)
	i = nx - 1 - 1;
	for (j = 2; j <= ny - 1; j += 1)
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npw = np - 1;
		npsw = nps - 1;

		fw = (x[npw] + x[npsw]) * 0.5;
		fe = (x[np] + x[nps]) * 0.5;

		xke[np] = fw - 4.0 * xp[np] + 3.0 * fe;

		fw = (y[npw] + y[npsw]) * 0.5;
		fe = (y[np] + y[nps]) * 0.5;

		yke[np] = fw - 4.0 * yp[np] + 3.0 * fe;
	}

	// Derivatives relatively to eta at the center of the north face (only inner real faces)
	for (j = 2; j <= ny - 2; j += 1)
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)
		{
			np = nx * (j - 1) + i;
			npn = np + nx;

			xen[np] = xp[npn] - xp[np];
			yen[np] = yp[npn] - yp[np];
		}
	}

	// Derivatives relatively to eta at the center of the north face (only south boundary of the domain)
	j = 2;
	for (i = 1; i <= nx - 1 - 1; i += 1)
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npw = np - 1;
		npsw = nps - 1;

		fs = (x[nps] + x[npsw]) * 0.5;
		fn = (x[np] + x[npw]) * 0.5;

		xen[nps] = -3.0 * fs + 4.0 * xp[np] - fn;

		fs = (y[nps] + y[npsw]) * 0.5;
		fn = (y[np] + y[npw]) * 0.5;

		yen[nps] = -3.0 * fs + 4.0 * yp[np] - fn;
	}

	// Derivatives relatively to eta at the center of the north face (only north boundary of the domain)
	j = ny - 1;
	for (i = 1; i <= nx - 1 - 1; i += 1)
	{
		np = nx * (j - 1) + i;
		nps = np - nx;
		npw = np - 1;
		npsw = nps - 1;

		fs = (x[nps] + x[npsw]) * 0.5;
		fn = (x[np] + x[npw]) * 0.5;

		xen[np] = fs - 4.0 * xp[np] + 3.0 * fn;

		fs = (y[nps] + y[npsw]) * 0.5;
		fn = (y[np] + y[npw]) * 0.5;

		yen[np] = fs - 4.0 * yp[np] + 3.0 * fn;
	}

	// Beta and J at the center of the east face (all real faces)
	for (j = 2; j <= ny - 1; j += 1)
	{
		for (i = 0; i <= nx - 1 - 1; i += 1)
		{
			np = nx * (j - 1) + i;

			betae[np] = xke[np] * xe[np] + yke[np] * ye[np];

			Je[np] = 1.0 / (xke[np] * ye[np] - xe[np] * yke[np]);
		}
	}

	// Beta and J at center of the north face (all real faces)
	for (j = 1; j <= ny - 1; j += 1)
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)
		{
			np = nx * (j - 1) + i;

			betan[np] = xk[np] * xen[np] + yk[np] * yen[np];

			Jn[np] = 1.0 / (xk[np] * yen[np] - xen[np] * yk[np]);
		}
	}

	// Jacobian J at the center of all real volumes
	for (j = 2; j <= ny - 1; j += 1)
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)
		{
			np = nx * (j - 1) + i;
			nps = np - nx;
			npw = np - 1;
			npsw = nps - 1;

			fw = (x[npw] + x[npsw]) * 0.5;
			fe = (x[np] + x[nps]) * 0.5;
			fs = (x[nps] + x[npsw]) * 0.5;
			fn = (x[np] + x[npw]) * 0.5;

			xkp = fe - fw; // (x_csi)_P
			xep = fn - fs; // (x_eta)_P

			fw = (y[npw] + y[npsw]) * 0.5;
			fe = (y[np] + y[nps]) * 0.5;
			fs = (y[nps] + y[npsw]) * 0.5;
			fn = (y[np] + y[npw]) * 0.5;

			ykp = fe - fw; // (y_csi)_P
			yep = fn - fs; // (y_eta)_P

			Jp[np] = 1.0 / (xkp * yep - xep * ykp);
		}
	}

	for (int i = 0; i < nx*ny; i++)
	{
		Je[i] = std::fabs(Je[i]);
		Jn[i] = std::fabs(Jn[i]);
		Jp[i] = std::fabs(Jp[i]);
	}

	metrics result;
	result.metrics_coord.xe = xe;
	result.metrics_coord.ye = ye;
	result.metrics_coord.xen = xen;
	result.metrics_coord.yen = yen;
	result.metrics_coord.xk = xk;
	result.metrics_coord.yk = yk;
	result.metrics_coord.xke = xke;
	result.metrics_coord.yke = yke;
	result.jacobian.Jp = Jp;
	result.jacobian.Je = Je;
	result.jacobian.Jn = Jn;
	result.metric_tensor.alphae = alphae;
	result.metric_tensor.gamman = gamman;
	result.metric_tensor.betae = betae;
	result.metric_tensor.betan = betan;
	return result;
}

radius get_radius(int coord, int nx, int ny, double* y, double* yp, double* re, double* rn, double* rp)
{
	// Calculates the radius at center of faces and at the centroid of each real volume.
	
	//coord = Kind of coord. system ( 1=cylindrical, 0 = cartesian)

	// Cartesian coordinate system
	//re - Radius at the center of east face of volume P (m)
	//rn - Radius at the center of north face of volume P (m)
	//rp - Radius at the center of volume P (m)

	int i, j, np, npe, npw, npn, nps, npsw, npse, npnw, npne;

	for (int i = 0; i < nx*ny; i++)
	{
		re[i] = 1.0;
		rn[i] = 1.0;
		rp[i] = 1.0;
	}

	// Cylindrical coordinate system
	if (coord == 1)
	{
		// Radius at the center of the east face
		for (i = 0; i <= nx - 1 - 1; i += 1)
		{
			for (j = 2; j <= ny - 1; j += 1)
			{
				np = nx * (j - 1) + i;
				nps = np - nx;

				re[np] = (y[np] + y[nps]) * 0.5;
			}
		}

		// Radius at the center of the north face
		for (i = 1; i <= nx - 1 - 1; i += 1)
		{
			for (j = 1; j <= ny - 1; j += 1)
			{
				np = nx * (j - 1) + i;
				npw = np - 1;

				rn[np] = (y[np] + y[npw]) * 0.5;
			}
		}

		// Radius of the center of the volume
		rp = yp;

		// Radius  of the center of the volume (south fictitious)
		j = 1;
		for (i = 1; i <= nx - 1 - 1; i += 1)
		{
			np = nx * (j - 1) + i;
			npn = np + nx;

			rp[np] = rn[np] - rp[npn] + rn[np];
			//rp[np] = rp[npn] - rp[npn+nx] + rp[npn];
		}

		// Radius  of the center of the volume (north fictitious)
		j = ny;
		for (i = 1; i <= nx - 1 - 1; i += 1)
		{
			np = nx * (j - 1) + i;
			nps = np - nx;

			rp[np] = rn[nps] - rp[nps] + rn[nps];
			//rp[np] = rp[nps] - rp[nps-nx] + rp[nps];
		}

		// Radius  of the center of the volume (west fictitious)
		i = 1 - 1;
		for (j = 2; j <= ny - 1; j += 1)
		{
			np = nx * (j - 1) + i;
			npe = np + 1;

			rp[np] = re[np] - rp[npe] + re[np];
			//rp[np] = rp[npe] - rp[npe+1] + rp[npe];
		}

		// Radius  of the center of the volume (east fictitious)
		i = nx - 1;
		for (j = 2; j <= ny - 1; j += 1)
		{
			np = nx * (j - 1) + i;
			npw = np - 1;

			rp[np] = re[npw] - rp[npw] + re[npw];
			//rp[np] = rp[npw] - rp[npw-1] + rp[npw];
		}

		// Radius  of the center of the volume (corner fictitious: SW)
		i = 1 - 1;
		j = 1;

		np = nx * (j - 1) + i;
		npn = np + nx;
		npe = np + 1;
		npne = npn + 1;

		//rp[np] = rp[npn] - rp[npn+nx] + rp[npn];
		rp[np] = 4.0 * y[np] - rp[npe] - rp[npn] - rp[npne];

		// Radius  of the center of the volume (corner fictitious: SE)
		i = nx - 1;
		j = 1;

		np = nx * (j - 1) + i;
		npn = np + nx;
		npw = np - 1;
		npnw = npn - 1;

		//rp[np] = rp[npn] - rp[npn+nx] + rp[npn];
		rp[np] = 4.0 * y[npw] - rp[npw] - rp[npn] - rp[npnw];

		// Radius  of the center of the volume (corner fictitious: NW)
		i = 1 - 1;
		j = ny;

		np = nx * (j - 1) + i;
		nps = np - nx;
		npe = np + 1;
		npse = nps + 1;

		//rp[np] = rp[nps] - rp[nps-nx] + rp[nps];
		rp[np] = 4.0 * y[nps] - rp[npe] - rp[nps] - rp[npse];

		// Radius  of the center of the volume (corner fictitious: NE)
		i = nx - 1;
		j = ny;

		np = nx * (j - 1) + i;
		nps = np - nx;
		npw = np - 1;
		npsw = nps - 1;

		//rp[np] = rp[nps] - rp[nps-nx] + rp[nps];
		rp[np] = 4.0 * y[npsw] - rp[npw] - rp[nps] - rp[npsw];
	}

	radius result;
	result.rp = rp;
	result.re = re;
	result.rn = rn;
	return result;
}