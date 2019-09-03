#include "stdafx.h"

// Contains user-INDEPENDENT subroutines related to the calculation
// of coefficients and source of the linear systems.

double* get_density_at_nodes(int nx, int ny, double Rg, double* p, double* T, double* ro)
{
	// Calculates the specific mass at nodes using p, T and the state equation.
	
	//double* ro;	// Specific mass (absolute density) at center of volume P (kg/m3);
	//ro = new double[nx*ny];

	// Auxiliary variables
	int i;

	for (i = 0; i <= nx * ny - 1; i += 1)
	{
		ro[i] = p[i] / (Rg*T[i]);
	}
	return ro;
}

specific_mass get_density_at_faces(int nx, int ny, double beta, double* ro, double* Uce, double* Vcn, double* roe, double* ron)
{
	// Interpolates the specific mass at center of volume faces using the specific mass at nodes
	// based on the scheme UDS with deferred correction to CDS.

	//double* roe;	// Specific mass at east face of volume P (kg/m3);
	//double* ron;	// Specific mass at north face of volume P (kg/m3);

	// Auxiliary variables
	int i, j, np, npe, npn;
	double ae, an;				// Coefficients of UDS scheme

	//roe = new double[nx*ny];
	//ron = new double[nx*ny];

	// Density at east face
	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 0; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			npe = np + 1;
			
			ae = 0.5*sign(Uce[np]);
			roe[np] = (0.5 + ae) * ro[np] + (0.5 - ae) * ro[npe] + beta * ae * (ro[npe] - ro[np]);
		}
	}

	// Density at north face
	for (j = 1; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			npn = np + nx;

			an = 0.5*sign(Vcn[np]);
			ron[np] = (0.5 + an) * ro[np] + (0.5 - an) * ro[npn] + beta * an * (ro[npn] - ro[np]);
		}
	}
	specific_mass result;
	result.roe = roe;
	result.ron = ron;
	return result;
}

double** get_u_coefficients(int nx, int ny, double dt, double* rp, double* re, double* rn, double* Jp, double* roe, double* ron, double* roa, double* Uce, double* Vcn, double** a)
{
	// Calculates the coefficients of the linear system for u for the real volumes.

	//double** a;   // Coefficients of the linear system for u (kg/s);

	// Auxiliary variables
	int i, j, np, nps, npw;
	double fmw, fme, fms, fmn, mpa;
	double ae, aw, an, as;

	// Initializating
	//a = new double*[nx*ny];
	//for (i = 0; i < nx*ny; i++)
	//{
	//	a[i] = new double[9];
	//}

	// Contribution to the coefficients due to the advection term
	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			nps = np - nx;
			npw = np - 1;

			fme = roe[np] * re[np] * Uce[np];
			fmw = roe[npw] * re[npw] * Uce[npw];
			fmn = ron[np] * rn[np] * Vcn[np];
			fms = ron[nps] * rn[nps] * Vcn[nps];

			mpa = roa[np] * rp[np] / Jp[np];

			as = 0.5*sign(Vcn[nps]);
			an = 0.5*sign(Vcn[np]);
			aw = 0.5*sign(Uce[npw]);
			ae = 0.5*sign(Uce[np]);

			a[np][0] = 0.0; // SW
			a[np][2] = 0.0; // SE
			a[np][6] = 0.0; // NW
			a[np][8] = 0.0; // NE

			a[np][1] = -fms * (0.5 + as); // S
			a[np][3] = -fmw * (0.5 + aw); // W
			a[np][5] = fme * (0.5 - ae);  // E
			a[np][7] = fmn * (0.5 - an);  // N

			a[np][4] = mpa / dt - (a[np][7] + a[np][1] + a[np][3] + a[np][5]);
		}
	}
	return a;
}

source_correction get_u_source(int nx, int ny, double beta, double dt, double* rp, double* re, double* rn, double* ye, double* yk, double* Jp, double* roe, double* ron, double* roa, double* p, double* Uce, double* Vcn, double* ua, double* u, double* cup, double* b)
{
	// Calculates the source of the linear system for u for the real volumes.

	//double* cup;		// Term of deferred correction for u (N);
	//double* b;			// Source of the linear system for u (N);

	// Auxiliary variables
	int i, j, np, nps, npn, npw, npe;
	double fmw, fme, fmn, fms, mpa;
	double as, an, aw, ae;

	// Initializating
	//cup = new double[nx*ny];
	//b = new double[nx*ny];

	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			nps = np - nx;
			npn = np + nx;
			npw = np - 1;
			npe = np + 1;

			fme = roe[np] * re[np] * Uce[np];
			fmw = roe[npw] * re[npw] * Uce[npw];
			fmn = ron[np] * rn[np] * Vcn[np];
			fms = ron[nps] * rn[nps] * Vcn[nps];

			mpa = roa[np] * rp[np] / Jp[np];

			as = 0.5*sign(Vcn[nps]);
			an = 0.5*sign(Vcn[np]);
			aw = 0.5*sign(Uce[npw]);
			ae = 0.5*sign(Uce[np]);

			// Contribution to b due to advection (UDS)
			b[np] = mpa * ua[np] / dt;

			// Contribution to b due to advection (Deferred correction)
			cup[np] = -beta *
				(fme * ae * (u[npe] - u[np])
					- fmw * aw * (u[np] - u[npw])
					+ fmn * an * (u[npn] - u[np])
					- fms * as * (u[np] - u[nps])
					);

			b[np] = b[np] + cup[np];

			// Contribution to b due to pressure term
			b[np] = b[np] + 0.5 * rp[np] * (
				  yk[np] * (p[npn] + p[np])
				- yk[nps] * (p[nps] + p[np])
				- ye[np] * (p[npe] + p[np])
				+ ye[npw] * (p[npw] + p[np])
				);
		}
	}
	source_correction result;
	result.cuvp = cup;
	result.b = b;
	return result;
}

double** get_v_coefficients(int nx, int ny, double dt, double* rp, double* re, double* rn, double* Jp, double* roe, double* ron, double* roa, double* Uce, double* Vcn, double** a)
{
	// Calculates the coefficients of the linear system for v  for the real volumes.

	//double** a;		// Coefficients of the linear system for v (kg/s);

	// Auxiliary variables
	int i, j, np, nps, npw;
	double fmw, fme, fms, fmn, mpa;
	double ae, aw, an, as;

	// Initializating
	//a = new double*[nx*ny];
	//for (i = 0; i < nx*ny; i++)
	//{
	//	a[i] = new double[9];
	//}

	// Contribution to the coefficients due to the advection term
	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			nps = np - nx;
			npw = np - 1;

			fme = roe[np] * re[np] * Uce[np];
			fmw = roe[npw] * re[npw] * Uce[npw];
			fmn = ron[np] * rn[np] * Vcn[np];
			fms = ron[nps] * rn[nps] * Vcn[nps];

			mpa = roa[np] * rp[np] / Jp[np];

			as = 0.5*sign(Vcn[nps]);
			an = 0.5*sign(Vcn[np]);
			aw = 0.5*sign(Uce[npw]);
			ae = 0.5*sign(Uce[np]);

			a[np][0] = 0.0; // SW
			a[np][2] = 0.0; // SE
			a[np][6] = 0.0; // NW
			a[np][8] = 0.0; // NE

			a[np][1] = -fms * (0.5 + as); // S
			a[np][3] = -fmw * (0.5 + aw); // W
			a[np][5] = fme * (0.5 - ae);  // E
			a[np][7] = fmn * (0.5 - an);  // N

			a[np][4] = mpa / dt - (a[np][7] + a[np][1] + a[np][3] + a[np][5]);
		}
	}
	return a;
}

source_correction get_v_source(int nx, int ny, double beta, double dt, double* rp, double* re, double* rn, double* xe, double* xk, double* Jp, double* roe, double* ron, double* roa, double* p, double* Uce, double* Vcn, double* va, double* v, double* cvp, double* b)
{
	// Calculates the source of the linear system for v for the real volumes.

	//double* cvp;		// Term of deferred correction for v (N);
	//double* b;			// Source of the linear system for v (N);

	// Auxiliary variables
	int i, j, np, nps, npn, npw, npe;
	double fmw, fme, fmn, fms, mpa;
	double as, an, aw, ae;

	// Initializating
	//cvp = new double[nx*ny];
	//b = new double[nx*ny];

	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			nps = np - nx;
			npn = np + nx;
			npw = np - 1;
			npe = np + 1;

			fme = roe[np] * re[np] * Uce[np];
			fmw = roe[npw] * re[npw] * Uce[npw];
			fmn = ron[np] * rn[np] * Vcn[np];
			fms = ron[nps] * rn[nps] * Vcn[nps];

			mpa = roa[np] * rp[np] / Jp[np];

			as = 0.5*sign(Vcn[nps]);
			an = 0.5*sign(Vcn[np]);
			aw = 0.5*sign(Uce[npw]);
			ae = 0.5*sign(Uce[np]);

			// Contribution to b due to advection (UDS)
			b[np] = mpa * va[np] / dt;

			// Contribution to b due to advection (Deferred correction)
			cvp[np] = -beta *
				(fme * ae * (v[npe] - v[np])
					- fmw * aw * (v[np] - v[npw])
					+ fmn * an * (v[npn] - v[np])
					- fms * as * (v[np] - v[nps])
					);

			b[np] = b[np] + cvp[np];

			// Contribution to b due to pressure term
			b[np] = b[np] + 0.5 * rp[np] * (
				  xe[np] * (p[np] + p[npe])
				- xe[npw] * (p[np] + p[npw])
				- xk[np] * (p[np] + p[npn])
				+ xk[nps] * (p[np] + p[nps])
				);
		}
	}
	source_correction result;
	result.cuvp = cvp;
	result.b = b;
	return result;
}

coeffs_source get_T_coefficients_and_source(int nx, int ny, double beta, double dt, double* rp, double* re, double* rn, double* xe, double* ye, double* xk, double* yk, double* Jp, double* roe, double* ron, double* roa, double* p, double* pa, double* cp, double* Uce, double* Vcn, double* u, double* v, double* T, double* Ta, double** a, double* b)
{
	// Calculates the coefficients and source of the linear system for T for the real volumes.

	//double** a;		// Coefficients of the linear system for T (J/(K.s))
	//double* b;		// Source of the linear system for T (J/s)

	// Auxiliary variables
	int i, j;
	int np, nps, npn, npw, npe;
	double fmw, fme, fms, fmn, mpa, pup, pvp;
	double as, an, aw, ae;

	// Initializating
	//b = new double[nx*ny];
	//a = new double*[nx*ny];
	//for (i = 0; i < nx*ny; i++)
	//{
	//	a[i] = new double[9];
	//}

	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			nps = np - nx;
			npn = np + nx;
			npw = np - 1;
			npe = np + 1;

			fme = roe[np] * re[np] * Uce[np];
			fmw = roe[npw] * re[npw] * Uce[npw];
			fmn = ron[np] * rn[np] * Vcn[np];
			fms = ron[nps] * rn[nps] * Vcn[nps];

			mpa = roa[np] * rp[np] / Jp[np];

			as = 0.5*sign(Vcn[nps]);
			an = 0.5*sign(Vcn[np]);
			aw = 0.5*sign(Uce[npw]);
			ae = 0.5*sign(Uce[np]);

			// Contribution to the COEFFICIENTS due to advection (UDS)
			a[np][0] = 0.0; // SW
			a[np][2] = 0.0; // SE
			a[np][6] = 0.0; // NW
			a[np][8] = 0.0; // NE

			a[np][1] = -fms * (0.5 + as) * cp[np]; // S
			a[np][3] = -fmw * (0.5 + aw) * cp[np]; // W
			a[np][5] = fme * (0.5 - ae) * cp[np];  // E
			a[np][7] = fmn * (0.5 - an) * cp[np];  // N

			a[np][4] = mpa * cp[np] / dt - (a[np][7] + a[np][1] + a[np][3] + a[np][5]);

			// Contribution to the SOURCE due to advection (UDS)
			b[np] = mpa * Ta[np] * cp[np] / dt;

			// Contribution to the SOURCE due to advection (Deferred correction)
			b[np] = b[np] - beta * cp[np] *
				(fme * ae * (T[npe] - T[np])
					- fmw * aw * (T[np] - T[npw])
					+ fmn * an * (T[npn] - T[np])
					- fms * as * (T[np] - T[nps])
					);

			// Contribution to the SOURCE due to pressure term
			pup = 0.5 *  rp[np] * (
				  yk[np] * (p[npn] + p[np])
				- yk[nps] * (p[nps] + p[np])
				- ye[np] * (p[npe] + p[np])
				+ ye[npw] * (p[npw] + p[np])
				);

			pvp = 0.5 *   rp[np] * (
				  xe[np] * (p[np] + p[npe])
				- xe[npw] * (p[np] + p[npw])
				- xk[np] * (p[np] + p[npn])
				+ xk[nps] * (p[np] + p[nps])
				);

			b[np] = b[np] + rp[np] / Jp[np] * (p[np] - pa[np]) / dt
				- u[np] * pup - v[np] * pvp;
		}
	}
	coeffs_source result;
	result.a = a;
	result.b = b;
	return result;
}

double** get_p_coefficients(int nx, int ny, double dt, double* rp, double* re, double* rn, double* Jp, double* Uce, double* Vcn, double* roe, double* ron, double* g, double* de, double* dn, double** a)
{
	// Calculates the coefficients of the linear system for pressure correction
	// g, ro, Uce and Vcn used in this subroutine are those calculated in the previous iteraction
	// de and dn must be calculated with the coef. of the linear sytem for u and v from which
	// u* and v* are obtained.

	//double** a;			// Coefficients of the linear system for pl (m.s)

	// Auxiliary variables
	int i, j, np, nps, npn, npw, npe;
	double as, an, aw, ae;
	double roem, rowm, ronm, rosm;

	// Initializating
	//a = new double*[nx*ny];
	//for (i = 0; i < nx*ny; i++)
	//{
	//	a[i] = new double[5];
	//}

	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			nps = np - nx;
			npn = np + nx;
			npw = np - 1;
			npe = np + 1;

			as = 0.5*sign(Vcn[nps]);
			an = 0.5*sign(Vcn[np]);
			aw = 0.5*sign(Uce[npw]);
			ae = 0.5*sign(Uce[np]);

			roem = roe[np];   // Density on east face
			rowm = roe[npw];  // Density on west face
			ronm = ron[np];   // Density on north face
			rosm = ron[nps];  // Density on south face

			// South
			a[np][0] = -rn[nps] * ((0.5 + as) * Vcn[nps] * g[nps] + rosm * dn[nps]);

			// West
			a[np][1] = -re[npw] * ((0.5 + aw) * Uce[npw] * g[npw] + rowm * de[npw]);

			// Center
			a[np][2] = g[np] * (rp[np] / (Jp[np] * dt)
				+ (0.5 + ae) * re[np] * Uce[np]
				- (0.5 - aw) * re[npw] * Uce[npw]
				+ (0.5 + an) * rn[np] * Vcn[np]
				- (0.5 - as) * rn[nps] * Vcn[nps])
				+ roem * re[np] * de[np]
				+ rowm * re[npw] * de[npw]
				+ ronm * rn[np] * dn[np]
				+ rosm * rn[nps] * dn[nps];

			// East
			a[np][3] = re[np] * ((0.5 - ae) * Uce[np] * g[npe] - roem * de[np]);

			// North
			a[np][4] = rn[np] * ((0.5 - an) * Vcn[np] * g[npn] - ronm * dn[np]);
		}
	}
	return a;
}

double* get_p_source(int nx, int ny, double dt, double* rp, double* re, double* rn, double* Jp, double* roe, double* ron, double* roem, double* ronm, double* ro, double* roa, double* Ucem, double* Uce, double* Vcnm, double* Vcn, double* b)
{
	// Calculates the source of the linear system of the pressure deviation.
	// ro, Uce and Vcn are the incorrect ones (obtained with p*)
	// rom, Ucem and Vcnm are those of the previous iteraction

	//double* b;			// Source of the linear system for pl (kg/s)

	// Auxiliary variables
	int i, j, np, nps, npn, npw, npe;
	double as, an, aw, ae; // UDS coefficients

	// Initializating
	//b = new double[nx*ny];

	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			nps = np - nx;
			npn = np + nx;
			npw = np - 1;
			npe = np + 1;

			as = 0.5*sign(Vcnm[nps]);
			an = 0.5*sign(Vcnm[np]);
			aw = 0.5*sign(Ucem[npw]);
			ae = 0.5*sign(Ucem[np]);

			b[np] = -rp[np] * (ro[np] - roa[np]) / (Jp[np] * dt)
				- roem[np] * re[np] * Uce[np]
				+ roem[npw] * re[npw] * Uce[npw]
				- ronm[np] * rn[np] * Vcn[np]
				+ ronm[nps] * rn[nps] * Vcn[nps]
				+ (roem[np] - roe[np]) * re[np] * Ucem[np]
				- (roem[npw] - roe[npw]) * re[npw] * Ucem[npw]
				+ (ronm[np] - ron[np]) * rn[np] * Vcnm[np]
				- (ronm[nps] - ron[nps]) * rn[nps] * Vcnm[nps];
		}
	}
	return b;
}

velocity_face get_velocities_at_internal_faces(int nx, int ny, double dt, double* rp, double* re, double* rn, double* xe, double* ye, double* xk, double* yk,
	double* xen, double* yen, double* xke, double* yke, double* Jp, double* cup, double* cvp,
	double** au, double** av, double* roa, double* p, double* u, double* v, double* uea, double* vea, double* una, double* vna,
	double* ue, double* ve, double* un, double* vn, double* Uce, double* Vcn)
{
	// Calculates the velocities at the interfaces between real volumes.
	// The interpolation scheme applied here is the one appropriate for the SIMPLEC method.

	//double* ue(nx*ny);	// Cartesian velocity u at center of east face of volume P (m/s);
	//double* ve(nx*ny);	// Cartesian velocity v at center of east face of volume P (m/s);
	//double* un(nx*ny);	// Cartesian velocity u at center of north face of volume P (m/s);
	//double* vn(nx*ny);	// Cartesian velocity v at center of north face of volume P (m/s);
	//double* Uce(nx*ny);	// Contravariant velocity U at east face of volume P (m2/s);
	//double* Vcn(nx*ny);	// Contravariant velocity V at north face of volume P (m2/s);

	// Auxiliary variables
	int i, j, np, npsw, npse, nps, npw, npe, npnw, npn, npne;
	int npee, npsee, npnee, npnn, npnnw, npnne;
	double mpa, mea, mna, aux;
	double sumup, sumun, sumvp, sumvn, sumue, sumve;
	double pue, pve, pun, pvn;

	// Calculation of Uce at the inner faces of the real domain
	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 2 - 1; i += 1)			// already fixed loop - verify
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

			npee = npe + 1;
			npsee = npse + 1;
			npnee = npne + 1;

			mpa = roa[np] * rp[np] / Jp[np];
			mea = roa[npe] * rp[npe] / Jp[npe];

			sumup = au[np][0] * u[npsw] + au[np][1] * u[nps] + au[np][2] * u[npse] + au[np][3] * u[npw]
				+ au[np][5] * u[npe] + au[np][6] * u[npnw] + au[np][7] * u[npn] + au[np][8] * u[npne];

			sumue = au[npe][0] * u[nps] + au[npe][1] * u[npse] + au[npe][2] * u[npsee] + au[npe][3] * u[np]
				+ au[npe][5] * u[npee] + au[npe][6] * u[npn] + au[npe][7] * u[npne] + au[npe][8] * u[npnee];

			sumvp = av[np][0] * v[npsw] + av[np][1] * v[nps] + av[np][2] * v[npse] + av[np][3] * v[npw]
				+ av[np][5] * v[npe] + av[np][6] * v[npnw] + av[np][7] * v[npn] + av[np][8] * v[npne];

			sumve = av[npe][0] * v[nps] + av[npe][1] * v[npse] + av[npe][2] * v[npsee] + av[npe][3] * v[np]
				+ av[npe][5] * v[npee] + av[npe][6] * v[npn] + av[npe][7] * v[npne] + av[npe][8] * v[npnee];

			aux = (p[npn] + p[npne] - p[nps] - p[npse]) / 4.0;

			pue = re[np] * (yke[np] * aux + ye[np] * (p[np] - p[npe]));

			pve = re[np] * (-xke[np] * aux + xe[np] * (p[npe] - p[np]));

			ue[np] = ((mpa + mea) * uea[np] / dt + cup[np] + cup[npe]
				- sumup - sumue + 2.0 * pue) / (au[np][4] + au[npe][4]);

			ve[np] = ((mpa + mea) * vea[np] / dt + cvp[np] + cvp[npe]
				- sumvp - sumve + 2.0 * pve) / (av[np][4] + av[npe][4]);

			Uce[np] = ue[np] * ye[np] - ve[np] * xe[np];
		}
	}

	// Calculation of Vcn at the inner faces of the real domain
	for (j = 2; j <= ny - 2; j += 1)					// already fixed loop - verify
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

			npnn = npn + nx;
			npnnw = npnn - 1;
			npnne = npnn + 1;

			mpa = roa[np] * rp[np] / Jp[np];
			mna = roa[npn] * rp[npn] / Jp[npn];

			sumup = au[np][0] * u[npsw] + au[np][1] * u[nps] + au[np][2] * u[npse] + au[np][3] * u[npw]
				+ au[np][5] * u[npe] + au[np][6] * u[npnw] + au[np][7] * u[npn] + au[np][8] * u[npne];

			sumun = au[npn][0] * u[npw] + au[npn][1] * u[np] + au[npn][2] * u[npe] + au[npn][3] * u[npnw]
				+ au[npn][5] * u[npne] + au[npn][6] * u[npnnw] + au[npn][7] * u[npnn] + au[npn][8] * u[npnne];

			sumvp = av[np][0] * v[npsw] + av[np][1] * v[nps] + av[np][2] * v[npse] + av[np][3] * v[npw]
				+ av[np][5] * v[npe] + av[np][6] * v[npnw] + av[np][7] * v[npn] + av[np][8] * v[npne];

			sumvn = av[npn][0] * v[npw] + av[npn][1] * v[np] + av[npn][2] * v[npe] + av[npn][3] * v[npnw]
				+ av[npn][5] * v[npne] + av[npn][6] * v[npnnw] + av[npn][7] * v[npnn] + av[npn][8] * v[npnne];

			aux = (p[npne] + p[npe] - p[npnw] - p[npw]) / 4.0;

			pun = rn[np] * (yk[np] * (p[npn] - p[np]) - yen[np] * aux);

			pvn = rn[np] * (xk[np] * (p[np] - p[npn]) + xen[np] * aux);

			un[np] = ((mpa + mna) * una[np] / dt + cup[np] + cup[npn]
				- sumup - sumun + 2.0 * pun) / (au[np][4] + au[npn][4]);

			vn[np] = ((mpa + mna) * vna[np] / dt + cvp[np] + cvp[npn]
				- sumvp - sumvn + 2.0 * pvn) / (av[np][4] + av[npn][4]);

			Vcn[np] = vn[np] * xk[np] - un[np] * yk[np];
		}
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

simplec_coeffs get_internal_simplec_coefficients(int nx, int ny, double* re, double* rn, double* xe, double* ye, double* xk, double* yk, double** au, double** av,
	double* due, double* dve, double* dun, double* dvn, double* de, double* dn)
{
	// Calculates SIMPLEC coefficients at the interface between real volumes.

	//double* due;	  // SIMPLEC coefficients for ue (m2.s/kg);
	//double* dve;	  // SIMPLEC coefficients for ve (m2.s/kg);
	//double* dun;	  // SIMPLEC coefficients for un (m2.s/kg);
	//double* dvn;	  // SIMPLEC coefficients for vn (m2.s/kg);
	//double* de;	  // SIMPLEC coefficients for Uce (m3.s/kg);
	//double* dn;	  // SIMPLEC coefficients for Vcn (m3.s/kg);

	// Auxiliary variables;
	int i, j, k, np, npe, npn;
	double sum_au_np;
	double sum_au_npe;
	double sum_au_npn;
	double sum_av_np;
	double sum_av_npe;
	double sum_av_npn;

	// Initializating
	//due = new double[nx*ny];
	//dve = new double[nx*ny];
	//dun = new double[nx*ny];
	//dvn = new double[nx*ny];
	//de  = new double[nx*ny];
	//dn  = new double[nx*ny];
	
	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 2 - 1; i += 1)			// already fixed loop - verify
		{
			sum_au_np = 0.0;
			sum_au_npe = 0.0;
			sum_av_np = 0.0;
			sum_av_npe = 0.0;

			np = nx * (j - 1) + i;
			npe = np + 1;

			//best
			//sum_au_np = std::accumulate(au[np], au[np] + 9, 0.0, std::plus<double>());
			//sum_au_npe = std::accumulate(au[npe], au[npe] + 9, 0.0, std::plus<double>());
			//sum_av_np = std::accumulate(av[np], av[np] + 9, 0.0, std::plus<double>());
			//sum_av_npe = std::accumulate(av[npe], av[npe] + 9, 0.0, std::plus<double>());
			//worst
			//std::for_each(au[np], au[np] + 9, [&](double val) { sum_au_np += val; });
			//std::for_each(au[npe], au[npe] + 9, [&](double val) { sum_au_npe += val; });
			//std::for_each(av[np], av[np] + 9, [&](double val) { sum_av_np += val; });
			//std::for_each(av[npe], av[npe] + 9, [&](double val) { sum_av_npe += val; });
			for (k = 0; k < 9; k++)
			{
				sum_au_np += au[np][k];
				sum_au_npe += au[npe][k];
				sum_av_np += av[np][k];
				sum_av_npe += av[npe][k];
			}

			due[np] = 2.0 * re[np] * ye[np] / (sum_au_np + sum_au_npe);
			dve[np] = 2.0 * re[np] * xe[np] / (sum_av_np + sum_av_npe);

			de[np] = ye[np] * due[np] + xe[np] * dve[np];
		}
	}

	for (j = 2; j <= ny - 2; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			sum_au_np = 0.0;
			sum_au_npn = 0.0;
			sum_av_np = 0.0;
			sum_av_npn = 0.0;

			np = nx * (j - 1) + i;
			npn = np + nx;

			//sum_au_np = std::accumulate(au[np], au[np] + 9, 0.0, std::plus<double>());
			//sum_au_npn = std::accumulate(au[npn], au[npn] + 9, 0.0, std::plus<double>());
			//sum_av_np = std::accumulate(av[np], av[np] + 9, 0.0, std::plus<double>());
			//sum_av_npn = std::accumulate(av[npn], av[npn] + 9, 0.0, std::plus<double>());
			//std::for_each(au[np], au[np] + 9, [&](double val) { sum_au_np += val; });
			//std::for_each(au[npn], au[npn] + 9, [&](double val) { sum_au_npn += val; });
			//std::for_each(av[np], av[np] + 9, [&](double val) { sum_av_np += val; });
			//std::for_each(av[npn], av[npn] + 9, [&](double val) { sum_av_npn += val; });
			for (k = 0; k < 9; k++)
			{
				sum_au_np += au[np][k];
				sum_au_npn += au[npn][k];
				sum_av_np += av[np][k];
				sum_av_npn += av[npn][k];
			}

			dun[np] = 2.0 * rn[np] * yk[np] / (sum_au_np + sum_au_npn);
			dvn[np] = 2.0 * rn[np] * xk[np] / (sum_av_np + sum_av_npn);

			dn[np] = xk[np] * dvn[np] + yk[np] * dun[np];
		}
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

pressure_specificmass get_pressure_density_correction_with_pl(int nx, int ny, double* pl, double* g, double* ro, double* p)
{
	// Corrects pressure and specific mass at all nodes (real+fictitious) using the
	// pressure deviation pl.

	// Auxiliary variables;
	int i;

	//std::transform(p.begin(), p.end(), pl.begin(), p.begin(), std::plus<double>());
	for (i = 0; i <= nx * ny - 1; i += 1)					// already fixed loop - verify
	{
		p[i] = p[i] + pl[i];
		ro[i] = ro[i] + pl[i]*g[i];
	}
	pressure_specificmass result;
	result.p = p;
	result.ro = ro;
	return result;
}

velocity get_u_v_at_real_nodes_with_pl(int nx, int ny, double* xe, double* ye, double* xk,
	double* yk, double* rp, double* pl, double** au, double** av,
	double* u, double* v)
{
	// Corrects u and v at all real nodes with the pressure deviation pl.
	// Input (u*,v*) ->  Output (u,v).
	// Auxiliary variables
	int i, j, k, np, nps, npn, npw, npe;
	double plup, plvp;
	double sum_au_np;
	double sum_av_np;
	
	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			sum_au_np = 0.0;
			sum_av_np = 0.0;

			np = nx * (j - 1) + i;
			nps = np - nx;
			npn = np + nx;
			npw = np - 1;
			npe = np + 1;

			plup = 0.5 *  rp[np] * (
				  yk[np] * (pl[npn] + pl[np])
				- yk[nps] * (pl[nps] + pl[np])
				- ye[np] * (pl[npe] + pl[np])
				+ ye[npw] * (pl[npw] + pl[np])
				);

			plvp = 0.5 *   rp[np] * (
				  xe[np] * (pl[np] + pl[npe])
				- xe[npw] * (pl[np] + pl[npw])
				- xk[np] * (pl[np] + pl[npn])
				+ xk[nps] * (pl[np] + pl[nps])
				);

			//sum_au_np = std::accumulate(au[np], au[np] + 9, 0.0, std::plus<double>());
			//sum_av_np = std::accumulate(av[np], av[np] + 9, 0.0, std::plus<double>());
			//std::for_each(au[np], au[np] + 9, [&](double val) { sum_au_np += val; });
			//std::for_each(av[np], av[np] + 9, [&](double val) { sum_av_np += val; });
			for (k = 0; k < 9; k++)
			{
				sum_au_np += au[np][k];
				sum_av_np += av[np][k];
			}

			u[np] = u[np] + plup / sum_au_np;
			v[np] = v[np] + plvp / sum_av_np;
		}
	}
	velocity result;
	result.u = u;
	result.v = v;
	return result;
}

velocity_face get_velocity_correction_at_faces_with_pl(int nx, int ny, double* due, double* dve, double* dun, double* dvn, double* de, double* dn, double* pl,
	double* ue, double* ve, double* un, double* vn, double* Uce, double* Vcn)
{
	// Corrects the velocities ue, ve, un, vn, Uce, Vcn at all faces using the pressure
	// deviation pl. Input (ue*, ve*, un*, vn*, Uce*, Vcn*) -> Output (ue, ve, un, vn, Uce, Vcn)
	// Auxiliary variables
	int i, j, np, npn, npe;
	
	for (j = 2; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 0; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			npe = np + 1;

			Uce[np] = Uce[np] + de[np] * (pl[np] - pl[npe]);

			ue[np] = ue[np] + due[np] * (pl[np] - pl[npe]);

			ve[np] = ve[np] + dve[np] * (pl[npe] - pl[np]);
		}
	}
	
	for (j = 1; j <= ny - 1; j += 1)					// already fixed loop - verify
	{
		for (i = 1; i <= nx - 1 - 1; i += 1)			// already fixed loop - verify
		{
			np = nx * (j - 1) + i;
			npn = np + nx;

			Vcn[np] = Vcn[np] + dn[np] * (pl[np] - pl[npn]);

			un[np] = un[np] + dun[np] * (pl[npn] - pl[np]);

			vn[np] = vn[np] + dvn[np] * (pl[np] - pl[npn]);
		}
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