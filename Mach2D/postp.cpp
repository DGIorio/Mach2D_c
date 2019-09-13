#include "stdafx.h"

// Contains subroutines for data post-processing.

std::string text_editor = "notepad ";				// Text editor name
std::string graph_viewer = "start 'a' ";			// Graph viewer name
std::string graph_generator = "./gnuplot/gnuplot ";	// Plotter name

void post_processing(int nx, int ny, int it, int sem_a, int sem_g, int w_g, int w_cam
	, std::string sim_id, std::chrono::duration<double, std::milli>, double RAM, double res, double lr, double rb, double* x, double* y, double* xe, double* ye, double* xk, double* yk, double* Jp, double Rg, double* cp
	, double* gcp, double* Uce, double* Vcn, double* de, double* dn, double* pl, double* bp, double* xp, double* yp, double* u, double* v, double* p, double* T, double* ro, double Cdfi)
{
	// Auxiliary variables
	int i;
	double* M;		// Mach number at the center of the volumes

	M = new double[nx*ny];

	// Writes grid parameters
	// TO DO

	// Calculates the boundary values and store them in the fictitious
	post_processing_boundaries(nx, ny, x, y, xp, yp, u, v, p, T);		 // InOutput: last six entries

	// Calculates the Mach number field
	for (i = 0; i < nx * ny; i++)
	{
		M[i] = sqrt((u[i] * u[i] + v[i] * v[i]) / (gcp[i] * Rg * T[i]));
	}

	// Calculates the specific mass field
	for (i = 0; i < nx * ny; i++)
	{
		ro[i] = p[i] / (Rg * T[i]);
	}

	// Write fields
	write_main_fields(nx, ny, sim_id, xp, yp, p, ro, T, u, v, M);
}

void write_main_fields(int nx, int ny, std::string sim_id, double* xp, double* yp, double* p, double* ro, double* T, double* u, double* v, double* M)
{
	// Prints the main fields.

	// Auxiliary variables
	int i, j, np;
	std::string str1;

	str1 = "./mach2d_output/mach2d_" + sim_id + "_" + "main_fields.dat";

	std::ofstream outputFile;
	CreateFolder("./mach2d_output");
	outputFile.open(str1, std::fstream::in | std::fstream::out | std::fstream::trunc);

	if (outputFile.is_open())
	{
		outputFile << "#" << std::setw(11) << "x" << " " << std::setw(11) << "y" << " " << std::setw(11) << "p" << " " << std::setw(11) << "ro" << " " << std::setw(11) << "T" << " " << std::setw(11) << "u" << " " << std::setw(11) << "v" << " " << std::setw(11) << "M" << std::endl;

		for (j = 1; j <= ny; j++)
			{
			for (i = 0; i < nx; i++)
				{
				np = (j - 1)*nx + i;

				outputFile << " " << std::setw(11) << std::scientific << std::setprecision(4)
					       << xp[np] << " " << std::setw(11) << yp[np] << " " << std::setw(11) << p[np] << " " << std::setw(11) << ro[np] << " "
						   << std::setw(11) << T[np] << " " << std::setw(11) << u[np] << " " << std::setw(11) << v[np] << " " << std::setw(11) << M[np] << std::endl;
				}
			outputFile << std::endl;
		}

		outputFile.close();
	}
	else std::cout << "Unable to open file to write main fields results" << std::endl;
}

void post_processing_boundaries(int nx, int ny, double* x, double* y, double* xp, double* yp, double* u, double* v, double* p, double* T) // InOutput: last six entries
{
	// Calculates the fields u, v, p and T over the boundaries and store them in the fictitious volumes.

	// Auxiliary variables
	int i, j, np, npn, npw, nps, npe, npsw, npse, npnw, npne;

	// west boundary (frontal symmetry line)
	i = 0;
	for (j = 2; j <= ny - 1; j++)
	{
		np = (j - 1)*nx + i;
		nps = np - nx;
		npe = np + 1;
		u[np] = (u[np] + u[npe]) * 0.5;
		v[np] = (v[np] + v[npe]) * 0.5;
		p[np] = (p[np] + p[npe]) * 0.5;
		T[np] = (T[np] + T[npe]) * 0.5;
		xp[np] = (x[nps] + x[np]) * 0.5;
		yp[np] = (y[nps] + y[np]) * 0.5;
	}

	// east boundary (exit)
	i = nx - 1;
	for (j = 2; j <= ny - 1; j++)
	{
		np = (j - 1)*nx + i;
		npw = np - 1;
		u[np] = (u[npw] + u[np]) * 0.5;
		v[np] = (v[npw] + v[np]) * 0.5;
		p[np] = (p[npw] + p[np]) * 0.5;
		T[np] = (T[npw] + T[np]) * 0.5;
		xp[np] = (x[npw - nx] + x[npw]) * 0.5;
		yp[np] = (y[npw - nx] + y[npw]) * 0.5;
	}

	// lower boundary (body wall)
	j = 1;
	for (i = 1; i <= nx - 1 - 1; i++)
	{
		np = (j - 1)*nx + i;
		npw = np - 1;
		npn = np + nx;
		u[np] = (u[npn] + u[np]) * 0.5;
		v[np] = (v[npn] + v[np]) * 0.5;
		p[np] = (p[npn] + p[np]) * 0.5;
		T[np] = (T[npn] + T[np]) * 0.5;
		xp[np] = (x[npw] + x[np]) * 0.5;
		yp[np] = (y[npw] + y[np]) * 0.5;
	}

	// SW corner
	j = 1;
	i = 0;
	np = (j - 1)*nx + i;
	npn = np + nx;
	npe = np + 1;
	npne = npn + 1;

	u[np] = (u[npe] + u[npn] + u[npne]) / 3.0;
	v[np] = (v[npe] + v[npn] + v[npne]) / 3.0;
	p[np] = (p[npe] + p[npn] + p[npne]) / 3.0;
	T[np] = (T[npe] + T[npn] + T[npne]) / 3.0;
	xp[np] = x[np];
	yp[np] = y[np];

	// SE corner
	j = 1;
	i = nx - 1;
	np = (j - 1)*nx + i;
	npn = np + nx;
	npw = np - 1;
	npnw = npn - 1;

	u[np] = (u[npw] + u[npn] + u[npnw]) / 3.0;
	v[np] = (v[npw] + v[npn] + v[npnw]) / 3.0;
	p[np] = (p[npw] + p[npn] + p[npnw]) / 3.0;
	T[np] = (T[npw] + T[npn] + T[npnw]) / 3.0;
	xp[np] = x[npw];
	yp[np] = y[npw];

	// upper boundary (far field)
	j = ny;
	for (i = 1; i <= nx - 1 - 1; i++)
	{
		np = (j - 1)*nx + i;
		nps = np - nx;
		u[np] = (u[np] + u[nps]) * 0.5;
		v[np] = (v[np] + v[nps]) * 0.5;
		p[np] = (p[np] + p[nps]) * 0.5;
		T[np] = (T[np] + T[nps]) * 0.5;
		xp[np] = (x[nps] + x[nps - 1]) * 0.5;
		yp[np] = (y[nps] + y[nps - 1]) * 0.5;
	}

	// NW corner
	j = ny;
	i = 0;
	np = (j - 1)*nx + i;
	nps = np - nx;
	npe = np + 1;
	npse = nps + 1;
	u[np] = (u[npe] + u[nps] + u[npse]) / 3.0;
	v[np] = (v[npe] + v[nps] + v[npse]) / 3.0;
	p[np] = (p[npe] + p[nps] + p[npse]) / 3.0;
	T[np] = (T[npe] + T[nps] + T[npse]) / 3.0;
	xp[np] = x[nps];
	yp[np] = y[nps];


	// NE corner
	j = ny;
	i = nx - 1;
	np = (j - 1)*nx + i;
	nps = np - nx;
	npw = np - 1;
	npsw = nps - 1;
	u[np] = (u[npw] + u[nps] + u[npsw]) / 3.0;
	v[np] = (v[npw] + v[nps] + v[npsw]) / 3.0;
	p[np] = (p[npw] + p[nps] + p[npsw]) / 3.0;
	T[np] = (T[npw] + T[nps] + T[npsw]) / 3.0;
	xp[np] = x[npsw];
	yp[np] = y[npsw];
}

double get_cdfi(int nx, int ny, int coord, double Rg, double PF, double TF, double UF, double rb, double* yk, double* rn, double* p, double Cdfi)
{
	// Calculates the pressure foredrag coefficient.

	// Auxiliary variables
	int i, j, np, npn;

	double ROF, QF, pn;

	// Far field density
	ROF = PF / (Rg * TF);

	// Far field dynamic pressure
	QF = ROF * UF * UF * 0.5;

	Cdfi = 0.0;

	j = 1;
	
	for (i = 1; i < nx - 1; i++)
	{
		np = nx * (j - 1) + i;

		npn = np + nx;

		pn = 0.5 * (p[np] + p[npn]);

		Cdfi = Cdfi + (pn - PF) * rn[np] * yk[np];
	}

	Cdfi = Cdfi * std::pow(2.0, coord) / (std::pow(rb, (coord + 1)) * QF);
	return Cdfi;
}