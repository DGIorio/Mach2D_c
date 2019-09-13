#include "stdafx.h"

// Goals:
//    Given the matrix of coefficients B of a linear system Bx = rhs 
//    1) Find the lower L and upper U triangular matrix of the factorization
//       B = LU - P. (P matrix is not important) 
//       B must have only five non-null diagonals.
//    2) Solve the linear system Bx = rhs iteratively using the LU factorization (MSI method)
//
//    Latest version: 28 Apr 2011
//
// Subroutines given by this module
//
//    1) fb2d5 
//       (ASSUMES THE COEFFICIENTS OF B MATRIX CORRESPONDING TO
//        BOUNDARY NEIGHBOORS OF FICTITIOUS VOLUMES ARE ZERO)
//
//    2) lu2d5 
//       (SETS ALPHA=0 AT FICTITIOUS VOLUMES INDEPENDENTLY OF ALPHA ELSEWHERE)
//       (the parameter alpha must be set inside this subroutine)

double* fb2d5(double* b, double* dl, double* du, int nj, int ni, int nij, double* var, double gamma, int nitm, double* rhs, double* r, double* w, double* z)
{
	// fb2d5 ASSUMES THE COEFFICIENTS OF B MATRIX CORRESPONDING TO
	// BOUNDARY NEIGHBOORS OF FICTITIOUS VOLUMES ARE ZERO
	//
	// fb2d5 solves the linear system b.var = rhs iteratively
	//
	// b     - matrix  of coefficients of the original system (5 diagonals)
	// dl    - lower triangular matrix 
	// du    - upper triangular matrix
	// nj    - number of volumes in the j direction
	// ni    - number of volumes in the i direction
	// nij   - nij=ni*nj
	// var   - vector of unkowns
	// gamma - L2( r_n ) < L2( r_0 ) * gamma
	// L2    - L2 norm
	// r_n   - residual vector: r_n = b var_n - rhs
	// nitm  - maximum number of iteractions
	// rhs   - vector with the source term of the linear system

	// Auxiliary variables
	int i;
	double sum_r;
	int nit;				 // iteraction counter
	int is, iw, ip, ie, in;
	double gammarp;

	//r - residual vector
	//w - vector used to solve linear system Lw = -r
	//z - vector used to solve linear system Uz =  w

	
	// Calculation of the residual vector r = A.var-rhs
	//--------------------------------------------------------------------------------
	ip = 0;
	ie = ip + 1;
	in = ip + ni;
	r[ip] = -rhs[ip] + b[ip * 5 + 2] * var[ip] + b[ip * 5 + 3] * var[ie]
		+ b[ip * 5 + 4] * var[in];
	//--------------------------------------------------------------------------------
	for (ip = 1; ip <= ni - 1 - 1; ip += 1)
	{
		iw = ip - 1;
		ie = ip + 1;
		in = ip + ni;
		r[ip] = -rhs[ip] + b[ip * 5 + 1] * var[iw] + b[ip * 5 + 2] * var[ip] + b[ip * 5 + 3] * var[ie]
			+ b[ip * 5 + 4] * var[in];
	}
	//--------------------------------------------------------------------------------
	ip = ni - 1;
	iw = ip - 1;
	in = ip + ni;
	r[ip] = -rhs[ip] + b[ip * 5 + 1] * var[iw] + b[ip * 5 + 2] * var[ip]
		+ b[ip * 5 + 4] * var[in];
	//--------------------------------------------------------------------------------
	for (ip = ni; ip <= nij - ni - 1; ip += 1)
	{
		iw = ip - 1;
		ie = ip + 1;
		is = ip - ni;
		in = ip + ni;
		switch ((ip + 1) % ni)
		{
		case 0:
			r[ip] = -rhs[ip] + b[ip * 5 + 0] * var[is]
				+ b[ip * 5 + 1] * var[iw] + b[ip * 5 + 2] * var[ip]
				+ b[ip * 5 + 4] * var[in];
			break;
		case 1:
			r[ip] = -rhs[ip] + b[ip * 5 + 0] * var[is]
				+ b[ip * 5 + 2] * var[ip] + b[ip * 5 + 3] * var[ie]
				+ b[ip * 5 + 4] * var[in];
			break;
		default:
			r[ip] = -rhs[ip] + b[ip * 5 + 0] * var[is]
				+ b[ip * 5 + 1] * var[iw] + b[ip * 5 + 2] * var[ip] + b[ip * 5 + 3] * var[ie]
				+ b[ip * 5 + 4] * var[in];
			break;
		}
	}
	//--------------------------------------------------------------------------------
	ip = nij - ni;
	ie = ip + 1;
	is = ip - ni;
	r[ip] = -rhs[ip] + b[ip * 5 + 0] * var[is]
		+ b[ip * 5 + 2] * var[ip] + b[ip * 5 + 3] * var[ie];
	//--------------------------------------------------------------------------------
	for (ip = nij - ni + 1; ip <= nij - 1 - 1; ip += 1)
	{
		iw = ip - 1;
		ie = ip + 1;
		is = ip - ni;
		r[ip] = -rhs[ip] + b[ip * 5 + 0] * var[is]
			+ b[ip * 5 + 1] * var[iw] + b[ip * 5 + 2] * var[ip] + b[ip * 5 + 3] * var[ie];
	}
	//--------------------------------------------------------------------------------
	ip = nij - 1;
	iw = ip - 1;
	is = ip - ni;
	r[ip] = -rhs[ip] + b[ip * 5 + 0] * var[is]
		+ b[ip * 5 + 1] * var[iw] + b[ip * 5 + 2] * var[ip];
	//--------------------------------------------------------------------------------
	sum_r = 0.0;
	//#pragma omp parallel for reduction(+:sum_r)
	for (i = 0; i < nij; i++)
	{
		sum_r += r[i] * r[i];
	}
	//gammarp = gamma * sqrt(std::inner_product(r, r + nij, r, 0.0));
	gammarp = gamma * sqrt(sum_r);
	//--------------------------------------------------------------------------------
	nit = 0;
	for (;;)
	{
		nit = nit + 1;

		// Solving linear system Lw = -r
		//------------------------------
		ip = 0;
		w[ip] = -r[ip] / dl[ip * 4 + 0];
		//------------------------------
		for (ip = 1; ip <= ni - 1 - 1; ip += 1)
		{
			w[ip] = -(r[ip] + dl[ip * 4 + 1] * w[ip - 1]) / dl[ip * 4 + 0];
		}
		//--------------------------------------------------
		ip = ni - 1;
		w[ip] = -(r[ip] + dl[ip * 4 + 1] * w[ip - 1] + dl[ip * 4 + 2] * w[ip - ni + 1]) / dl[ip * 4 + 0];
		//---------------------------------------------------------------------
		for (ip = ni; ip <= nij - 1; ip += 1)
		{
			w[ip] = -(r[ip] + dl[ip * 4 + 1] * w[ip - 1] + dl[ip * 4 + 2] * w[ip - ni + 1] + dl[ip * 4 + 3] * w[ip - ni]) / dl[ip * 4 + 0];
		}
		//--------------------------------------------------------------------------------------------

		// Solving linear system Uz = w
		//-----------------------------
		ip = nij - 1;
		z[ip] = w[ip];
		//------------
		for (ip = nij - 2; ip >= nij - ni + 2 - 1; ip -= 1)
		{
			z[ip] = w[ip] - du[ip * 3 + 0] * z[ip + 1];
		}
		//----------------------------------
		ip = nij - ni;
		z[ip] = w[ip] - du[ip * 3 + 0] * z[ip + 1] - du[ip * 3 + 1] * z[ip + ni - 1];
		//-----------------------------------------------------
		//for (ip = nij - ni - 1 + 1; ip-- > 0;)
		for (ip = nij - ni - 1; ip >= 0; ip -= 1)
		{
			z[ip] = w[ip] - du[ip * 3 + 0] * z[ip + 1] - du[ip * 3 + 1] * z[ip + ni - 1] - du[ip * 3 + 2] * z[ip + ni];
		}
		//----------------------------------------------------------------------------

		// Updating solution
		//------------------
		for (ip = 0; ip <= nij - 1; ip += 1)
		{
			var[ip] = z[ip] + var[ip];
		}

		// Calculation of the residual vector r = A.var-rhs
		//--------------------------------------------------------------------------------
		ip = 0;
		ie = ip + 1;
		in = ip + ni;
		r[ip] = -rhs[ip] + b[ip * 5 + 2] * var[ip] + b[ip * 5 + 3] * var[ie]
			+ b[ip * 5 + 4] * var[in];
		//--------------------------------------------------------------------------------
		for (ip = 1; ip <= ni - 1 - 1; ip += 1)
		{
			iw = ip - 1;
			ie = ip + 1;
			in = ip + ni;
			r[ip] = -rhs[ip] + b[ip * 5 + 1] * var[iw] + b[ip * 5 + 2] * var[ip] + b[ip * 5 + 3] * var[ie]
				+ b[ip * 5 + 4] * var[in];
		}
		//--------------------------------------------------------------------------------
		ip = ni - 1;
		iw = ip - 1;
		in = ip + ni;
		r[ip] = -rhs[ip] + b[ip * 5 + 1] * var[iw] + b[ip * 5 + 2] * var[ip]
			+ b[ip * 5 + 4] * var[in];
		//--------------------------------------------------------------------------------
		for (ip = ni; ip <= nij - ni - 1; ip += 1)
		{
			iw = ip - 1;
			ie = ip + 1;
			is = ip - ni;
			in = ip + ni;
			switch ((ip + 1) % ni)
			{
			case 0:
				r[ip] = -rhs[ip] + b[ip * 5 + 0] * var[is]
					+ b[ip * 5 + 1] * var[iw] + b[ip * 5 + 2] * var[ip]
					+ b[ip * 5 + 4] * var[in];
				break;
			case 1:
				r[ip] = -rhs[ip] + b[ip * 5 + 0] * var[is]
					+ b[ip * 5 + 2] * var[ip] + b[ip * 5 + 3] * var[ie]
					+ b[ip * 5 + 4] * var[in];
				break;
			default:
				r[ip] = -rhs[ip] + b[ip * 5 + 0] * var[is]
					+ b[ip * 5 + 1] * var[iw] + b[ip * 5 + 2] * var[ip] + b[ip * 5 + 3] * var[ie]
					+ b[ip * 5 + 4] * var[in];
				break;
			}
		}
		//--------------------------------------------------------------------------------
		ip = nij - ni;
		ie = ip + 1;
		is = ip - ni;
		r[ip] = -rhs[ip] + b[ip * 5 + 0] * var[is]
			+ b[ip * 5 + 2] * var[ip] + b[ip * 5 + 3] * var[ie];
		//--------------------------------------------------------------------------------
		for (ip = nij - ni + 1; ip <= nij - 1 - 1; ip += 1)
		{
			iw = ip - 1;
			ie = ip + 1;
			is = ip - ni;
			r[ip] = -rhs[ip] + b[ip * 5 + 0] * var[is]
				+ b[ip * 5 + 1] * var[iw] + b[ip * 5 + 2] * var[ip] + b[ip * 5 + 3] * var[ie];
		}
		//--------------------------------------------------------------------------------
		ip = nij - 1;
		iw = ip - 1;
		is = ip - ni;
		r[ip] = -rhs[ip] + b[ip * 5 + 0] * var[is]
			+ b[ip * 5 + 1] * var[iw] + b[ip * 5 + 2] * var[ip];
		//--------------------------------------------------------------------------------
		sum_r = 0.0;
		//#pragma omp parallel for reduction(+:sum_r)
		for (i = 0; i < nij; i++)
		{
			sum_r += r[i] * r[i];
		}
		//if (sqrt(std::inner_product(r, r+nij, r, 0.0)) < gammarp)
		if (sqrt(sum_r) < gammarp)
		{
			break;
		}
		if (nit == nitm)
		{
			break;
		}
	}
	return var;
}

diagonal_matrix lu2d5(double* b, int nij, int ni, int nj, double* dl, double* du)
{
	// lu2d5 SETS ALPHA=0 AT FICTITIOUS VOLUMES INDEPENDENTLY OF ALPHA OTHERWISE
	// lu2d5 calculates the entries of the L and U matrix in the decomposition LU = B + P
	// where B is the matrix of the coefficients of the linear system to be solved and P is a
	// matrix whose entries are not necessary.
	//
	// b     - matrix  of coefficients of the original system (has 5 diagonals)
	// dl    - lower triangular matrix 
	// du    - upper triangular matrix
	// nj    - number of volumes in the j direction
	// ni    - number of volumes in the i direction
	// nij   - nij=ni*nj

	// Auxiliary variables
	int is, ise, iw, ip;
	double alpha, phi1, phi4, alpha0;
	int mod_ip_ni;

	alpha0 = 0.5;
	alpha = alpha0;

	//------------------------------------------------------------------------------------------------------
	ip = 0;
	dl[ip * 4 + 0] = b[ip * 5 + 2];
	du[ip * 3 + 0] = b[ip * 5 + 3] / dl[ip * 4 + 0];
	du[ip * 3 + 2] = b[ip * 5 + 4] / dl[ip * 4 + 0];
	//------------------------------------------------------------------------------------------------------
	for (ip = 1; ip <= ni - 1; ip += 1)
	{
		iw = ip - 1;
		dl[ip * 4 + 1] = b[ip * 5 + 1];
		dl[ip * 4 + 0] = b[ip * 5 + 2] - dl[ip * 4 + 1] * du[iw * 3 + 0];
		du[ip * 3 + 0] = b[ip * 5 + 3] / dl[ip * 4 + 0];
		du[ip * 3 + 1] = -dl[ip * 4 + 1] * du[iw * 3 + 2] / dl[ip * 4 + 0];
		du[ip * 3 + 2] = b[ip * 5 + 4] / dl[ip * 4 + 0];
	}
	//------------------------------------------------------------------------------------------------------
	for (ip = ni; ip <= nij - ni - 1; ip += 1)
	{
		iw = ip - 1;
		is = ip - ni;
		ise = is + 1;
		mod_ip_ni = (ip + 1) % ni;
		if (mod_ip_ni == 0 || mod_ip_ni == 1)
		{
			alpha = 0.0;
		}
		else
		{
			alpha = alpha0;
		}
		dl[ip * 4 + 3] = b[ip * 5 + 0] / (1.0 - alpha * du[is * 3 + 0] * du[ise * 3 + 0]);
		dl[ip * 4 + 2] = -dl[ip * 4 + 3] * du[is * 3 + 0];
		dl[ip * 4 + 1] = (b[ip * 5 + 1] - dl[ip * 4 + 3] * du[is * 3 + 1]) / (1.0 + 2.0 * alpha * du[iw * 3 + 1]);
		phi1 = dl[ip * 4 + 2] * du[ise * 3 + 0];
		phi4 = dl[ip * 4 + 1] * du[iw * 3 + 1];
		dl[ip * 4 + 0] = b[ip * 5 + 2] - dl[ip * 4 + 1] * du[iw * 3 + 0] - dl[ip * 4 + 2] * du[ise * 3 + 1] - dl[ip * 4 + 3] * du[is * 3 + 2] + 2.0*alpha*(phi1 + phi4);
		du[ip * 3 + 0] = (b[ip * 5 + 3] - dl[ip * 4 + 2] * du[ise * 3 + 2] - 2.0*alpha*phi1) / dl[ip * 4 + 0];
		du[ip * 3 + 1] = -dl[ip * 4 + 1] * du[iw * 3 + 2] / dl[ip * 4 + 0];
		du[ip * 3 + 2] = (b[ip * 5 + 4] - alpha * phi4) / dl[ip * 4 + 0];
	}
	//------------------------------------------------------------------------------------------------------
	ip = nij - ni;
	iw = ip - 1;
	is = ip - ni;
	ise = is + 1;
	dl[ip * 4 + 3] = b[ip * 5 + 0];
	dl[ip * 4 + 2] = -dl[ip * 4 + 3] * du[is * 3 + 0];
	dl[ip * 4 + 1] = (b[ip * 5 + 1] - dl[ip * 4 + 3] * du[is * 3 + 1]);
	dl[ip * 4 + 0] = b[ip * 5 + 2] - dl[ip * 4 + 1] * du[iw * 3 + 0] - dl[ip * 4 + 2] * du[ise * 3 + 1] - dl[ip * 4 + 3] * du[is * 3 + 2];
	du[ip * 3 + 0] = (b[ip * 5 + 3] - dl[ip * 4 + 2] * du[ise * 3 + 2]) / dl[ip * 4 + 0];
	du[ip * 3 + 1] = -dl[ip * 4 + 1] * du[iw * 3 + 2] / dl[ip * 4 + 0];
	//------------------------------------------------------------------------------------------------------
	for (ip = nij - ni + 1; ip <= nij - 1 - 1; ip += 1)
	{
		iw = ip - 1;
		is = ip - ni;
		ise = is + 1;
		dl[ip * 4 + 3] = b[ip * 5 + 0];
		dl[ip * 4 + 2] = -dl[ip * 4 + 3] * du[is * 3 + 0];
		dl[ip * 4 + 1] = b[ip * 5 + 1] - dl[ip * 4 + 3] * du[is * 3 + 1];
		dl[ip * 4 + 0] = b[ip * 5 + 2] - dl[ip * 4 + 1] * du[iw * 3 + 0] - dl[ip * 4 + 2] * du[ise * 3 + 1] - dl[ip * 4 + 3] * du[is * 3 + 2];
		du[ip * 3 + 0] = (b[ip * 5 + 3] - dl[ip * 4 + 2] * du[ise * 3 + 2]) / dl[ip * 4 + 0];
	}
	//------------------------------------------------------------------------------------------------------
	ip = nij - 1;
	iw = ip - 1;
	is = ip - ni;
	ise = is + 1;
	dl[ip * 4 + 3] = b[ip * 5 + 0];
	dl[ip * 4 + 2] = -dl[ip * 4 + 3] * du[is * 3 + 0];
	dl[ip * 4 + 1] = b[ip * 5 + 1] - dl[ip * 4 + 3] * du[is * 3 + 1];
	dl[ip * 4 + 0] = b[ip * 5 + 2] - dl[ip * 4 + 1] * du[iw * 3 + 0] - dl[ip * 4 + 2] * du[ise * 3 + 1] - dl[ip * 4 + 3] * du[is * 3 + 2];
	//------------------------------------------------------------------------------------------------------
	diagonal_matrix result;
	result.lower = dl;
	result.upper = du;
	return result;
}