#pragma once

// inout var
double* fb2d9(double** b, double** dl, double** du, int nj, int ni, int nij, double* var, double gamma, int nitm, double* rhs, double* r, double* w, double* z);

diagonal_matrix lu2d9(double** b, int nij, int ni, int nj, double** dl, double** du);