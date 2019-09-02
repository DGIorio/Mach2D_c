#pragma once
double* get_newton_interpolation_coefficients(int n, double* xn, double* fn, double* c);

double get_newton_interpolation_f0(int n, double x, double* xn, double* c);

double get_newton_interpolation_f1(int n, double x, double* xn, double* c);

interpolation get_ni_f0_vec_d1(int n, double* xn, double* fn, double xc);

interpolation get_ni_f0_vec_d2(int n, double* xn, double* fn, double xc);

int searcher(int n, double* x, double xc);
