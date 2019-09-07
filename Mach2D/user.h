#pragma once

coordinates get_grid_boundary_power_law(int nx, int ny, double akn, double aks, double la, double lb, double lr, double rb, double* x, double* y);

double f_body_shape(double x, double lr, double rb);

double* set_cp(int nx, int ny, double* cp);

double* set_gamma(int nx, int ny, double Rg, double* cp, double* gcp);

// talvez mudar a, b para inout
coeffs_source set_bcu(int nx, int ny, double UF, double* xk, double* yk, double* alphae, double* betae, double* u, double* v, double* Uce, double* Vcw, double* a, double* b);

coeffs_source set_bcv(int nx, int ny, double* xk, double* yk, double* u, double* v, double* Uce, double* Vcw, double* a, double* b);

coeffs_source set_bcT(int nx, int ny, double TF, double* Uce, double* Vcw, double* T, double* alphae, double* betae, double* betan, double* gamman, double* a, double* b);

coeffs_source set_bcp(int nx, int ny, double* alphae, double* betae, double* betan, double* gamman, double* Uce, double* Vcw, double* p, double* a, double* b);

double* get_Vcw(int nx, int ny, double* xe, double* xk, double* ue, double* ve, double* Vcw);

// inout u, v
velocity get_u_v_extrapolation_to_fictitious(int nx, int ny, double UF, double* xk, double* yk, double* alphae, double* betae, double* Uce, double* Vcw, double* u, double* v);

// inout T
double* get_T_extrapolation_to_fictitious(int nx, int ny, double TF, double* alphae, double* betae, double* betan, double* gamman, double* Uce, double* Vcw, double* T);

// inout p
double* get_p_extrapolation_to_fictitious(int nx, int ny, double PF, double* alphae, double* betae, double* betan, double* gamman, double* Uce, double* Vcw, double* p);

// inout f
double* get_extrapolation_to_corners(int nx, int ny, double* f);

// inout due, dve, dun, dvn, de, dn
simplec_coeffs get_boundary_simplec_coefficients(int nx, int ny, double* xe, double* ye, double* due, double* dve, double* dun, double* dvn, double* de, double* dn);

velocity_face get_velocities_at_boundary_faces(int nx, int ny, double UF, double* xe, double* ye, double* xk, double* yk, double* u, double* v,
												double* ue, double* ve, double* un, double* vn, double* Uce, double* Vcn);

initial_conditions get_initial_conditions(int nx, int ny, double PF, double TF, double Rg, double GF, double MF, double UF, double* xe, double* ye, double* xk, double* yk, double* alphae, double* betae, double* betan, double* gamman,
	double* p, double* pl, double* T, double* ro, double* roe, double* ron, double* u, double* v, double* ue, double* ve, double* un, double* vn, double* Uce, double* Vcn, double* Vcw);

// Inout vars (may add an &, read about references types in C++)
// double* &T = 1