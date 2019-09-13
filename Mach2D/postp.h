#pragma once

void post_processing(int nx, int ny, int it, int sem_a, int sem_g, int w_g, int w_cam
	, std::string sim_id, std::chrono::duration<double, std::milli>, double RAM, double res, double lr, double rb, double* x, double* y, double* xe, double* ye, double* xk, double* yk, double* Jp, double Rg, double* cp
	, double* gcp, double* Uce, double* Vcn, double* de, double* dn, double* pl, double* bp, double* xp, double* yp, double* u, double* v, double* p, double* T, double* ro, double Cdfi);

void write_main_fields(int nx, int ny, std::string sim_id, double* xp, double* yp, double* p, double* ro, double* T, double* u, double* v, double* M);

void post_processing_boundaries(int nx, int ny, double* x, double* y, double* xp, double* yp, double* u, double* v, double* p, double* T);

double get_cdfi(int nx, int ny, int coord, double Rg, double PF, double TF, double UF, double rb, double* yk, double* rn, double* p, double Cdfi);