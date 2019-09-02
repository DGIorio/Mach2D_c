#pragma once

double* get_density_at_nodes(int nx, int ny, double Rg, double* p, double* T, double* ro);

specific_mass get_density_at_faces(int nx, int ny, double beta, double* ro, double* Uce, double* Vcn, double* roe, double* ron);

double** get_u_coefficients(int nx, int ny, double dt, double* rp, double* re, double* rn, double* Jp, double* roe, double* ron, double* roa, double* Uce, double* Vcn, double** a);

source_correction get_u_source(int nx, int ny, double beta, double dt, double* rp, double* re, double* rn, double* ye, double* yk, double* Jp, double* roe, double* ron, double* roa, double* p, double* Uce, double* Vcn, double* ua, double* u, double* cup, double* b);

double** get_v_coefficients(int nx, int ny, double dt, double* rp, double* re, double* rn, double* Jp, double* roe, double* ron, double* roa, double* Uce, double* Vcn, double** a);

source_correction get_v_source(int nx, int ny, double beta, double dt, double* rp, double* re, double* rn, double* xe, double* xk, double* Jp, double* roe, double* ron, double* roa, double* p, double* Uce, double* Vcn, double* va, double* v, double* cvp, double* b);

coeffs_source get_T_coefficients_and_source(int nx, int ny, double beta, double dt, double* rp, double* re, double* rn, double* xe, double* ye, double* xk, double* yk, double* Jp, double* roe, double* ron, double* roa, double* p, double* pa, double* cp, double* Uce, double* Vcn, double* u, double* v, double* T, double* Ta, double** a, double* b);

double** get_p_coefficients(int nx, int ny, double dt, double* rp, double* re, double* rn, double* Jp, double* Uce, double* Vcn, double* roe, double* ron, double* g, double* de, double* dn, double** a);

double* get_p_source(int nx, int ny, double dt, double* rp, double* re, double* rn, double* Jp, double* roe, double* ron, double* roem, double* ronm, double* ro, double* roa, double* Ucem, double* Uce, double* Vcnm, double* Vcn, double* b);

velocity_face get_velocities_at_internal_faces(int nx, int ny, double dt, double* rp, double* re, double* rn, double* xe, double* ye, double* xk, double* yk,
												double* xen, double* yen, double* xke, double* yke, double* Jp, double* cup, double* cvp,
												double** au, double** av, double* roa, double* p, double* u, double* v, double* uea, double* vea, double* una, double* vna,
												double* ue, double* ve, double* un, double* vn, double* Uce, double* Vcn);

simplec_coeffs get_internal_simplec_coefficients(int nx, int ny, double* re, double* rn, double* xe, double* ye, double* xk, double* yk, double** au, double** av,
												double* due, double* dve, double* dun, double* dvn, double* de, double* dn);

// inout ro, p
pressure_specificmass get_pressure_density_correction_with_pl(int nx, int ny, double* pl, double* g, double* ro, double* p);

// inout u, v
velocity get_u_v_at_real_nodes_with_pl(int nx, int ny, double* xe, double* ye, double* xk,
								double* yk, double* rp, double* pl, double** au, double** av,
								double* u, double* v);

// inout ue, ve, un, vn, Uce, Vcn
velocity_face get_velocity_correction_at_faces_with_pl(int nx, int ny, double* due, double* dve, double* dun, double* dvn, double* de, double* dn, double* pl,
														double* ue, double* ve, double* un, double* vn, double* Uce, double* Vcn);

