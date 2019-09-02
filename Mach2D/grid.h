#pragma once

// inout x, y
coordinates set_grid(int kg, int nx, int ny, double a1, double* x, double* y);
// inout x, y
coordinates get_uniform_grid(int nx, int ny, double* x, double* y);
// inout x, y
coordinates get_power_grid(double a1, int nx, int ny, double* x, double* y);
// inout x, y
coordinates get_pg_grid(double a1, int nx, int ny, double* x, double* y);

double get_gp_ratio(int n, double r);

coordinates get_real_centroids_xy(int opt, int nx, int ny, double* x, double* y, double* xp, double* yp);

metrics get_metrics(int nx, int ny, double* x, double* y, double* xp, double* yp,
					double* xe, double* ye, double* xen, double* yen, double* xk, double* yk, double* xke, double* yke,
					double* Jp, double* Je, double* Jn, double* alphae, double* gamman, double* betae, double* betan);

radius get_radius(int coord, int nx, int ny, double* y, double* yp, double* re, double* rn, double* rp);