#pragma once

//typedef std::vector<std::string> string_vec_1D;
//typedef std::vector<int> int_vec_1D;
//typedef std::vector<double> double*;
//typedef std::vector<std::vector<double>> double**;

struct coordinates
{ 
	double* x;
	double* y;
};

struct radius
{
	double* rp;
	double* re;
	double* rn;
};

struct specific_mass
{
	double* roe;
	double* ron;
};

struct pressure_specificmass
{
	double* p;
	double* ro;
};

struct source_correction
{
	double* cuvp;
	double* b;
};

struct metrics_coord
{
	double* xe;
	double* ye;
	double* xen;
	double* yen;
	double* xk;
	double* yk;
	double* xke;
	double* yke;
};

struct jacobian
{
	double* Jp;
	double* Je;
	double* Jn;
};

struct metric_tensor
{
	double* alphae;
	double* gamman;
	double* betae;
	double* betan;
};

struct metrics
{
	metrics_coord metrics_coord;
	jacobian jacobian;
	metric_tensor metric_tensor;
};

struct velocity
{
	double* u;
	double* v;
};

struct coeffs_source
{
	double* a;
	double* b;
};

struct simplec_coeffs
{
	double* due;
	double* dve;
	double* dun;
	double* dvn;
	double* de;
	double* dn;
};

struct velocity_face
{
	double* ue;
	double* ve;
	double* un;
	double* vn;
	double* Uce;
	double* Vcn;
};

struct initial_conditions
{
	double UF;
	double* p;
	double* pl;
	double* T;
	double* ro;
	double* roe;
	double* ron;
	velocity velocity;
	velocity_face velocity_face;
	double* Vcw;
};

struct diagonal_matrix
{
	double* lower;
	double* upper;
};

struct interpolation
{
	double fc;
	int status;
};

template <class T>
inline int
sign(T v) {
	return (v >= T(0)) - (v < T(0));
}