#pragma once

#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <iterator>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <windows.h>
#include <stdio.h>
#include <math.h>

#include "typedefs.h"
#include "main.h"
#include "data.h"
#include "grid.h"
#include "coef.h"
#include "user.h"
#include "solvers.h"
#include "msi2d5-v03.h"
#include "msi2d9-v03.h"
#include "newton_interpolation.h"
#include "postp.h"

//Extern variables
extern int nx;
extern int ny;
extern int nxy;
extern int kg;
extern int kcm;
extern int coord;
extern double la;
extern double lb;
extern double a1;
extern double akn;
extern double aks;

extern double* x;
extern double* y;
extern double* xp;
extern double* yp;
extern double* xe;
extern double* ye;
extern double* xen;
extern double* yen;
extern double* xk;
extern double* yk;
extern double* xke;
extern double* yke;
extern double* Jp;
extern double* Je;
extern double* Jn;
extern double* alphae;
extern double* gamman;
extern double* betae;
extern double* betan;
extern double* re;
extern double* rn;
extern double* rp;

extern int it;
extern int iti;
extern int it1;
extern int it2;
extern int itb1;
extern int itb2;
extern int itm;
extern int imax;
extern int itmax;
extern int itimax;
extern int nitm_u;
extern int nitm_p;
extern int wlf;
extern int sem_a;
extern int sem_g;
extern int w_g;
extern int w_cam;
extern int it_stop;
extern int nit_res;
extern int res_check;
extern double tol_u;
extern double tol_p;
extern double beta;
extern double beta1;
extern double beta2;
extern double dt;
extern double dt1;
extern double dt2;
extern double res_u;
extern double res_v;
extern double res_T;
extern double res_p;
extern double res;
extern double res_mean;
extern double tol_res;
extern double norm;
extern double tcpu1;
extern double tcpu2;
extern double tcpu;
extern double RAM;

extern double* bu;
extern double* bv;
extern double* bt;
extern double* bp;
extern double* pl;
extern double* due;
extern double* dve;
extern double* dun;
extern double* dvn;
extern double* de;
extern double* dn;

extern double* au;
extern double* av;
extern double* at;
extern double* ap;
extern double* dl9;
extern double* du9;
extern double* dl5;
extern double* du5;

extern double lr;
extern double rb;

extern int const tfid;
extern int const lid;
extern int const rid;

extern std::string date_s;
extern std::string time_s;
extern std::string sim_id;
extern std::string input_file_parameters;

extern double GF;
extern double Rg;
extern double PF;
extern double TF;
extern double MF;
extern double UF;
extern double Cdfi;

extern double* p;
extern double* pa;
extern double* T;
extern double* Ta;
extern double* ro;
extern double* roa;
extern double* roe;
extern double* ron;
extern double* g;
extern double* u;
extern double* ua;
extern double* v;
extern double* va;
extern double* ue;
extern double* un;
extern double* ve;
extern double* vn;
extern double* Uce;
extern double* Vcn;
extern double* cp;
extern double* gcp;
extern double* cup;
extern double* cvp;
extern double* uea;
extern double* vea;
extern double* una;
extern double* vna;
extern double* Ucea;
extern double* Vcna;
extern double* Vcw;

extern double* r;
extern double* w;
extern double* z;