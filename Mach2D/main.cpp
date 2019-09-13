#include "stdafx.h"

int main()
{
	// Reads parameters, allocates variables and sets vectors to zero
	allocate_initialize_variables();

	int i;
	coordinates xy;
	coordinates xpyp;
	source_correction cup_bu;
	source_correction cvp_bv;
	coeffs_source au_bu;
	coeffs_source av_bv;
	coeffs_source ap_bp;
	coeffs_source at_bt;
	simplec_coeffs simplec_coefficients;
	diagonal_matrix d5;
	diagonal_matrix d9;
	velocity u_v;
	velocity_face velocities_face;
	pressure_specificmass p_ro;
	specific_mass roe_ron;


	// DELETAR REMOVER
	int n_p, ii, jj, kk;
	std::ofstream DGIoutputFile;
	CreateFolder("./mach2d_output");
	DGIoutputFile.open("./mach2d_output/dgi_output.txt", std::fstream::in | std::fstream::out | std::fstream::trunc);
	std::streamsize ss = std::cout.precision();


	// Generates the north and south boundaries of the grid based on a power-law model.					OK
	xy = get_grid_boundary_power_law(nx, ny, akn, aks, la, lb, lr, rb, x, y);
	x = xy.x;
	y = xy.y;

	// Generates the grid according to kg option														OK
	xy = set_grid(kg, nx, ny, a1, x, y);
	x = xy.x;
	y = xy.y;

	// Calculates the centroids of each real volume														OK
	xpyp =  get_real_centroids_xy(kcm, nx, ny, x, y, xp, yp);
	xp = xpyp.x;
	yp = xpyp.y;

	// Calculates the components of the metric tensor and the Jacobian									OK
	metrics grid_metrics = get_metrics(nx, ny, x, y, xp, yp,
										xe, ye, xen, yen, xk, yk, xke, yke, Jp,  Je, Jn, alphae, gamman, betae, betan);
	xe = grid_metrics.metrics_coord.xe;
	ye = grid_metrics.metrics_coord.ye;
	xen = grid_metrics.metrics_coord.xen;
	yen = grid_metrics.metrics_coord.yen;
	xk = grid_metrics.metrics_coord.xk;
	yk = grid_metrics.metrics_coord.yk;
	xke = grid_metrics.metrics_coord.xke;
	yke = grid_metrics.metrics_coord.yke;
	Jp = grid_metrics.jacobian.Jp;
	Je = grid_metrics.jacobian.Je;
	Jn = grid_metrics.jacobian.Jn;
	alphae = grid_metrics.metric_tensor.alphae;
	gamman = grid_metrics.metric_tensor.gamman;
	betae = grid_metrics.metric_tensor.betae;
	betan = grid_metrics.metric_tensor.betan;

	// All the radii are equal 1 if planar flow is choose. If axisymmetric flow, radii are calculated	OK
	radius grid_radius = get_radius(coord, nx, ny, y, yp, re, rn, rp);
	rp = grid_radius.rp;
	re = grid_radius.re;
	rn = grid_radius.rn;

	// Guess the initial conditions																		OK
	initial_conditions init_cond = get_initial_conditions(nx, ny, PF, TF, Rg, GF, MF, UF, xe, ye, xk, yk, alphae, betae, betan, gamman,
															p, pl, T, ro, roe, ron, u, v, ue, ve, un, vn, Uce, Vcn, Vcw);
	UF = init_cond.UF;		// can not comment
	p = init_cond.p;
	pl = init_cond.pl;
	T = init_cond.T;
	ro = init_cond.ro;
	roe = init_cond.roe;
	ron = init_cond.ron;
	u = init_cond.velocity.u;
	v = init_cond.velocity.v;
	ue = init_cond.velocity_face.ue;
	ve = init_cond.velocity_face.ve;
	un = init_cond.velocity_face.un;
	vn = init_cond.velocity_face.vn;
	Uce = init_cond.velocity_face.Uce;
	Vcn = init_cond.velocity_face.Vcn;
	Vcw = init_cond.Vcw;

	// Calculations of the thermal properties
	cp = set_cp(nx, ny, cp);

	// Beginning time evolution cycle

	// Openning file of residuals
	std::ofstream outputFile;
	CreateFolder("./mach2d_output");
	outputFile.open("./mach2d_output/mach2d_" + sim_id + "_residual.dat", std::fstream::in | std::fstream::out | std::fstream::trunc);

	// Initialising CPU time
	auto tcpu1 = std::chrono::steady_clock::now();

	for (it = 1; it <= itmax; it += 1)
	{
		// Changing the time step
		if (it <= it1) dt = dt1;
		if (it > it1 && it < it2) dt = dt1 + (dt2 - dt1)*(it - it1) / (it2 - it1);
		if (it >= it2) dt = dt2;

		// Changing the beta value
		if (it <= itb1) beta = beta1;
		if (it > itb1 && it < itb2) beta = beta1 + (beta2 - beta1)*(it - itb1) / (itb2 - itb1);
		if (it >= itb2) beta = beta2;

		// Fields update
		roa = ro;
		Ta = T;
		pa = p;
		ua = u;
		va = v;
		uea = ue;
		vea = ve;
		una = un;
		vna = vn;

		// Starting the iterative cycle for the instant t
		for (iti = 1; iti <= itimax; iti++)
		{
			// Updating variables
			Ucea = Uce;
			Vcna = Vcn;

			// Calculations of the thermal properties
			cp = set_cp(nx, ny, cp);

			// Coefficients of the linear system for u (real volumes)							au OK
			au = get_u_coefficients(nx, ny, dt, rp, re, rn, Jp
									, roe, ron, roa, Uce, Vcn
									, au);

			// Source of the linear system for u (real volumes)									cup, bu OK
			cup_bu = get_u_source(nx, ny, beta, dt, rp, re, rn, ye, yk
								, Jp, roe, ron, roa, p, Uce, Vcn, ua, u, cup, bu);
			cup = cup_bu.cuvp;
			bu = cup_bu.b;

			// Coefficients and source of the linear system of u for the fictitious volumes		au, bu OK
			au_bu = set_bcu(nx, ny, UF, xk, yk, alphae, betae, u, v, Uce, Vcw, au, bu);
			au = au_bu.a;
			bu = au_bu.b;

			// Coefficients of the linear system for v (real volumes)							av OK
			av = get_v_coefficients(nx, ny, dt, rp, re, rn, Jp
									, roe, ron, roa, Uce, Vcn
									, av);

			// Source of the linear system for v (real volumes)									cvp OK								bv - five differences (because xe and xk)
			cvp_bv = get_v_source(nx, ny, beta, dt, rp, re, rn, xe, xk
								, Jp, roe, ron, roa, p, Uce, Vcn, va, v, cvp, bv);
			cvp = cvp_bv.cuvp;
			bv = cvp_bv.b;

			// Coefficients and source of the linear system of v (fictitious volumes)			av OK								bv - same differences
			av_bv = set_bcv(nx, ny, xk, yk, u, v, Uce, Vcw, av, bv);
			av = av_bv.a;
			bv = av_bv.b;

			// Calculates SIMPLEC coefficients at the internal real faces						due, dve, dun, dvn, de, dn OK
			simplec_coefficients = get_internal_simplec_coefficients(nx, ny, re, rn, xe, ye, xk, yk, au, av,
																		due, dve, dun, dvn, de, dn);
			due = simplec_coefficients.due;
			dve = simplec_coefficients.dve;
			dun = simplec_coefficients.dun;
			dvn = simplec_coefficients.dvn;
			de = simplec_coefficients.de;
			dn = simplec_coefficients.dn;

			// Calculates the SIMPLEC coefficients at the boundary faces						due, dve, dun, dvn, de, dn OK
			simplec_coefficients = get_boundary_simplec_coefficients(nx, ny, xe, ye
									, due, dve, dun, dvn, de, dn); // InOutput: last six
			due = simplec_coefficients.due;
			dve = simplec_coefficients.dve;
			dun = simplec_coefficients.dun;
			dvn = simplec_coefficients.dvn;
			de = simplec_coefficients.de;
			dn = simplec_coefficients.dn;

			// Calculates g																		g OK
			for (i = 0; i < nx * ny; i++)
			{
				g[i] = 1.0 / (Rg*T[i]);
			}

			// Calculates the coefficients of the linear system for pressure correction			ap OK								ap - first two wrong (0 to e-311) OK
			// g, ro, Uce and Vcn used in this subroutine are those calculated in the previous iteraction
			// de and dn must be calculated with the coef. of the linear sytem for u and v from which
			// u* and v* are obtained.
			ap = get_p_coefficients(nx, ny, dt, rp, re, rn, Jp, Uce, Vcn, roe, ron, g, de, dn, ap);

			// Solves the linear system for u

			// LU decomposition																	dl9, du9 OK							du9 - one value rounded
			d9 = lu2d9(au, nxy, nx, ny, dl9, du9);
			dl9 = d9.lower;
			du9 = d9.upper;

			// Linear system solution															u OK
			u = fb2d9(au, dl9, du9, ny, nx, nxy, u, tol_u, nitm_u, bu, r, w, z);

			// Norm initialization
			norm = 0.0;

			// Calculates the norm of the residual of the linear system for u					norm OK								precision issue
			norm = norm_l1_9d(nx, ny, u, bu, au, norm);

			// Relative residual calculation													res_u OK							precision issue
			//res_u = norm / sum(abs(bu));
			res_u = norm / norm_l1_b(nx, ny, bu);

			// Solves the linear system for v

			// LU decomposition																	dl9, du9 OK							du9 - one value rounded
			d9 = lu2d9(av, nxy, nx, ny, dl9, du9);
			dl9 = d9.lower;
			du9 = d9.upper;

			// Linear system solution															v OK								rounding e-20
			v = fb2d9(av, dl9, du9, ny, nx, nxy, v, tol_u, nitm_u, bv, r, w, z);

			// Norm initialization
			norm = 0.0;

			// Calculates the norm of the residual of the linear system for v					norm OK								small precision issue (e-15)
			norm = norm_l1_9d(nx, ny, v, bv, av, norm);

			// Relative residual calculation													res_v OK							small precision issue (e-5)
			//res_v = norm / sum(abs(bv));
			res_v = norm / norm_l1_b(nx, ny, bv);

			// Uses velocities at nodes to calculate velocities at internal real faces			ue, ve, un, vn, Ucem, Vcn OK		ve,vn - precision e-16
			velocities_face = get_velocities_at_internal_faces(nx, ny, dt, rp, re, rn, xe, ye, xk, yk
								, xen, yen, xke, yke, Jp, cup, cvp, au, av, roa, p, u, v, uea, vea, una, vna
								, ue, ve, un, vn, Uce, Vcn);	//InOut
			ue = velocities_face.ue;
			ve = velocities_face.ve;
			un = velocities_face.un;
			vn = velocities_face.vn;
			Uce = velocities_face.Uce;
			Vcn = velocities_face.Vcn;

			// Calculates the velocities at boundary faces based on the boundary conditions		ue, ve, un, vn, Ucem, Vcn OK		ve,vn - precision e-16
			velocities_face = get_velocities_at_boundary_faces(nx, ny, UF, xe, ye, xk, yk, u, v
																, ue, ve, un, vn, Uce, Vcn);	//InOut
			ue = velocities_face.ue;
			ve = velocities_face.ve;
			un = velocities_face.un;
			vn = velocities_face.vn;
			Uce = velocities_face.Uce;
			Vcn = velocities_face.Vcn;

			// ro, Uce and Vcn are the incorrect ones (obtained with p*)						bp OK								bp - precision e-14
			// rom, Ucem and Vcnm are those of the previous iteraction
			bp = get_p_source(nx, ny, dt, rp, re, rn, Jp, roe, ron, roe, ron, ro, roa
								, Ucea, Uce, Vcna, Vcn, bp);

			for (itm = 1; itm <= imax; itm++)
			{
				// Coefficients of the linear system of pl for the fictitious volumes			ap, bp OK							bp - precision e-14
				ap_bp = set_bcp(nx, ny, alphae, betae, betan, gamman, Uce, Vcw, p, ap, bp);
				ap = ap_bp.a;
				bp = ap_bp.b;

				// Solves the linear system for pl

				// LU decomposition																dl5, du5 OK							du5 - precision e-14
				d5 = lu2d5(ap, nxy, nx, ny, dl5, du5);
				dl5 = d5.lower;
				du5 = d5.upper;

				// Solution of the linear system												pl OK								pl - precision e-10
				pl = fb2d5(ap, dl5, du5, ny, nx, nxy, pl, tol_p, nitm_p, bp, r, w, z);
			}

			// Norm initialization
			norm = 0.0;

			// Calculates the norm of the residuals												norm OK
			norm = norm_l1_5d(nx, ny, pl, bp, ap, norm);

			// Residual calculation
			res_p = norm;

			// Pressure correction  (SIMPLEC method)											ro, p OK
			p_ro = get_pressure_density_correction_with_pl(nx, ny, pl, g, ro, p); // InOutput: last two
			p = p_ro.p;
			ro = p_ro.ro;


			// Extrapolation of p to fictitious													p OK
			p = get_p_extrapolation_to_fictitious(nx, ny, PF, alphae, betae, betan, gamman, Uce, Vcw, p);

			// Density at nodes with the state equation											ro OK
			ro = get_density_at_nodes(nx, ny, Rg, p, T, ro);

			// Velocity correction (SIMPLEC method)												u, v OK								v - precision e-14
			u_v = get_u_v_at_real_nodes_with_pl(nx, ny, xe, ye, xk, yk, rp, pl, au, av, u, v);
			u = u_v.u;
			v = u_v.v;

			// Extrapolation is based on the boundary conditions								u, v OK								v - precision e-14
			u_v = get_u_v_extrapolation_to_fictitious(nx, ny, UF, xk, yk, alphae, betae, Uce, Vcw, u, v);
			u = u_v.u;
			v = u_v.v;

			// Velocity correction at faces (SIMPLEC method)									ue,ve,un,vn,Uce,Vce OK				ue,ve,un,vn,Uce,Vce - precision e-4 (un,vn a bit weird but ok)
			velocities_face = get_velocity_correction_at_faces_with_pl(nx, ny
										, due, dve, dun, dvn, de, dn, pl
										, ue, ve, un, vn, Uce, Vcn);      // InOutput
			ue = velocities_face.ue;
			ve = velocities_face.ve;
			un = velocities_face.un;
			vn = velocities_face.vn;
			Uce = velocities_face.Uce;
			Vcn = velocities_face.Vcn;

			// Calculates density at faces using the corrected density and velocities at nodes
			roe_ron = get_density_at_faces(nx, ny, beta, ro, Uce, Vcn, roe, ron);
			roe = roe_ron.roe;															//		roe, ron OK							roe, ron - precision e-4
			ron = roe_ron.ron;

			// Coefficients and source of the temperature linear system							bt OK								bt - precision e+5
			at_bt = get_T_coefficients_and_source(nx, ny, beta, dt, rp, re, rn
										, xe, ye, xk, yk, Jp, roe, ron
										, roa, p, pa, cp, Uce, Vcn, u, v, T, Ta
										, at, bt);
			at = at_bt.a;
			bt = at_bt.b;

			// Calculating the contravariant velocity Vcw on the east boundary					Vcw OK								Vcw - precision e-8
			Vcw = get_Vcw(nx, ny, xe, xk, ue, ve, Vcw);

			// Coefficients of the linear system of T for the fictitious volumes				at, bt OK							at, bt- precision e+5
			at_bt = set_bcT(nx, ny, TF, Uce, Vcw, T, alphae, betae, betan, gamman, at, bt);
			at = at_bt.a;
			bt = at_bt.b;

			// Solves the linear system for T

			// LU decomposition																	dl9, du9 OK							dl9, du9 - precision e+4, e-1
			d9 = lu2d9(at, nxy, nx, ny, dl9, du9);
			dl9 = d9.lower;
			du9 = d9.upper;

			// Linear system solution															T OK								T - precision e+1
			T = fb2d9(at, dl9, du9, ny, nx, nxy, T, tol_u, nitm_u, bt, r, w, z);

			// Norm initialization
			norm = 0.0;

			// Calculates the norm of the residual of the linear system for T					norm better							norm - seems better, lower
			norm = norm_l1_9d(nx, ny, T, bt, at, norm);

			// Relative residual calculation													res_T better						res_T - seems better, lower
			//res_t = norm / sum(abs(bt));
			res_T = norm / norm_l1_b(nx, ny, bt);

			// Density at nodes with the state equation											ro OK								ro - precision e-1
			ro = get_density_at_nodes(nx, ny, Rg, p, T, ro);

			// Calculates density at faces using the corrected density and velocities at nodes	roe, ron OK							roe, ron - precision e-1 and first real (ron) but OK
			roe_ron = get_density_at_faces(nx, ny, beta, ro, Uce, Vcn, roe, ron);
			roe = roe_ron.roe;
			ron = roe_ron.ron;

			res = res_u + res_v + res_T + res_p;											  //res and other OK					all lower or good

			if (isnan(res))
			{
				std::cout << "NaN found. Stopping..." << std::endl;

				outputFile << "NaN found. Stopping..." << std::endl;

				exit(0);
			}
		}

		// Periodic data print
		if (it == 1)
		{
			Cdfi = get_cdfi(nx, ny, coord, Rg, PF, TF, UF, rb, yk, rn, p, Cdfi);

			std::cout << std::setw(10) << "it" << " " << std::setw(24) << "Residuals" << " " << std::setw(24) << "Cdfi" << " " << std::setw(10) << "it_stop" << std::endl;

			outputFile << std::setw(10) << "it" << " " << std::setw(24) << "Residuals" << " " << std::setw(24) << "Cdfi" << " " << std::setw(10) << "it_stop" << std::endl;

			outputFile << std::setw(10) << it << " " << std::setw(24) << std::scientific << std::setprecision(16) << res << " "  << std::setw(24) << Cdfi << " " << std::setw(10) << it_stop << std::endl;
		}

		if (it % wlf == 0)
		{
			Cdfi = get_cdfi(nx, ny, coord, Rg, PF, TF, UF, rb, yk, rn, p, Cdfi);

			std::cout << std::setw(10) << it << " " << std::setw(24) << std::scientific << std::setprecision(16) << res << " " << std::setw(24) << Cdfi << " " << std::setw(10) << it_stop << std::endl;

			outputFile << std::setw(10) << it << " " << std::setw(24) << std::scientific << std::setprecision(16) << res << " " << std::setw(24) << Cdfi << " " << std::setw(10) << it_stop << std::endl;
		}

		// Updating residual mean
		res_mean += res;

		// Checking for convergence
		if (it % nit_res == 0)
		{
			res_mean = res_mean / (double)(nit_res);

			if (res_mean < tol_res  &&  res_check == 0)
			{
				it_stop = 2 * it;

				res_check = 1;
			}

			if (it == it_stop) break;

			res_mean = 0.0;
		}
	}

	// Finishing cpu time
	auto tcpu2 = std::chrono::steady_clock::now();
	std::chrono::duration<double, std::milli> tcpu = tcpu2 - tcpu1;
	std::cout << "Elapsed time:" << " " << tcpu.count() << std::endl;

	// Closing file of residuals
	outputFile.close();

	// Calculating the RAM memory used
	//RAM = get_RAM(sim_id);

	// Calculating gamma from cp
	gcp = set_gamma(nx, ny, Rg, cp, gcp);

	// Post processing data
	post_processing(nx, ny, it, sem_a, sem_g, w_g, w_cam
		, sim_id, tcpu, RAM, res, lr, rb, x, y, xe, ye, xk, yk, Jp, Rg, cp
		, gcp, Uce, Vcn, de, dn, pl, bp, xp, yp, u, v, p, T, ro, Cdfi);

	// DELETAR REMOVER
	DGIoutputFile.close();
}
