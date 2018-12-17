#include <cmath>
#include <algorithm>
#include "numerov_solver.h"

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#endif

Numerov_solver::Numerov_solver()
{}

std::vector<double> Numerov_solver::solve_left(Logarithmic_mesh &mesh,
	std::vector<double> &v, size_t i_lr, std::vector<double> &init_cond, double E)
{
	size_t l = init_cond.size();
	std::vector<double> res(mesh.r.size(),0);
	double g, g1, g2;
	double A = mesh.A();

	size_t offset = res.size() - l - 1;
	for(size_t i = 0; i < l; i++){
		res[offset - i] = init_cond[i];
	}

	for(size_t i = offset; i >= i_lr; i--){
		g = mesh.drx[i]*mesh.drx[i]*(E - v[i]) - A*A/4.;
		g1 = mesh.drx[i+1]*mesh.drx[i+1]*(E - v[i+1]) - A*A/4.;
		g2 = mesh.drx[i+2]*mesh.drx[i+2]*(E - v[i+2]) - A*A/4.;
		res[i] = (2*res[i+1]*(1 - 5*g1/12) - res[i+2]*(1 + g2/12))/(1 + g/12);
	}

	return res;
}

std::vector<double> Numerov_solver::solve_right(Logarithmic_mesh &mesh,
	std::vector<double> &v, size_t i_lr, std::vector<double> &init_cond, double E)
{
	size_t l = init_cond.size();
	std::vector<double> res(mesh.r.size(),0);
	double g, g1, g2;

	for(size_t i = 0; i < l; i++){
		res[i] = init_cond[i];
	}

	double A = mesh.A();
	for(size_t i = l; i <= i_lr; i++){
		g = mesh.drx[i]*mesh.drx[i]*(E - v[i]) - A*A/4.;
		g1 = mesh.drx[i-1]*mesh.drx[i-1]*(E - v[i-1]) - A*A/4.;
		g2 = mesh.drx[i-2]*mesh.drx[i-2]*(E - v[i-2]) - A*A/4.;
		res[i] = (2*res[i-1]*(1 - 5*g1/12) - res[i-2]*(1 + g2/12))/(1 + g/12);
	}

	return res;
}

int sign(double f)
{
    int res = 1;
    if(f < 0){
	    res = -1;
    }
    return res;
}

size_t Numerov_solver::find_inversion_point(Logarithmic_mesh &mesh,
	std::vector<double> &v, double e_trial)
{
	size_t inv = mesh.r.size() - 3;
	double pot;
	bool done = false;
	int s = sign(v[inv] - e_trial);
	for(size_t i = mesh.r.size() - 3; i > 1   && !done; i--){
		pot = v[i];
		if(s != sign(pot - e_trial)){
			inv = i;
			done = true;
		}
	}

	return inv;
}

int count_nodes(std::vector<double> fun, size_t i_inv)
{
	int n_nodes = 0;
	for(size_t i = 0; i < i_inv - 1; i++){
		n_nodes += (1 - sign(fun[i + 1])/sign(fun[i]))/2;
	}
	return n_nodes;
}

double Numerov_solver::variational_energy_correction(Logarithmic_mesh &mesh,
	std::vector<double> &v, std::vector<double> &fun, size_t i_inv, double e_trial)
{
	double f_inv, f_m1, f_p1, cusp_val, df;
	double A = mesh.A();
	f_inv = 1 + (mesh.drx[i_inv]*mesh.drx[i_inv]*(e_trial - v[i_inv]) -
	A*A/4.)/12.;
	f_m1 = 1 + (mesh.drx[i_inv-1]*mesh.drx[i_inv-1]*(e_trial - v[i_inv-1]) -
	A*A/4.)/12.;
	f_p1 = 1 + (mesh.drx[i_inv+1]*mesh.drx[i_inv+1]*(e_trial - v[i_inv+1]) -
	A*A/4.)/12.;
	cusp_val = (fun[i_inv - 1]*f_m1 + fun[i_inv + 1]*f_p1 +
		10*fun[i_inv]*f_inv)/12.;
	df = f_inv*fun[i_inv] - f_inv*cusp_val;
	return 6.*df*cusp_val;
}

std::vector<double> Numerov_solver::solve(Logarithmic_mesh &mesh,
	std::vector<double> &v, std::vector<double> &l_init,
	std::vector<double> &r_init, double& en, int n_nodes)
{


    double e_max = v.back(), e_min = e_max;
	for(size_t i = 1; i < mesh.r.size(); i++){
		if(e_min > v[i] + 0.25/mesh.r2[i]){
			e_min = v[i] + 0.25/mesh.r2[i];
		}
	}
    double e_trial = en;
	if(e_trial < e_min){
		e_trial = e_min;
	}
#ifdef DEBUG
	std::ofstream out_file("numerov.debug", std::fstream::out | std::fstream::app);
	out_file << std::setprecision(18);
	out_file << std::string(80, '=') << std::endl;
	out_file << e_min << " < " << en << " < " << e_max << std::endl;
#endif
    bool done = false;
    std::vector<double> res(mesh.r.size(),0), tmp(mesh.r.size(),0);

    size_t i_inv = 0, it1 = 0, it2 = 0;
    int n_tmp = -1;

    double scale = 0;
    double de = 0;



    // Loop until we have found a good enough estimate of the energy
    while(!done){
	    // Rough estimate of the energy
	    // Make sure we have the correct number of nodes
	    it2 = 0;
	    while(n_tmp != n_nodes && it2 < 5){
		    i_inv = find_inversion_point(mesh, v, e_trial);
		    tmp = this->solve_right(mesh, v, i_inv, l_init, e_trial);
		    n_tmp = count_nodes(tmp, i_inv);
			de = 0.;
			if(n_tmp > n_nodes){
				e_max = e_trial;
				de = 0.5*(e_min - e_trial);
		    }else if(n_tmp < n_nodes){
				e_min = e_trial;
				de = 0.5*(e_max - e_trial);
		    }
			e_trial += de;
#ifdef DEBUG
			out_file << "# Nodes : " << n_tmp << ", required :  " << n_nodes << std::endl;
			out_file << "Inner iteration " << it2 << std::endl;
#endif
			it2++;
		}

		i_inv = find_inversion_point(mesh, v, e_trial);
	    // Correct number of nodes, now match values at mesh.r[i_inv]
	    for(size_t i = 0; i <= i_inv; i++){
		    res[i] = tmp[i];
	    }
	    tmp = this->solve_left(mesh, v, i_inv, r_init, e_trial);
		if(std::abs(tmp[i_inv]) > 0){
			scale = res[i_inv]/tmp[i_inv];
			for(size_t i = i_inv; i < tmp.size(); i++){
				res[i] = scale*tmp[i];
			}
		}

		// Finer corrections to the energy
	    // Match slope at mesh.r[i_inv]
		de = variational_energy_correction(mesh, v, res, i_inv, e_trial);
	    if(std::abs(de) < 1E-14){
		    done = true;
	    }

#ifdef DEBUG
		out_file << e_min << " < " << e_trial << " < "<< e_max << std::endl;
		out_file << "Outer iteration  " << it1;
		out_file << " E_T = " << e_trial;
		out_file << " de = " << de << std::endl;
		out_file << " Emax - Emin = " << e_max - e_min << std::endl;
#endif

		if(de > 0){
			e_min = e_trial;
		}else{
			e_max = e_trial;
		}

		if(std::abs(e_max - e_min) < 5E-14){
			if(de < 0){
				e_min -= 1.5;
			}else{
				e_max += 1.5;
			}
		}

		e_trial = 0.5*(e_max + e_min);
		n_tmp = -1;
	    it1++;
    }



    en = e_trial;
#ifdef DEBUG
		out_file.close();
#endif
    return res;
}
