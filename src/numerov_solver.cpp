#include <cmath>
#include "numerov_solver.h"

double empty_pot(double r)
{
	return 0.*r;
}

Numerov_solver::Numerov_solver()
{
	this->v_at = &empty_pot;
	this->v_eff = &empty_pot;
}

void Numerov_solver::set_v_eff(double (*new_v_eff)(double r))
{
	this->v_eff = new_v_eff;
}

void Numerov_solver::set_v_at(double (*new_v_at)(double r))
{
	this->v_at = new_v_at;
}

std::vector<double> Numerov_solver::solve_left(Logarithmic_mesh &mesh, int i_lr, std::vector<double> &init_cond, double E)
{
	unsigned int l = init_cond.size();
	std::vector<double> res(mesh.r.size(),0);
	double g, g1, g2;

	for(unsigned int i = res.size() - 1; i > res.size() - 1 - l; i--){
		res[i] = init_cond[res.size() - 1 - i];
	}

	for(int i = res.size() - 1 - l; i >= i_lr; i--){
		g = mesh.drx[i]*mesh.drx[i]*(E - v_eff(mesh.r[i])) - mesh.A*mesh.A/4.;
		g1 = mesh.drx[i+1]*mesh.drx[i+1]*(E - v_eff(mesh.r[i+1])) - mesh.A*mesh.A/4.;
		g2 = mesh.drx[i+2]*mesh.drx[i+2]*(E - v_eff(mesh.r[i+2])) - mesh.A*mesh.A/4.;
		res[i] = (2*res[i+1]*(1 - 5*g1/12) - res[i+2]*(1 + g2/12))/(1 + g/12);
	}

	return res;
}

std::vector<double> Numerov_solver::solve_right(Logarithmic_mesh &mesh, int i_lr, std::vector<double> &init_cond, double E)
{
	unsigned int l = init_cond.size();
	std::vector<double> res(mesh.r.size(),0);
	double g, g1, g2;

	for(unsigned int i = 0; i < l; i++){
		res[i] = init_cond[i];
	}

	for(int i = l; i <= i_lr; i++){
		g = mesh.drx[i]*mesh.drx[i]*(E - v_eff(mesh.r[i])) - mesh.A*mesh.A/4.;
		g1 = mesh.drx[i-1]*mesh.drx[i-1]*(E - v_eff(mesh.r[i-1])) - mesh.A*mesh.A/4.;
		g2 = mesh.drx[i-2]*mesh.drx[i-2]*(E - v_eff(mesh.r[i-2])) - mesh.A*mesh.A/4.;
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

int Numerov_solver::find_inversion_point(Logarithmic_mesh &mesh, double e_trial)
{
	int inv = mesh.r.size() - 1;
	double pot;
	bool done = false;
	for(unsigned int i = 0; i < mesh.r.size() && !done; i++){
		pot = v_at(mesh.r[i]);
		if(pot > e_trial){
			inv = i;
			done = true;
		}
	}
	return inv;
}

int count_nodes(std::vector<double> &fun, int i_inv)
{
	int n_nodes = 0;
	for(int i = 0; i < i_inv; i++){
		n_nodes += 0.5*(1 - sign(fun[i + 1])/sign(fun[i]));
	}
	return n_nodes;
}

double Numerov_solver::variational_energy_correction(Logarithmic_mesh &mesh, std::vector<double> &fun, int i_inv, double e_trial)
{
	double f_inv, f_m1, f_p1, cusp_val, df;
	f_inv = 1 + (mesh.drx[i_inv]*mesh.drx[i_inv]*(e_trial - v_eff(mesh.r[i_inv])) - mesh.A*mesh.A/4.)/12.;
	f_m1 = 1 + (mesh.drx[i_inv-1]*mesh.drx[i_inv-1]*(e_trial - v_eff(mesh.r[i_inv-1])) - mesh.A*mesh.A/4.)/12.;
	f_p1 = 1 + (mesh.drx[i_inv+1]*mesh.drx[i_inv+1]*(e_trial - v_eff(mesh.r[i_inv+1])) - mesh.A*mesh.A/4.)/12.;
	cusp_val = (fun[i_inv - 1]*f_m1 + fun[i_inv + 1]*f_p1 + 10*fun[i_inv]*f_inv)/12.;
	df = f_inv*(fun[i_inv]/cusp_val - 1);
	return 12.*df*cusp_val*cusp_val/mesh.drx[i_inv];
}

std::vector<double> Numerov_solver::solve(Logarithmic_mesh &mesh, std::vector<double> &l_init, std::vector<double> &r_init, double &en, int n_nodes)
{
    double e_trial = en;
    double e_max = en + 2, e_min = en - 2;
    bool done = false;
    std::vector<double> res(mesh.r.size(),0), tmp(mesh.r.size(),0);

    int i_inv = mesh.r.size();
    int n_tmp = -1;

    double scale = 0;
    double de = 0;



    // Loop until we have found a good enough estimate of the energy
    while(!done){
	    // Rough estimate of the energy
	    // MAke sure we have the correct number of nodes
	    while(n_tmp != n_nodes){	    
		    i_inv = find_inversion_point(mesh, e_trial);
		    tmp = this->solve_right(mesh, i_inv, l_init, e_trial);
		    n_tmp = count_nodes(tmp, i_inv);
		    if(n_tmp > n_nodes){
			    e_max = e_trial;
			    e_trial = 0.5*(e_trial + e_min);
		    }else if(n_tmp < n_nodes){
			    e_min = e_trial;
			    e_trial = 0.5*(e_trial + e_max);
		    }
		    if(std::abs(e_max - e_min) < 1E-10){
			    e_min = e_max - 4;
			    e_trial = 0.5*(e_max + e_min);
		    }
	    }
	    // Correct number of nodes, now match values at mesh.r[i_inv]
	    for(int i = 0; i <= i_inv; i++){
		    res[i] = tmp[i];
	    }
	    tmp = this->solve_left(mesh, i_inv, r_init, e_trial);
	    scale = res[i_inv]/tmp[i_inv];
	    for(unsigned int i = i_inv + 1; i < res.size(); i++){
		    res[i] = scale*tmp[i];
	    }
	    // Finer corrections to the energy
	    // Match slope at mesh.r[i_inv]
	    de = variational_energy_correction(mesh, res, i_inv, e_trial);
	    if(std::abs(de) < 1E-10){
		    done = true;
	    }else{
		    e_trial += de;
		    n_tmp = -1; 
	    }
    }
    en = e_trial;
    return res;
}
