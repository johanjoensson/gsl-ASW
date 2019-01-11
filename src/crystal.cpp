#include "crystal.h"
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include "spherical_fun.h"
#include "GSLpp/special_functions.h"
#include "ewald_int.h"

namespace std {
	template<>
	struct hash<GSL::Vector>{
		size_t operator()(const GSL::Vector &v) const
		{
			return std::hash<double>()((v*v).norm()) ^
			std::hash<double>()(v.norm());
		}
	};
	template<>
	struct hash<Atom>{
		size_t operator()(const Atom &a) const
		{
			return std::hash<GSL::Vector>()(a.pos) ^
			std::hash<size_t>()(a.get_Z());
		}
	};
}


bool comp_norm(GSL::Vector& a, GSL::Vector& b)
{
    return a.norm() < b.norm();
}

double bisect_q(double tol, double kappa, double eta, lm l, double q_min,
	double q_max)
{
	double q = 0;
	double ql = q_min, qu = q_max;
	while(qu - ql > tol){
		q = (qu + ql)/2;
		if((-q*q + eta*(l.l*GSL::log(q) + GSL::log(1 + kappa*kappa) -
		GSL::log(q*q + kappa*kappa) - GSL::log(tol)) + 1).val > 0){
			ql = q;
		}else{
			qu = q;
		}
	}

	return q;
}

double Crystal::calc_eta() const
{
	return GSL::exp(GSL::log(6.5) + 2./3*GSL::log(4*M_PI/3) -
			2./3*GSL::log(this->volume)).val;
}

size_t Crystal::calc_nk(double tol, double kappa, lm l)
{

	double eta = calc_eta();
	double q = 1;
	while((-q*q + eta*(l.l*GSL::log(q) + GSL::log(1 + kappa*kappa) -
	GSL::log(q*q + kappa*kappa) - GSL::log(tol)) + 1).val > 0)
	{
		q += sqrt(eta);
	}
	q = bisect_q(tol, kappa, eta, l, q - sqrt(eta), q);

	return static_cast<size_t>(4*M_PI/3 * GSL::pow_int(q, 3)/
	GSL::pow_int(2*M_PI, 3)*this->volume);
}

double bisect_r(double tol, double kappa, double eta, lm l, double r_min,
	double r_max)
{
	double r = 0;
	double rl = r_min, ru = r_max;

	Ewald_integral I;
	I.set_kappa(kappa);
	I.set_ewald_param(eta);

	while(std::abs(ru - rl) > tol){
		r = (ru + rl)/2;
		if((l.l*GSL::log(r) + GSL::log(I.ewald_int(l, r)) -
		GSL::log(I.ewald_int(l, 1.)) - GSL::log(tol)).val > 0){
			rl = r;
		}else{
			ru = r;
		}
	}

	return r;
}

size_t Crystal::calc_nr(double tol, double kappa, lm l)
{
	double eta = calc_eta();

	Ewald_integral I;
	I.set_kappa(kappa);
	I.set_ewald_param(eta);

	double r = 1;
	while((l.l*GSL::log(r) + GSL::log(I.ewald_int(l, r)) -
	GSL::log(I.ewald_int(l, 1.)) - GSL::log(tol)).val > 0)
	{
		r += 2/sqrt(eta);
	}
	r = bisect_r(tol, kappa, eta, l, r - 2/sqrt(eta), r);

	return static_cast<size_t>(4*M_PI/3 * GSL::pow_int(r, 3)/this->volume);

}

double Crystal::calc_Rmax(double tol, double kappa, lm l)
{
	double eta = calc_eta();

	Ewald_integral I;
	I.set_kappa(kappa);
	I.set_ewald_param(eta);

	double r = 1;
	while((l.l*GSL::log(r) + GSL::log(I.ewald_int(l, r)) -
	GSL::log(I.ewald_int(l, 1.)) - GSL::log(tol)).val > 0)
	{
		r += 2./sqrt(eta);
	}
	r = bisect_r(tol, kappa, eta, l, r - 2./sqrt(eta), r);

	return r;

}
void Crystal::set_Rn(double Rmax)
{
    GSL::Vector a1(this->lat.lat[0]*this->lat.scale);
    GSL::Vector a2(this->lat.lat[1]*this->lat.scale);
    GSL::Vector a3(this->lat.lat[2]*this->lat.scale);
    GSL::Vector b1(this->lat.r_lat[0]/this->lat.scale);
    GSL::Vector b2(this->lat.r_lat[1]/this->lat.scale);
    GSL::Vector b3(this->lat.r_lat[2]/this->lat.scale);

	int N1 = static_cast<int>(b1.norm()/(2*M_PI)*Rmax);
	int N2 = static_cast<int>(b2.norm()/(2*M_PI)*Rmax);
	int N3 = static_cast<int>(b3.norm()/(2*M_PI)*Rmax);

	// std::unordered_set<GSL::Vector> tmp;
	for(int n1 = -N1; n1 <= N1; n1++){
		for(int n2 = -N2; n2 <= N2; n2++){
			for(int n3 = -N3; n3 <= N3; n3++){
				Rn_vecs.push_back(n1*a1 + n2*a2 + n3*a3);
			}
		}
	}
	// Rn_vecs.assign(tmp.begin(), tmp.end());
	std::sort(Rn_vecs.begin(), Rn_vecs.end(), comp_norm);
}

double Crystal::calc_Kmax(double tol, double kappa, lm l)
{
	double eta = calc_eta();

	double q = 1;
	while((-q*q + eta*(l.l*GSL::log(q) + GSL::log(1 + kappa*kappa) -
	GSL::log(q*q + kappa*kappa) - GSL::log(tol)) + 1).val > 0)
	{
		q += sqrt(eta);
	}
	q = bisect_q(tol, kappa, eta, l, q - sqrt(eta), q);

	return q;
}

void Crystal::set_Kn(double Kmax)
{
    GSL::Vector a1(this->lat.lat[0]*this->lat.scale);
    GSL::Vector a2(this->lat.lat[1]*this->lat.scale);
    GSL::Vector a3(this->lat.lat[2]*this->lat.scale);
    GSL::Vector b1(this->lat.r_lat[0]/this->lat.scale);
    GSL::Vector b2(this->lat.r_lat[1]/this->lat.scale);
    GSL::Vector b3(this->lat.r_lat[2]/this->lat.scale);

	int N1 = static_cast<int>(a1.norm()/(2*M_PI)*Kmax);
	int N2 = static_cast<int>(a2.norm()/(2*M_PI)*Kmax);
	int N3 = static_cast<int>(a3.norm()/(2*M_PI)*Kmax);


	for(int n1 = -N1; n1 <= N1; n1++){
		for(int n2 = -N2; n2 <= N2; n2++){
			for(int n3 = -N3; n3 <= N3; n3++){
				Kn_vecs.push_back(n1*b1 + n2*b2 + n3*b3);
			}
		}
	}
	std::sort(Kn_vecs.begin(), Kn_vecs.end(), comp_norm);
}

double Crystal::calc_volume() 
{
    GSL::Vector tmp, tmp1 = lat.lat[1], tmp2 = lat.lat[2];
    tmp = GSL::cross(tmp1, tmp2);
    return std::abs(GSL::dot(lat.lat[0], tmp)*GSL::pow_int(lat.scale, 3));
}

double Crystal::calc_bz_volume() 
{
    GSL::Vector tmp;
    tmp = GSL::cross(lat.r_lat[1], lat.r_lat[2]);
	std::cout << tmp << std::endl;
    return std::abs(GSL::dot(lat.r_lat[0], tmp)/GSL::pow_int(lat.scale, 3));
}

Crystal::Crystal()
 : Rn_vecs(), Kn_vecs(), lat(), atoms(), volume(), bz_volume()
{
}


Crystal::Crystal(const double& a)
 : Rn_vecs(), Kn_vecs(), lat({a, 0., 0.},{0., a, 0.},{0., 0., a}), atoms(),
   volume(0), bz_volume(0)
{
    volume = lat.volume*GSL::pow_int(lat.scale, 3);

    bz_volume = lat.bz_volume/GSL::pow_int(lat.scale, 3);
}

Crystal::Crystal(const double& a, const double& b, const double& c)
 : Rn_vecs(), Kn_vecs(), lat({a, 0., 0.},{0., b, 0.},{0., 0., c}), atoms(),
   volume(0), bz_volume(0)
{
	GSL::Vector a1(3), a2(3), a3(3);
	a1[0] = a;
	a2[1] = b;
	a3[2] = c;
	lat = Lattice(a1, a2, a3);

    volume = lat.volume*GSL::pow_int(lat.scale, 3);

    bz_volume = lat.bz_volume/GSL::pow_int(lat.scale, 3);
}


Crystal::Crystal(const GSL::Vector& a, const GSL::Vector& b,
	const GSL::Vector& c)
 : Rn_vecs(), Kn_vecs(), lat(a, b, c), atoms(), volume(0), bz_volume(0)
{
    volume = lat.volume*GSL::pow_int(lat.scale, 3);

    bz_volume = lat.bz_volume/GSL::pow_int(lat.scale, 3);

}

void Crystal::add_atoms(const std::vector<Atom>& v)
{
    for(Atom a : v){
        this->atoms.push_back(a);
    }
}

GSL::Vector& Crystal::get_Rn(const size_t& i)
{
    return Rn_vecs[i];
}

bool comp_norm_at(Atom& a, Atom& b)
{
    return a.pos.norm() < b.pos.norm();
}

std::vector<std::vector<Atom>> Crystal::calc_nearest_neighbours()
{
	std::vector<std::vector<Atom>> res(this->atoms.size());
	std::vector<Atom> tmp;
	// Container for storing all neighbours, in all shells
	std::unordered_set<Atom> pre_res;
	Atom tmp_at;
	GSL::Vector ri(3), rj(3);
	for(size_t i = 0; i < this->atoms.size(); i++){

		tmp.clear();
		pre_res.clear();
		ri = this->atoms[i].pos;
		tmp.push_back(this->atoms[i]);
		// Calculate all neighbours inside the cell
		for(size_t j = 0; j < this->atoms.size(); j++){
			if(i != j){
				rj = this->atoms[j].pos;
				tmp.push_back(this->atoms[j]);
			}
		}

		for(Atom a : tmp){
			// Add all lattice vectors
			for(GSL::Vector R : this->Rn_vecs){
				// Insert all atoms in the system, except fot the original one
				tmp_at = a;
				tmp_at.set_pos(a.pos + R - ri);
				if(tmp_at.get_pos() != GSL::Vector(3)){
					pre_res.insert(tmp_at);
				}
			}
		}
		res[i].assign(pre_res.begin(), pre_res.end());
		std::sort(res[i].begin(), res[i].end(), comp_norm_at);
	}

	return res;
}
