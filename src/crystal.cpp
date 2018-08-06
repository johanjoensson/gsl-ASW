#include "crystal.h"
#include <algorithm>
#include <unordered_set>
#include "../../GSL-lib/src/special_functions.h"
#include "ewald_int.h"

namespace std {
	template<>
	struct hash<GSL::Vector>{
		size_t operator()(const GSL::Vector &v) const
		{
			return std::hash<double>()((v*v).norm()) ^ std::hash<double>()(v.norm());
		}
	};
}


bool comp_norm(GSL::Vector& a, GSL::Vector& b)
{
    return a.norm() < b.norm();
}

double bisect_q(double tol, double kappa, double eta, lm l, double q_min, double q_max)
{
	double q = 0;
	double ql = q_min, qu = q_max;
	while(qu - ql > tol){
		q = (qu + ql)/2;
		if((-q*q + eta*(l.l*GSL::log(q) + GSL::log(1 + kappa*kappa) - GSL::log(q*q + kappa*kappa) - GSL::log(tol)) + 1).val > 0){
			ql = q;
		}else{
			qu = q;
		}
	}

	return q;
}

size_t Crystal::calc_nk(double tol, double kappa, lm l)
{
	double eta = GSL::exp(GSL::log(6.5) + 2./3*GSL::log(4*M_PI/3) -
	2./3*GSL::log(this->volume)).val;

	int sign = 1.;

	if(-eta*GSL::log(tol).val < 0){
		sign = -1;
	}

	double q = 1;
	while((-q*q + eta*(l.l*GSL::log(q) + GSL::log(1 + kappa*kappa) - GSL::log(q*q + kappa*kappa) - GSL::log(tol)) + 1).val*sign > 0)
	{
		q += sqrt(eta);
	}
	q = bisect_q(tol, kappa, eta, l, q - sqrt(eta), q);

	return std::ceil(4*M_PI/3 * GSL::pow_int(q, 3)/GSL::pow_int(2*M_PI, 3)*this->volume);
}

double bisect_r(double tol, double kappa, double eta, lm l, double r_min, double r_max)
{
	double r = 0;
	double rl = r_min, ru = r_max;

	Ewald_integral I;
	I.set_kappa(kappa);
	I.set_ewald_param(eta);

	while(ru - rl > tol){
		r = (ru + rl)/2;
		if((l.l*GSL::log(r) + GSL::log(I.ewald_int(l, r)) - GSL::log(I.ewald_int(l, 1.)) - GSL::log(tol)).val > 0){
			rl = r;
		}else{
			ru = r;
		}
	}

	return r;
}

size_t Crystal::calc_nr(double tol, double kappa, lm l)
{
	double eta = GSL::exp(GSL::log(6.5) + 2./3*GSL::log(4*M_PI/3) -
	2./3*GSL::log(this->volume)).val;

	Ewald_integral I;
	I.set_kappa(kappa);
	I.set_ewald_param(eta);

	int sign = 1.;

	if(-GSL::log(tol).val < 0){
		sign = -1;
	}

	double r = 1;
	while((l.l*GSL::log(r) + GSL::log(I.ewald_int(l, r)) - GSL::log(I.ewald_int(l, 1.)) - GSL::log(tol)).val*sign > 0)
	{
		r += 2/sqrt(eta);
	}
	r = bisect_r(tol, kappa, eta, l, r - 2/sqrt(eta), r);

	return std::ceil(4*M_PI/3 * GSL::pow_int(r, 3)/this->volume);

}

void Crystal::calc_Rn(size_t num)
{
    std::cout << "Calculating lattice vectors" << std::endl;
    GSL::Vector a1(3);
    GSL::Vector a2(3);
    GSL::Vector a3(3);

    a1 = this->lat.lat[0]*this->lat.scale;
    a2 = this->lat.lat[1]*this->lat.scale;
    a3 = this->lat.lat[2]*this->lat.scale;

    std::vector<GSL::Vector> possibles;
    std::unordered_set<GSL::Vector> tmp, res;
    possibles = {a1, a2, a3};
    size_t n = tmp.size(), prev = 0;
	bool resort = true;
    std::cout << "0%| ";
	bool done = false;
    while(n < num || !done){
        n = tmp.size();
		if(resort){
			std::sort(possibles.begin(), possibles.end(), comp_norm);
			resort = !resort;
		}
        tmp.insert(possibles[0]);
        for(std::unordered_set<GSL::Vector>::const_iterator it = tmp.begin(); it != tmp.end(); it++){
            /* Only consider vectors not already found */
            if(tmp.find(possibles[0] + *it) == tmp.end()){
                possibles.push_back(possibles[0] + *it);
            }
        }
		if(possibles[0].norm() != possibles[1].norm()){
			if(n >= num){
				done = true;
			}
			resort = true;
		}
        possibles.erase(possibles.begin());
        if(n - prev > 0.10*num){
            std::cout << "==" << std::flush;
            prev = n;
        }
    }
    std::cout << "== |100%" << std::endl << std::endl;
    /*
       Above we only find linear combinations with positive integer coefficients
       we can generate all other linear combinations with the same length by
       simply mirroring our vectors along the three planes defned by the three
       lattice vectors.
    */
	GSL::Vector e1(3), e2(3), e3(3);
	e1[0] = 1.0;
	e2[1] = 1.0;
	e3[2] = 1.0;
	GSL::Matrix neg_x(3,3), neg_y(3,3), neg_z(3,3);
	neg_x[0] = -e1;
	neg_x[1] = e2;
	neg_x[2] = e3;

	neg_y[0] = e1;
	neg_y[1] = -e2;
	neg_y[2] = e3;

	neg_z[0] = e1;
	neg_z[1] = e2;
	neg_z[2] = -e3;
    for(GSL::Vector R : tmp){
        res.insert(R);
		if(R[0] > 0){
			res.insert(neg_x*R);
			if(R[1] > 0){
				res.insert(neg_y*R);
				res.insert(neg_y*neg_x*R);
				if(R[2] > 0){
					res.insert(neg_z*R);
					res.insert(neg_z*neg_y*R);
					res.insert(neg_z*neg_y*neg_x*R);
				}
			}else if(R[2] > 0){
				res.insert(neg_z*R);
				res.insert(neg_z*neg_x*R);
			}
		}else if(R[1] > 0){
			res.insert(neg_y*R);
			if(R[2] > 0){
				res.insert(neg_z*R);
				res.insert(neg_z*neg_y*R);
			}
		}else if(R[2] > 0){
			res.insert(neg_z*R);
		}
    }
	for(GSL::Vector R : res){
		Rn_vecs.push_back(R);
	}
    std::sort(Rn_vecs.begin(), Rn_vecs.end(), comp_norm);
    std::cout << "Number of lattice vectors : " << Rn_vecs.size() <<
	std::endl << std::endl;
}

void Crystal::calc_Kn(size_t num)
{
    std::cout << "Calculating reciprocal lattice vectors" << std::endl;
    GSL::Vector b1(3);
    GSL::Vector b2(3);
    GSL::Vector b3(3);

    b1 = this->lat.r_lat[0]/this->lat.scale;
    b2 = this->lat.r_lat[1]/this->lat.scale;
    b3 = this->lat.r_lat[2]/this->lat.scale;

    std::vector<GSL::Vector> possibles;
    std::unordered_set<GSL::Vector> tmp, res;
    possibles = {b1, b2, b3};
    size_t n = tmp.size(), prev = 0;
	bool resort = true, done = false;
    std::cout << "0%| ";
    while(n < num || !done){
        n = tmp.size();
		if(resort){
			std::sort(possibles.begin(), possibles.end(), comp_norm);
			resort = ! resort;
		}
        tmp.insert(possibles[0]);
        for(std::unordered_set<GSL::Vector>::const_iterator it = tmp.begin(); it != tmp.end(); it++){
            if(tmp.find(possibles[0] + *it) == tmp.end()){
                possibles.push_back(possibles[0] + *it);
            }
        }
		if(possibles[0].norm() != possibles[1].norm()){
			if(n >= num){
				done = true;
			}
			resort = true;
		}
        possibles.erase(possibles.begin());
        if((n-prev) > 0.10*num){
            std::cout << "==" << std::flush;
            prev = n;
        }
    }
    std::cout << "== |100%" << std::endl << std::endl;
    /*
       Above we only find linear combinations with positive integer coefficients
       we can generate all other linear combinations with the same length by
       simply mirroring our vectors in different planes
       lattice vectors.
    */
	GSL::Vector e1(3), e2(3), e3(3);
	e1[0] = 1.0;
	e2[1] = 1.0;
	e3[2] = 1.0;
	GSL::Matrix neg_x(3,3), neg_y(3,3), neg_z(3,3);
	neg_x[0] = -e1;
	neg_x[1] = e2;
	neg_x[2] = e3;

	neg_y[0] = e1;
	neg_y[1] = -e2;
	neg_y[2] = e3;

	neg_z[0] = e1;
	neg_z[1] = e2;
	neg_z[2] = -e3;
    for(GSL::Vector G : tmp){
        res.insert(G);
		if(G[0] > 0){
			res.insert(neg_x*G);
			if(G[1] > 0){
				res.insert(neg_y*G);
				res.insert(neg_y*neg_x*G);
				if(G[2] > 0){
					res.insert(neg_z*G);
					res.insert(neg_z*neg_y*G);
					res.insert(neg_z*neg_y*neg_x*G);
				}
			}
		}else if(G[1] > 0){
			res.insert(neg_y*G);
			if(G[2] > 0){
				res.insert(neg_z*G);
				res.insert(neg_z*neg_y*G);
			}
		}else if(G[2] > 0){
			res.insert(neg_z*G);
		}
    }
	for(GSL::Vector G : res){
		Kn_vecs.push_back(G);
	}
    std::sort(Kn_vecs.begin(), Kn_vecs.end(), comp_norm);
    std::cout << "Number of reciprocal lattice vectors : " << Kn_vecs.size() <<
	std::endl << std::endl;
}

double Crystal::calc_volume()
{
    GSL::Vector tmp;
    tmp = GSL::cross(lat.lat[1], lat.lat[2]);
    return std::abs(GSL::dot(lat.lat[0], tmp)*GSL::pow_int(lat.scale, 3));
}

double Crystal::calc_bz_volume()
{
    GSL::Vector tmp;
    tmp = GSL::cross(lat.r_lat[1], lat.r_lat[2]);
    return std::abs(GSL::dot(lat.r_lat[0], tmp)/GSL::pow_int(lat.scale, 3));
}

Crystal::Crystal()
 : Rn_vecs(), Kn_vecs(), lat(), atoms()
{
}

Crystal::Crystal(double& a)
 : Rn_vecs(), Kn_vecs(), atoms()
{
	GSL::Vector a1(3), a2(3), a3(3);
	a1[0] = a;
	a2[1] = a;
	a3[2] = a;
	lat = Lattice(a1, a2, a3);

    volume = lat.volume*GSL::pow_int(lat.scale, 3);

    bz_volume = lat.bz_volume/GSL::pow_int(lat.scale, 3);
}

Crystal::Crystal(const double& a)
 : Rn_vecs(), Kn_vecs(), atoms()
{
	GSL::Vector a1(3), a2(3), a3(3);
	a1[0] = a;
	a2[1] = a;
	a3[2] = a;
	lat = Lattice(a1, a2, a3);

    volume = lat.volume*GSL::pow_int(lat.scale, 3);

    bz_volume = lat.bz_volume/GSL::pow_int(lat.scale, 3);
}

Crystal::Crystal(const double& a, const double& b, const double& c)
 : Rn_vecs(), Kn_vecs(), atoms()
{
	GSL::Vector a1(3), a2(3), a3(3);
	a1[0] = a;
	a2[1] = b;
	a3[2] = c;
	lat = Lattice(a1, a2, a3);

    volume = lat.volume*GSL::pow_int(lat.scale, 3);

    bz_volume = lat.bz_volume/GSL::pow_int(lat.scale, 3);
}

Crystal::Crystal(double& a, double& b, double& c)
 : Rn_vecs(), Kn_vecs(), atoms()
{
	GSL::Vector a1(3), a2(3), a3(3);
	a1[0] = a;
	a2[1] = b;
	a3[2] = c;
	lat = Lattice(a1, a2, a3);

    volume = lat.volume*GSL::pow_int(lat.scale, 3);

    bz_volume = lat.bz_volume/GSL::pow_int(lat.scale, 3);
}

Crystal::Crystal(const GSL::Vector& a, const GSL::Vector& b, const GSL::Vector& c)
 : Rn_vecs(), Kn_vecs(), lat(a, b, c), atoms()
{
    volume = lat.volume*GSL::pow_int(lat.scale, 3);

    bz_volume = lat.bz_volume/GSL::pow_int(lat.scale, 3);

}

Crystal::Crystal(GSL::Vector& a, GSL::Vector& b, GSL::Vector& c)
 : Rn_vecs(), Kn_vecs(), lat(a, b, c), atoms()
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
    return a.get_pos().norm() < b.get_pos().norm();
}

std::vector<std::vector<Atom>> Crystal::calc_nearest_neighbours()
{
	std::vector<std::vector<Atom>> res(this->atoms.size());
	std::vector<Atom> tmp;
	GSL::Vector ri(3), rj(3);
	for(size_t i = 0; i < this->atoms.size(); i++){
		tmp.clear();
		ri = this->atoms[i].get_pos();
		tmp.push_back(this->atoms[i]);
		tmp.back().set_pos(GSL::Vector(3));
		for(size_t j = i + 1; j < this->atoms.size(); j++){
			rj = this->atoms[j].get_pos();
			tmp.push_back(this->atoms[j]);
			tmp.back().set_pos(ri - rj);
		}
		for(Atom a : tmp){
			if(a.get_pos() != GSL::Vector(3)){
				res[i].push_back(a);
			}
			for(GSL::Vector R : this->Rn_vecs){
				res[i].push_back(a);
				res[i].back().set_pos(a.get_pos() + R);
			}
		}
		std::sort(res[i].begin(), res[i].end(), comp_norm_at);
	}
	return res;
}
