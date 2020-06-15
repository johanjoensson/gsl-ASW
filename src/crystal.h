#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <vector>
#include <algorithm>
#include <unordered_set>
#include <set>
#include "lattice.h"
#include "site.h"
#include "GSLpp/vector.h"
#include "GSLpp/matrix.h"

template<size_t dim>
using Neighbours =  std::vector<Site_t<dim>>;

class Empty_t{
};

template<size_t dim, class Atom = Empty_t>
class Crystal_t {
	Lattice_t<dim> lat_m;
	std::vector<Site_t<dim>> sites_m;
	std::vector<size_t> atom_index_m;
	std::vector<Atom> atoms_m;


	std::vector<GSL::Vector> R_m, K_m;
	std::array<size_t, dim> size_m;
public:
	Crystal_t(const Lattice_t<dim>& lat):lat_m(lat), sites_m(), atom_index_m(), atoms_m(), R_m(), K_m(), size_m(){}
	void add_sites(const std::vector<GSL::Vector>&);

	void add_basis(const std::vector<Atom>&);

	void set_Rn(const double Rmax);
	void set_Kn(const double Kmax);
	void set_size(const std::array<size_t, dim>& size){size_m = size;}
	std::vector<Neighbours<dim>> calc_nearest_neighbours() const;
	std::vector<Neighbours<dim>> calc_nearest_neighbours(const size_t n_shells) const;
	std::vector<std::vector<Neighbours<dim>>> determine_nn_shells(const std::vector<Neighbours<dim>>& nn) const;

	Lattice_t<dim> lat() const {return lat_m;}
	std::vector<Site_t<dim>>& sites(){return sites_m;}
	std::vector<Site_t<dim>> sites() const {return sites_m;}

	Site_t<dim>& site(const size_t i) {return sites_m[i];}
	Site_t<dim> site(const size_t i) const {return sites_m[i];}

	std::vector<Atom>& atoms() {return atoms_m;}
	std::vector<Atom> atoms() const {return atoms_m;}
	std::array<size_t, dim>& size(){return size_m;}

	Atom& atom(const Site_t<dim>& s) {return atoms_m[atom_index_m[s.index()]];}
	Atom atom(const Site_t<dim>& s) const{return atoms_m[atom_index_m[s.index()]];}

	size_t& atom_index(const Site_t<dim>& s) {return atom_index_m[s.index()];}
	size_t atom_index(const Site_t<dim>& s) const {return atom_index_m[s.index()];}

	std::vector<GSL::Vector>& Rn_vecs(){return R_m;}
	std::vector<GSL::Vector> Rn_vecs() const {return R_m;}
	std::vector<GSL::Vector>& Kn_vecs(){return K_m;}
	std::vector<GSL::Vector> Kn_vecs() const {return K_m;}

	double volume() const {return lat().lat().det();}
};

template<size_t dim, class Atom>
void Crystal_t<dim, Atom>::add_sites(const std::vector<GSL::Vector>& positions)
{
	for(size_t i = 0; i < positions.size(); i++){
		sites_m.push_back(Site_t<dim>(i, positions[i]*lat_m.lat(), size_m));
	}
}

template<size_t dim, class Atom>
void Crystal_t<dim, Atom>::add_basis(const std::vector<Atom>& basis)
{
	for(const auto unit : basis){
		atoms_m.push_back(unit);
		atom_index_m.push_back(atoms_m.size() - 1);
	}
}

inline bool comp_norm(const GSL::Vector& a, const GSL::Vector& b)
{
    return a.norm<double>() < b.norm<double>();
}

template<size_t dim>
bool comp_norm_site(const Site_t<dim>& a, const Site_t<dim>& b)
{
	const GSL::Vector va = a.pos(), vb = b.pos();
	return comp_norm(va, vb);
}

struct Vector_comp_norm{
	bool operator()(const GSL::Vector& a, const GSL::Vector& b){
		return comp_norm(a, b);
	}

};

template<size_t dim, class Atom>
void Crystal_t<dim, Atom>::set_Rn(const double Rmax)
{
	std::array<int, dim> N, n, zero;
	N.fill(0);
	zero.fill(0);
	GSL::Matrix a(lat_m.lat()), b(lat_m.recip_lat());
	GSL::Vector tmp(dim);
	std::set<GSL::Vector, Vector_comp_norm> res;

	// Calculate limits
	for(size_t i = 0; i < dim; i++){
		N[i] = static_cast<int>(std::ceil(b[i].norm<double>()/(2*M_PI)*Rmax));
		n[i] = -N[i];
	}

	// temporary vector for storing linear combinations of new and old vector
	// std::unordered_set<GSL::Vector, GSL::Vector_hasher_t<double, gsl_vector, std::allocator<double>>> r_tmp;
	std::set<GSL::Vector, Vector_comp_norm>r_tmp;

	while(n != N){
/*
		if(n == zero){
			n.back()++;
			continue;
		}
*/
		tmp.assign(n.begin(), n.end());
		r_tmp.insert(tmp*a);
		n.back()++;
		for(size_t i = dim - 1; i > 0; i--){
			if(n[i] > N[i]){
				n[i - 1]++;
				n[i] = -N[i];
			}
		}
	}
	tmp.assign(n.begin(), n.end());
	r_tmp.insert(tmp*a);

	R_m.assign(r_tmp.begin(), r_tmp.end());
	// std::sort(R_m.begin(), R_m.end(), comp_norm);
}

template<size_t dim, class Atom>
void Crystal_t<dim, Atom>::set_Kn(const double Kmax)
{
	std::array<int, dim> N, n, zero;
	N.fill(0);
	zero.fill(0);
	GSL::Vector tmp(dim);
	GSL::Matrix a(lat_m.lat()), b(lat_m.recip_lat());
	std::set<GSL::Vector, Vector_comp_norm> res;

	// Calculate limits
	for(size_t i = 0; i < dim; i++){
		N[i] = static_cast<int>(std::ceil(a[i].norm<double>()/(2*M_PI)*Kmax));
		n[i] = -N[i];
	}

	// temporary vector for storing linear combinations of new and old vector
	std::set<GSL::Vector, Vector_comp_norm> k_tmp;
	while(n != N){
/*
		if(n == zero){
			n.back()++;
			continue;
		}
*/
		tmp.assign(n.begin(), n.end());
		k_tmp.insert(tmp*a);
		n.back()++;
		for(size_t i = dim - 1; i > 0; i--){
			if(n[i] > N[i]){
				n[i - 1]++;
				n[i] = -N[i];
			}
		}
	}
	tmp.assign(n.begin(), n.end());
	k_tmp.insert(tmp*a);

	K_m.assign(k_tmp.begin(), k_tmp.end());
	// std::sort(K_m.begin(), K_m.end(), comp_norm);
}

template<size_t dim, class Atom>
std::vector<Neighbours<dim>> Crystal_t<dim, Atom>::calc_nearest_neighbours() const
{
	std::vector<Neighbours<dim>> res(sites_m.size());
	GSL::Vector ri, rj/*, zero_v(dim)*/;
	for(size_t i = 0; i < sites_m.size(); i++){
		for(size_t j = i + 1; j < sites_m.size(); j++){
			// Add all lattice vectors
			for(const auto& R : R_m){
				// Insert all atoms in the system
				res[i].push_back(Site_t<dim>(j, sites_m[i].pos() - sites_m[j].pos() + R, size_m));
				res[j].push_back(Site_t<dim>(i, sites_m[j].pos() - sites_m[i].pos() - R, size_m));

			}
		}
		std::sort(res[i].begin(), res[i].end(), comp_norm_site<dim>);
	}
	return res;
}

template<size_t dim, class Atom>
std::vector<Neighbours<dim>> Crystal_t<dim, Atom>:: calc_nearest_neighbours(const size_t n_steps) const
{
	std::vector<Neighbours<dim>> res(sites_m.size());
	std::unordered_set<Site_t<dim>, Site_t_hasher<dim>> unique_maker;
	std::array<size_t, dim> stop, current, new_coords, flips, flip_stop({2,0});
	stop.fill(n_steps);
	bool periodic = (R_m.size() != 0), add;

	GSL::Matrix a = lat_m.lat();
	GSL::Vector R, rp, zerov(dim);
	for(size_t i = 0; i < sites().size(); i++){
		current.fill(0);
		while(current != stop){
			// Start by incrementing current
			current.back()++;
			for(auto tmp = ++current.rbegin(), tmp_stop = ++stop.rbegin(); tmp != current.rend(); tmp++, tmp_stop++){
				if(*(tmp - 1) > *(tmp_stop - 1)){
					(*tmp)++;
					*(tmp - 1) = 0;
				}
			}

			flips.fill(0);
			while(flips != flip_stop){
				R = GSL::Vector(dim, 0);;
				add = false;
				for(size_t j = 0; j < dim; j++){
					new_coords[j] = sites_m[i].coord()[j];
					if(current[j] == 0){
						continue;
					}
					if(flips[j] == 0){
						if(periodic){
							new_coords[j] = (sites_m[i].coord()[j] + current[j]) % size_m[j];
							rp.copy(a[j]);
							rp *= static_cast<double>((sites_m[i].coord()[j] + current[j]) / size_m[j]);
							R += rp;
							add = true;
						}else if(sites_m[i].coord()[j] + current[j] < size_m[j]){
							new_coords[j] = sites_m[i].coord()[j] + current[j];
							add = true;
						}
					}else{
						if(periodic){
							new_coords[j] = ((current[j]/size_m[j] + 1)*size_m[j] + sites_m[i].coord()[j] - current[j]) % size_m[j];
							rp.copy(a[j]);
							rp *= static_cast<double>(current[j]/size_m[j] + 1 - (((current[j]/size_m[j] + 1)*size_m[j] + sites_m[i].coord()[j] - current[j])/size_m[j]));
							R -= rp;
							add = true;
						}else if(sites_m[i].coord()[j] >= current[j]){
							new_coords[j] = sites_m[i].coord()[j] - current[j];
							add = true;
						}
					}
				}
				if(add){
					Site_t<dim> tmp(new_coords, zerov, size_m);
					rp = sites_m[tmp.index()].pos();
					rp += R;
					rp -= sites_m[i].pos();
					std::cout << rp << std::endl;
					tmp.set_pos(rp);
					res[i].push_back(tmp);
				}

				flips.back()++;
				for(size_t idx = dim - 1; idx > 0; idx--){
					if(flips[idx] > 1){
						flips[idx - 1]++;
						flips[idx] = 0;
					}
				}
			}
		}
		unique_maker = std::unordered_set<Site_t<dim>, Site_t_hasher<dim>> (res[i].begin(), res[i].end());
		res[i].assign(unique_maker.begin(), unique_maker.end());
		std::sort(res[i].begin(), res[i].end(), comp_norm_site<dim>);

	}
	return res;
}

template<size_t dim, class Atom>
std::vector<std::vector<Neighbours<dim>>> Crystal_t<dim, Atom>::determine_nn_shells(const std::vector<Neighbours<dim>>& nn) const
{
	std::vector<std::vector<Neighbours<dim>>> res(nn.size());
	size_t shell_idx;
	double old_dist;
	Site_t<dim> tmp;
	for(size_t site_idx = 0; site_idx < nn.size(); site_idx++){
		tmp = nn[site_idx][0];
		old_dist = tmp.pos(). template norm<double>();
		res[site_idx] = std::vector<Neighbours<dim>>(1, Neighbours<dim>());
		shell_idx = 0;
		for(size_t neighbour_idx = 0; neighbour_idx < nn[site_idx].size(); neighbour_idx++){
			tmp = nn[site_idx][neighbour_idx];
			if(std::abs( tmp.pos(). template norm<double>() - old_dist ) > 1e-6){
				old_dist = tmp.pos(). template norm<double>();
				shell_idx++;
				res[site_idx].push_back(Neighbours<dim>());
			}
			res[site_idx][shell_idx].push_back( tmp );
		}
	}

	return res;
}

#endif //CRYSTAL_H
