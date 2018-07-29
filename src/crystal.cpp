#include "crystal.h"
#include <algorithm>
#include <unordered_set>
#include "../../GSL-lib/src/special_functions.h"

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

void Crystal::calc_Rn(size_t num)
{
    std::cout << "Calculating lattice vectors" << std::endl;
    GSL::Vector a1(3);
    GSL::Vector a2(3);
    GSL::Vector a3(3);

    a1 = this->lattice[0]*this->scale;
    a2 = this->lattice[1]*this->scale;
    a3 = this->lattice[2]*this->scale;

    std::vector<GSL::Vector> possibles;
    std::unordered_set<GSL::Vector> tmp;
    possibles = {a1, a2, a3};
    size_t n = tmp.size(), prev = 0;
    std::cout << "0%| ";
    while(n < num){
        n = tmp.size();
        std::sort(possibles.begin(), possibles.end(), comp_norm);
        tmp.insert(possibles[0]);
        for(std::unordered_set<GSL::Vector>::const_iterator it = tmp.begin(); it != tmp.end(); it++){
            /* Only consider vectors not already found */
            if(tmp.find(possibles[0] + *it) == tmp.end()){
                possibles.push_back(possibles[0] + *it);
            }
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
    for(GSL::Vector L : tmp){
        Rn_vecs.push_back(L);
        Rn_vecs.push_back(L - 2*(GSL::dot(L, a2))*a2/(a2.norm()*a2.norm()));
        Rn_vecs.push_back(L - 2*(GSL::dot(L, a3))*a3/(a3.norm()*a3.norm()));
        Rn_vecs.push_back(L - 2*(GSL::dot(L, a1))*a1/(a1.norm()*a1.norm()));
        Rn_vecs.push_back(L - 2*(GSL::dot(L, a2+a1))*(a2+a1)/((a2+a1).norm()*(a2+a1).norm()));
        Rn_vecs.push_back(L - 2*(GSL::dot(L, a2+a3))*(a2+a3)/((a2+a3).norm()*(a2+a3).norm()));
        Rn_vecs.push_back(L - 2*(GSL::dot(L, a1+a3))*(a1+a3)/((a1+a3).norm()*(a1+a3).norm()));
        Rn_vecs.push_back(L - 2*(GSL::dot(L, a1+a2+a3))*(a1+a2+a3)/((a1+a2+a3).norm()*(a1+a2+a3).norm()));
    }
    std::sort(Rn_vecs.begin(), Rn_vecs.end(), comp_norm);
    std::cout << "Number of lattice vectors : " << Rn_vecs.size() << std::endl;
}

void Crystal::calc_Kn(size_t num)
{
    std::cout << "Calculating reciprocal lattice vectors" << std::endl;
    GSL::Vector a1(3), b1(3);
    GSL::Vector a2(3), b2(3);
    GSL::Vector a3(3), b3(3);

    a1 = this->lattice[0]*this->scale;
    a2 = this->lattice[1]*this->scale;
    a3 = this->lattice[2]*this->scale;

    b1 = 2*M_PI*GSL::cross(a2, a3)/volume;
    b2 = 2*M_PI*GSL::cross(a3, a1)/volume;
    b3 = 2*M_PI*GSL::cross(a1, a2)/volume;

    std::vector<GSL::Vector> possibles;
    std::unordered_set<GSL::Vector> tmp;
    possibles = {b1, b2, b3};
    size_t n = 0, prev = 0;
    std::cout << "0%| ";
    while(n < num){
        n = tmp.size();
        std::sort(possibles.begin(), possibles.end(), comp_norm);
        tmp.insert(possibles[0]);
        for(std::unordered_set<GSL::Vector>::const_iterator it = tmp.begin(); it != tmp.end(); it++){
            if(tmp.find(possibles[0] + *it) == tmp.end()){
                possibles.push_back(possibles[0] + *it);
            }
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
    for(GSL::Vector G : tmp){
        Kn_vecs.push_back(G);
        Kn_vecs.push_back(G - 2*(GSL::dot(G, b2))*b2/(b2.norm()*b2.norm()));
        Kn_vecs.push_back(G - 2*(GSL::dot(G, b3))*b3/(b3.norm()*b3.norm()));
        Kn_vecs.push_back(G - 2*(GSL::dot(G, b1))*b1/(b1.norm()*b1.norm()));
        Kn_vecs.push_back(G - 2*(GSL::dot(G, b2+b1))*(b2+b1)/((b2+b1).norm()*(b2+b1).norm()));
        Kn_vecs.push_back(G - 2*(GSL::dot(G, b2+b3))*(b2+b3)/((b2+b3).norm()*(b2+b3).norm()));
        Kn_vecs.push_back(G - 2*(GSL::dot(G, b1+b3))*(b1+b3)/((b1+b3).norm()*(b1+b3).norm()));
        Kn_vecs.push_back(G - 2*(GSL::dot(G, b1+a2+b3))*(b1+b2+b3)/((b1+b2+b3).norm()*(b1+b2+b3).norm()));
    }
    std::sort(Kn_vecs.begin(), Kn_vecs.end(), comp_norm);
    std::cout << "Number of reciprocal lattice vectors : " << Kn_vecs.size() << std::endl;
}

double Crystal::calc_volume()
{
    GSL::Vector tmp;
    tmp = GSL::cross(lattice[1], lattice[2]);
    return std::abs(GSL::dot(lattice[0], tmp)*GSL::pow_int(scale, 3));
}

Crystal::Crystal()
 : Rn_vecs(), Kn_vecs(), scale(), lattice(), atoms(), volume()
{
}

Crystal::Crystal(double& a)
 : Rn_vecs(), Kn_vecs(), scale(a), lattice(3,3), atoms()
{
    lattice[0][0] = 1.0;
    lattice[1][1] = 1.0;
    lattice[2][2] = 1.0;

    volume = calc_volume();

    calc_Rn(15);
    calc_Kn(15);
}

Crystal::Crystal(const double& a)
 : Rn_vecs(), Kn_vecs(), scale(a), lattice(3,3), atoms()
{
    lattice[0][0] = 1.0;
    lattice[1][1] = 1.0;
    lattice[2][2] = 1.0;

    volume = calc_volume();

    calc_Rn(15);
    calc_Kn(15);
}

Crystal::Crystal(const double& a, const double& b, const double& c)
 : Rn_vecs(), Kn_vecs(), scale(a), lattice(3,3), atoms()
{
    lattice[0][0] = 1.0;
    lattice[1][1] = b/a;
    lattice[2][2] = c/a;

    volume = calc_volume();

    calc_Rn(15);
    calc_Kn(15);
}

Crystal::Crystal(double& a, double& b, double& c)
 : Rn_vecs(), Kn_vecs(), scale(a), lattice(3,3), atoms()
{
    lattice[0][0] = 1.0;
    lattice[1][1] = b/a;
    lattice[2][2] = c/a;

    volume = calc_volume();

    calc_Rn(15);
    calc_Kn(15);
}

Crystal::Crystal(const GSL::Vector& a, const GSL::Vector& b, const GSL::Vector& c)
 : Rn_vecs(), Kn_vecs(), scale(), lattice(3,3), atoms()
{
    scale = a.norm();
    lattice[0] = (1./scale)*a;
    lattice[1] = (1./scale)*b;
    lattice[2] = (1./scale)*c;

    volume = calc_volume();

    calc_Rn(15);
    calc_Kn(15);
}

Crystal::Crystal(GSL::Vector& a, GSL::Vector& b, GSL::Vector& c)
 : Rn_vecs(), Kn_vecs(), scale(), lattice(3,3), atoms()
{
    scale = a.norm();
    lattice[0] = a/scale;
    lattice[1] = b/scale;
    lattice[2] = c/scale;

    volume = calc_volume();

    calc_Rn(15);
    calc_Kn(15);

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
