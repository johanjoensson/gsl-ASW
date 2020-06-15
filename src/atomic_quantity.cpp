#include "atomic_quantity.h"
#include "utils.h"
#include "GSLpp/basic_math.h"

Atomic_quantity::Atomic_quantity(const Crystal_t<3, Atom>& cr, std::vector<Logarithmic_mesh>& at_meshes_n)
 : atoms(cr.atoms()), at_meshes(at_meshes_n), val(cr.atoms().size(), std::vector<double>(0,0))
{
    for(size_t at_i = 0; at_i < cr.atoms().size(); at_i++){
        val[at_i] = std::vector<double>(at_meshes[at_i].size(), 0.);
    }
}

double Atomic_quantity::operator()(const GSL::Vector& r)
{
    double res = 0.;
    GSL::Vector ri(3);
    size_t t = 1;
    for(size_t i = 0; i < atoms.size(); i++){
        Atom at{atoms[i]};
        ri = r - at.get_pos();
        if(ri.norm<double>() <= at.get_AS()){
            t = 1;

            // Linear interpolation to find value at point ri
            while(ri.norm<double>() > at_meshes[i].r(t) && t < at_meshes[i].size()){
                t++;
            }
            if(t < at_meshes[i].size()){
                res += lerp(ri.norm<double>(), at_meshes[i].r(t-1), at_meshes[i].r(t), val[i][t-1], val[i][t]);
            }
        }
    }

    return res;
}

Potential::Potential(const Crystal_t<3, Atom>& crystal, std::vector<Logarithmic_mesh>& at_meshes_n,
    std::function<double(const size_t, const double)> atomic_potential)
 : Atomic_quantity(crystal, at_meshes_n), electrostatic(crystal.atoms().size()), exchange_correlation(crystal.atoms().size()),
  at_pot(atomic_potential), MT_0(), xc_fun()
{
    for(size_t i = 0; i < atoms.size(); i++){
        electrostatic[i] = std::vector<double>(at_meshes[i].size(), 0.);
        exchange_correlation[i] = std::vector<double>(at_meshes[i].size(), 0.);
    }
}

void Potential::set_xc_fun(XC_FUN xc_func)
{
    xc_fun.set_xc(xc_func);
}



double atomic_potential(const size_t Z, const double r)
{
    return -2.*static_cast<double>(Z)/r;
}

double mod_at_pot(const size_t Z, const double r)
{
    double alpha = 0.53625;
    double x = r*std::cbrt(Z)/0.88534;
    double phi = 1./((1 + alpha*x)*(1 + alpha*x));
    return -(1 + static_cast<double>(Z - 1)*phi)/r;
}

double Xi0(size_t j, const std::vector<Atom>& sites, double r)
{
    double res = 0, a = 0;
    for(size_t i = 0; i < sites.size(); i++){
        if(i != j){
            a = (sites[j].pos - sites[i].pos).norm<double>();
            res += (a + r)/a*atomic_potential(sites[i].get_Z(), a + r);
        }
    }
    return res;
}

void Potential::initial_pot(double vol)
{
    size_t nel = 0;
    for(const auto& at : atoms){
        nel += at.get_Z();
    }
    std::vector<double> rho;
    for(size_t i = 0; i < atoms.size(); i++){
        rho = std::vector<double>(at_meshes[i].size(), static_cast<double>(nel)/vol);
        exchange_correlation[i] = xc_fun.exc(rho);
        for(size_t j = 0; j < at_meshes[i].size(); j++){
            electrostatic[i][j] =
            /*mod_at_pot(5*sites[i].get_Z(), sites[i].mesh.r[j]);*/
            at_pot(atoms[i].get_Z(), at_meshes[i].r(j));
        }
    }
    std::vector<double> alpha;
    double r = 0., r1 = 0., r2 = 0., drx = 0., drx1 = 0., drx2 = 0.;
    for(size_t j = 0; j < atoms.size(); j++){
        alpha = std::vector<double>(at_meshes[j].size(), 0.);
        alpha[0] = 0.;
        alpha[1] = 1./6 * (4*Xi0(j, atoms, at_meshes[j].r(0))*at_meshes[j].drx(0)
    + (Xi0(j, atoms, at_meshes[j].r(1)) + Xi0(j, atoms, -at_meshes[j].r(1)))*at_meshes[j].drx(1));
        electrostatic[j][0] += Xi0(j, atoms, 0);
        electrostatic[j][1] += alpha[1] * 1/at_meshes[j].r(1);
        for(size_t rj = 2; rj < at_meshes[j].size(); rj++){
            r = at_meshes[j].r(rj);
            r1 = at_meshes[j].r(rj-1);
            r2 = at_meshes[j].r(rj-2);
            drx = at_meshes[j].drx(rj);
            drx1 = at_meshes[j].drx(rj-1);
            drx2 = at_meshes[j].drx(rj-2);
            alpha[rj] = alpha[rj - 2];
            alpha[rj] += 1./6 * drx2 * (Xi0(j, atoms, r2) + Xi0(j, atoms, -r2));
            alpha[rj] += 4./6 * drx1 * (Xi0(j, atoms, r1) + Xi0(j, atoms, -r1));
            alpha[rj] += 1./6 * drx * (Xi0(j, atoms, r) + Xi0(j, atoms, -r));

            electrostatic[j][rj] += alpha[rj]/r;
        }
    }

    for(size_t i = 0; i < atoms.size(); i++){
        for(size_t j = 0; j < at_meshes[i].size(); j++){
            val[i][j] = electrostatic[i][j]  + exchange_correlation[i][j] ;
        }
    }


    // Calculate MT_0 as average potential over all atomic spheres
    double MT_0_s = 0, areas = 0;
    for(size_t i = 0; i < atoms.size(); i++){
        std::cout << "R-AS[" << i << "] = " << atoms[i].get_AS() <<"\n";
        MT_0_s += val[i].back()*GSL::pow_int(atoms[i].get_AS(), 2);
        areas += GSL::pow_int(atoms[i].get_AS(), 2);
    }
    this->MT_0 = MT_0_s/areas;

	std::cout << "MT0 = " << this->MT_0 << std::endl;

    // Make potential relative to MT_0
    // for(size_t i = 0; i < sites.size(); i++){
    //     for(size_t j = 0; j < sites[i].mesh.size(); j++){
    //         val[i][j] -= this->MT_0;
    //     }
    // }
}

Density::Density(const Crystal_t<3, Atom>& cr, std::vector<Logarithmic_mesh>& at_meshes_n)
 : Atomic_quantity(cr, at_meshes_n), valence(cr.atoms().size()), core(cr.atoms().size())
{}
