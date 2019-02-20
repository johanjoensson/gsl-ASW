#include "atomic_quantity.h"
#include "utils.h"
#include "GSLpp/basic_math.h"

Atomic_quantity::Atomic_quantity(const std::vector<Atom>& atoms)
 : sites(atoms), val(atoms.size(), std::vector<double>(0,0))
{
    for(size_t i = 0; i < sites.size(); i++){
        val[i] = std::vector<double>(sites[i].mesh.r.size(), 0.);
    }
}

double Atomic_quantity::operator()(const GSL::Vector& r)
{
    double res = 0.;
    GSL::Vector ri(3);
    Atom at;
    size_t t = 1;
    for(size_t i = 0; i < sites.size(); i++){
        at = sites[i];
        ri = r - at.get_pos();
        if(ri.norm<double>() <= at.get_AS()){
            t = 1;

            // Linear interpolation to find value at point ri
            while(ri.norm<double>() > at.mesh.r[t] && t < at.mesh.r.size()){
                t++;
            }
            if(t < at.mesh.r.size()){
                res += lerp(ri.norm<double>(), at.mesh.r[t-1], at.mesh.r[t], val[i][t-1], val[i][t]);
            }
        }
    }

    return res;
}

Potential::Potential(std::vector<Atom>& atoms,
    std::function<double(const size_t, const double)> atomic_potential)
 : Atomic_quantity(atoms), electrostatic(atoms.size()), exchange_correlation(atoms.size()),
  at_pot(atomic_potential), xc_fun(), MT_0(0)
{
    for(size_t i = 0; i < sites.size(); i++){
        electrostatic[i] = std::vector<double>(sites[i].mesh.r.size(), 0.);
        exchange_correlation[i] = std::vector<double>(sites[i].mesh.r.size(), 0.);
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

double Xi0(size_t j, std::vector<Atom> sites, double r)
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

void Potential::initial_pot(size_t nel, double vol)
{
    std::vector<double> rho;
    for(size_t i = 0; i < sites.size(); i++){
        rho = std::vector<double>(sites[i].mesh.r.size(), static_cast<double>(nel)/vol);
        exchange_correlation[i] = xc_fun.exc(rho);
        for(size_t j = 0; j < sites[i].mesh.r.size(); j++){
            electrostatic[i][j] =
            /*mod_at_pot(5*sites[i].get_Z(), sites[i].mesh.r[j]);*/
            at_pot(sites[i].get_Z(), sites[i].mesh.r[j]);
        }
    }
    std::vector<double> alpha;
    double r = 0., r1 = 0., r2 = 0., drx = 0., drx1 = 0., drx2 = 0.;
    for(size_t j = 0; j < sites.size(); j++){
        alpha = std::vector<double>(sites[j].mesh.r.size(), 0.);
        alpha[0] = 0.;
        alpha[1] = 1./6 * (4*Xi0(j, sites, sites[j].mesh.r[0])*sites[j].mesh.drx[0]
    + (Xi0(j, sites, sites[j].mesh.r[1]) + Xi0(j, sites, -sites[j].mesh.r[1]))*sites[j].mesh.drx[1]);
        electrostatic[j][0] += Xi0(j, sites, 0);
        electrostatic[j][1] += alpha[1] * 1/sites[j].mesh.r[1];
        for(size_t rj = 2; rj < sites[j].mesh.r.size(); rj++){
            r = sites[j].mesh.r[rj];
            r1 = sites[j].mesh.r[rj-1];
            r2 = sites[j].mesh.r[rj-2];
            drx = sites[j].mesh.drx[rj];
            drx1 = sites[j].mesh.drx[rj-1];
            drx2 = sites[j].mesh.drx[rj-2];
            alpha[rj] = alpha[rj - 2];
            alpha[rj] += 1./6 * drx2 * (Xi0(j, sites, r2) + Xi0(j, sites, -r2));
            alpha[rj] += 4./6 * drx1 * (Xi0(j, sites, r1) + Xi0(j, sites, -r1));
            alpha[rj] += 1./6 * drx * (Xi0(j, sites, r) + Xi0(j, sites, -r));

           electrostatic[j][rj] += alpha[rj]/r;
        }
    }
    for(size_t i = 0; i < sites.size(); i++){
        for(size_t j = 0; j < sites[i].mesh.r.size(); j++){
            val[i][j] = electrostatic[i][j] + exchange_correlation[i][j];
        }
    }


    // Calculate MT_0 as average potential over all atomic spheres
    double MT_0_s = 0, areas = 0;
    for(size_t i = 0; i < sites.size(); i++){
        MT_0_s += val[i].back()*GSL::pow_int(sites[i].get_AS(), 2);
        areas += GSL::pow_int(sites[i].get_AS(), 2);
    }
    this->MT_0 = MT_0_s/areas;

	std::cout << "MT0 = " << this->MT_0 << std::endl;

    // Make potential relative to MT_0
    for(size_t i = 0; i < sites.size(); i++){
        for(size_t j = 0; j < sites[i].mesh.r.size(); j++){
            val[i][j] -= this->MT_0;
        }
    }
}

Density::Density(std::vector<Atom>& atoms)
 : Atomic_quantity(atoms), valence(atoms.size()), core(atoms.size())
{}
