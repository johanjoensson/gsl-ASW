#include "atomic_quantity.h"
#include "utils.h"
#include "GSLpp/basic_math.h"
#include <iomanip>

Atomic_quantity::Atomic_quantity(const Crystal_t<3, Atom>& cr_n, const std::vector<Logarithmic_mesh>& at_meshes_n)
 : cr(cr_n), at_meshes(at_meshes_n), val(cr_n.atoms().size(), std::vector<double>(0,0))
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
    for(const auto& site : cr.sites()){
        size_t at_i = cr.atom_index(site);
        Atom at{cr.atom(site)};
        ri = r - site.pos();
        if(ri.norm<double>() <= at.get_AS()){
            t = 1;

            // Linear interpolation to find value at point ri
            while(ri.norm<double>() > at_meshes[at_i].r(t) && t < at_meshes[at_i].size()){
                t++;
            }
            if(t < at_meshes[at_i].size()){
                res += lerp(ri.norm<double>(), at_meshes[at_i].r(t-1), at_meshes[at_i].r(t), val[at_i][t-1], val[at_i][t]);
            }
        }
    }

    return res;
}

Potential::Potential(const Crystal_t<3, Atom>& crystal, const std::vector<Logarithmic_mesh>& at_meshes_n,
    std::function<double(const size_t, const double)> atomic_potential)
 : Atomic_quantity(crystal, at_meshes_n), electrostatic(crystal.atoms().size()), exchange_correlation(crystal.atoms().size()),
  at_pot(atomic_potential), MT_0(), xc_fun()
{}

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

double Potential::Xi0(const Site_t<3>& j, const double r)
{
    double res = 0, a = 0;
    for(auto site : cr.sites()){
        if(site != j){
            a = (j.pos() - site.pos()).norm<double>();
            res += (a + r)/a*atomic_potential(cr.atom(site).get_Z(), a + r);
        }
    }
    return res;
}

void Potential::initial_pot(double vol)
{
    size_t nel = 0;
    for(const auto& at : cr.atoms()){
        nel += at.get_Z();
    }
    std::vector<double> rho;
    for(size_t i = 0; i < cr.atoms().size(); i++){
        electrostatic[i] = std::vector<double>(at_meshes[i].size());
        val[i] = std::vector<double>(at_meshes[i].size());
        // exchange_correlation[i] = std::vector<double>(at_meshes[i].size(), 0.);
        rho = std::vector<double>(at_meshes[i].size(), static_cast<double>(nel)/vol);
        exchange_correlation[i] = xc_fun.exc(rho);
        for(size_t j = 0; j < at_meshes[i].size(); j++){
            electrostatic[i][j] =
            /*mod_at_pot(5*sites[i].get_Z(), sites[i].mesh.r[j]);*/
            at_pot(cr.atoms(i).get_Z(), at_meshes[i].r(j));
        }
    }
    std::vector<bool> done(cr.atoms().size(), false);
    double r = 0., r1 = 0., r2 = 0., drx = 0., drx1 = 0., drx2 = 0.;
    for(const auto& site : cr.sites()){
        size_t at_i = cr.atom_index(site);
        if(done[at_i]){
            continue;
        }
        std::vector<double> alpha(at_meshes[at_i].size());
        alpha[0] = 0.;
        alpha[1] = 1./6 * (4*Xi0(site, at_meshes[at_i].r(0))*at_meshes[at_i].drx(0)
    + (Xi0(site, at_meshes[at_i].r(1)) + Xi0(site, -at_meshes[at_i].r(1)))*at_meshes[at_i].drx(1));
        electrostatic[at_i][0] += Xi0(site, 0);
        electrostatic[at_i][1] += alpha[1] * 1/at_meshes[at_i].r(1);
        for(size_t ri = 2; ri < at_meshes[at_i].size(); ri++){
            r = at_meshes[at_i].r(ri);
            r1 = at_meshes[at_i].r(ri-1);
            r2 = at_meshes[at_i].r(ri-2);
            drx = at_meshes[at_i].drx(ri);
            drx1 = at_meshes[at_i].drx(ri-1);
            drx2 = at_meshes[at_i].drx(ri-2);
            alpha[ri] = alpha[ri - 2];
            alpha[ri] += 1./6 * drx2 * (Xi0(site, r2) + Xi0(site, -r2));
            alpha[ri] += 4./6 * drx1 * (Xi0(site, r1) + Xi0(site, -r1));
            alpha[ri] += 1./6 * drx * (Xi0(site, r) + Xi0(site, -r));

            // electrostatic[at_i][ri] += alpha[ri]/r;
        }
        done[at_i] = true;
        std::cout << "Done setting up electrostatic potential of atom " << at_i << "\n";
    }

    for(size_t i = 0; i < cr.atoms().size(); i++){
        for(size_t j = 0; j < at_meshes[i].size(); j++){
            val[i][j] = electrostatic[i][j]  + exchange_correlation[i][j] ;
        }
    }


    // Calculate MT_0 as average potential over all atomic spheres
    double MT_0_s = 0, areas = 0;
    for(size_t i = 0; i < cr.atoms().size(); i++){
        MT_0_s += val[i].back()*GSL::pow_int(cr.atoms(i).get_AS(), 2);
        areas += GSL::pow_int(cr.atoms(i).get_AS(), 2);
        std::cout << std::fixed << std::setprecision(12) << "MT_0 contribution = " << val[i].back()*GSL::pow_int(cr.atoms(i).get_AS(), 2) << "\n";
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

Density::Density(const Crystal_t<3, Atom>& cr_n, std::vector<Logarithmic_mesh>& at_meshes_n)
 : Atomic_quantity(cr_n, at_meshes_n), valence(cr_n.atoms().size()), core(cr_n.atoms().size())
{}
