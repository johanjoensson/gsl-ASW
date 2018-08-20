#include "atomic_quantity.h"
#include "utils.h"
#include "../../GSL-lib/src/basic_math.h"

Atomic_quantity::Atomic_quantity(const std::vector<Atom> atoms)
 : sites(atoms), val(atoms.size())
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
        if(ri.norm() <= at.get_AS()){
            t = 1;

            // Linear interpolation to find value at point ri
            while(ri.norm() > at.mesh.r[t] && t < at.mesh.r.size()){
                t++;
            }
            if(t < at.mesh.r.size()){
                res += lerp(ri.norm(), at.mesh.r[t-1], at.mesh.r[t], val[i][t-1], val[i][t]);
            }
        }
    }

    return res;
}

Potential::Potential(std::vector<Atom>& atoms)
 : Atomic_quantity(atoms), electrostatic(atoms.size()), exchange_correlation(atoms.size())
{}

double atomic_potential(const int Z, const double r)
{
    return -2.*Z/r;
}


void Potential::initial_pot()
{
    for(size_t i = 0; i < sites.size(); i++){
        for(size_t j = 0; j < sites[i].mesh.r.size(); j++){
            electrostatic[i][j] = atomic_potential(sites[i].get_Z(), sites[i].mesh.r[j]);
            exchange_correlation[i][j] = 0.;
            val[i][j] = electrostatic[i][j] + exchange_correlation[i][j];
        }
    }
    std::vector<double> A;
    double xi0 = 0, a = 0., r = 0.;
    for(size_t i = 0; i < sites.size(); i++){
        A = std::vector<double>(sites[i].mesh.r.size(), 0.);
        for(size_t ri = 0; ri < sites[i].mesh.r.size(); ri++){
            r = sites[i].mesh.r[ri];
            for(size_t j = 0; j < sites.size(); j++){
                if(i != j){
                    
                    a = (sites[i].pos - sites[j].pos).norm();
                    xi0 = (a + r)/a*atomic_potential();
                    A[ri] = ;
                }
            }
        }
    }


    double MT_0 = 0, vols = 0;
    for(size_t i = 0; i < sites.size(); i++){
        // Add v[rs]*4*Pi*rs^2, the value of the potential integrated over the
        //entire atomic sphere, to V_MT
        MT_0 += val[i][sites[i].mesh.r.size() - 1]*4*M_PI*GSL::pow_int(2,sites[i].get_AS());
        vols += 4*M_PI*GSL::pow_int(2,sites[i].get_AS());
    }
    // Divide V_MT by the total atomic areas in the crystal, thus getting an
    // average of the potential on all spheres
    MT_0 /= vols;
}

Density::Density(std::vector<Atom>& atoms)
 : Atomic_quantity(atoms), valence(atoms.size()), core(atoms.size())
{}
