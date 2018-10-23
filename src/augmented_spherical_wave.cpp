#include "augmented_spherical_wave.h"
#include "spherical_fun.h"
#include "numerov_solver.h"
#include "structure_const.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>


Augmented_spherical_wave::Augmented_spherical_wave()
 : center(), off_centers(), kappa(),  n(), l(), s(), H(), J()
{}

Augmented_spherical_wave::Augmented_spherical_wave(double kappa, unsigned int n,
     lm l, spin s, Atom& center, std::vector<Atom>& off_centers)
 : center(center), off_centers(off_centers), kappa(kappa),
  n(n), l(l), s(s), H(n, l, kappa, center.pos, center.mesh), J(off_centers.size())
{
    // Set up off-center spheres
    int l_low = 2;
    if(center.Z >= 21){
        l_low = 3;
    }
    Numerov_solver sol;
    Augmented_Bessel tmp;
    std::vector<Atom>::iterator it = find(off_centers.begin(), off_centers.end(),
    center);
    size_t i;// = std::distance(off_centers.begin(), it);
    for(Atom at : off_centers){
        it = find(off_centers.begin(), off_centers.end(), at);
        i = std::distance(off_centers.begin(), it);
        if(at != center){
            for(int l = 0; l <= std::min((int)n - 1, l_low + 1); l++){
                for(int m = -l; m <= l; m++){
                    lm lp = {l, m};
                    tmp = Augmented_Bessel(n, lp, kappa, at.pos, at.mesh);
                    J[i].insert(tmp);
                }
            }
        }
    }
}

std::vector<double> v_eff(const Logarithmic_mesh &mesh,
    const std::vector<double>& v, const lm& l)
{
    std::vector<double> res(v.size(), 0.);
    for(size_t i = 0; i < v.size(); i++){
        res[i] = v[i] + l.l*(l.l + 1)/mesh.r2[i];
    }
    return res;
}

void Augmented_spherical_wave::set_up(Potential &v)
{
    double en = -1.*center.get_Z()*center.get_Z()/(n*n) + v.MT_0;

    // On-center expression
    std::vector<Atom>::iterator it = find(v.sites.begin(), v.sites.end(),
    center);
    size_t i = std::distance(v.sites.begin(), it);
    std::vector<double> v_tot = v_eff(center.mesh, v.val[i], l);

    H.update(v_tot, en, core_state);

#ifdef DEBUG
	std::ofstream out_file;
	out_file.open("check_Hankel.dat");
	out_file << "# r\trH" << std::endl;
	for(size_t i = 0; i < H.val.size(); i++){
		out_file << std::setprecision(8) << center.mesh.r[i] << " " << v_tot[i] << " " << H.val[i] << std::endl;
	}
	out_file.close();
#endif
    std::unordered_set<Augmented_Bessel> tmp;
    for(Atom at : off_centers){
        it = find(v.sites.begin(), v.sites.end(), at);
        i = std::distance(v.sites.begin(), it);

        if(at != center){
            en =H.EH;

            for(Augmented_Bessel Jil : J[i]){
                v_tot = v_eff(Jil.mesh, v.val[i], Jil.l);
                Jil.update(v_tot, en, core_state);
                tmp.insert(Jil);
            }
            J[i].swap(tmp);
            tmp.clear();
        }
    }
}

double Augmented_spherical_wave::operator()(const GSL::Vector &r)
{
    double res = 0.;
    // Start with on-center value
    res += H(r)* cubic_harmonic(l, r).val;

    // Then do off-center contributions
    GSL::Vector R(3);
    Structure_constant B;
    Atom at;
    for(size_t j = 0; j < off_centers.size(); j++){
        at = off_centers[j];
        if(at.pos != center.pos){
            for(Augmented_Bessel Jj : J[j]){
                R = Jj.center - center.pos;
                B = Structure_constant( Jj.l, this->l, R);
                res += Jj(r)*cubic_harmonic(Jj.l, r - Jj.center).val*B.val;
            }
        }
    }

    return res;
}
