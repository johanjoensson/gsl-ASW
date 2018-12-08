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

Augmented_spherical_wave::Augmented_spherical_wave(double kappa_n, int n_n,
     lm l_n, spin s_n, const Atom& center_n, const std::vector<Atom>& off_centers_n)
 : center(center_n), off_centers(off_centers_n), kappa(kappa_n),
   n(n_n), l(l_n), s(s_n), H(n_n, l_n, kappa_n, center_n.pos, center_n.mesh),
   J(off_centers_n.size())
{
    // Set up off-center spheres
    int l_low = 2;
    Numerov_solver sol;
    Augmented_Bessel tmp;
    std::vector<Atom>::iterator it = find(off_centers.begin(), off_centers.end(),
    center);
    size_t i;// = std::distance(off_centers.begin(), it);
    for(Atom at : off_centers){
        it = find(off_centers.begin(), off_centers.end(), at);
        i = static_cast<size_t>(std::distance(off_centers.begin(), it));
        if(at != center){
            if(at.Z >= 20){
                l_low = 3;
            }else{
                l_low = 2;
            }
            for(int l_s = 0; l_s <= l_low + 1; l_s++){
                for(int m = -l_s; m <= l_s; m++){
                    lm lp = {l_s, m};
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
    double en = -1.*static_cast<double>(center.get_Z()*center.get_Z())/
                    static_cast<double>(n*n) + v.MT_0;

    // On-center expression
    std::vector<Atom>::iterator it = find(v.sites.begin(), v.sites.end(),
    center);
    size_t i = static_cast<size_t>(std::distance(v.sites.begin(), it));
    std::vector<double> v_tot = v_eff(center.mesh, v.val[i], l);

    H.update(v_tot, en, core_state);

#ifdef DEBUG
	std::ofstream out_file;
	out_file.open("check_Hankel.dat");
	out_file << "# r\trH" << std::endl;
	for(size_t j = 0; j < H.val.size(); j++){
		out_file << std::setprecision(8) << center.mesh.r[j] << " " << v_tot[j] << " " << H.val[j] << std::endl;
	}
	out_file.close();
#endif
    std::unordered_set<Augmented_Bessel> tmp;
    for(Atom at : off_centers){
        it = find(v.sites.begin(), v.sites.end(), at);
        i = static_cast<size_t>(std::distance(v.sites.begin(), it));

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
    int l_low = 2;
    for(size_t j = 0; j < off_centers.size(); j++){
        at = off_centers[j];
        if(at.pos != center.pos){
            for(Augmented_Bessel Jj : J[j]){
                R = Jj.center - center.pos;
                if(at.Z > 20){
                    l_low = 3;
                }
                B = Structure_constant( l_low, Jj.l, this->l, R);
                res += Jj(r)*cubic_harmonic(Jj.l, r - Jj.center).val*B.val;
            }
        }
    }

    return res;
}
