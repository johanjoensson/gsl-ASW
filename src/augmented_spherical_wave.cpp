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
    for(size_t i = 0; i < off_centers.size(); i++){
        Atom at = off_centers[i];
        if(at != center){
            if(at.Z >= 20){
                l_low = 3;
            }else{
                l_low = 2;
            }
            for(int l_s = 0; l_s <= std::min(l_low + 1, n - 1); l_s++){
                for(int m_s = -l_s; m_s <= l_s; m_s++){
                    lm lp = {l_s, m_s};
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
	out_file.open("check_Hankel.dat", std::fstream::out|std::fstream::app);
	out_file << "# r\tV(r)\trH(r)" << "\n";
	for(size_t j = 0; j < H.val.size(); j++){
		out_file << std::setprecision(8) << center.mesh.r[j] << " " << v_tot[j] << " " << H.val[j] << "\n";
	}
    out_file << "\n\n";
	out_file.close();
#endif
    for(Atom at : off_centers){
        std::unordered_set<Augmented_Bessel> Ji_tmp;
        it = find(v.sites.begin(), v.sites.end(), at);
        i = static_cast<size_t>(std::distance(v.sites.begin(), it));
        for(Augmented_Bessel Jil : J[i]){
            v_tot = v_eff(Jil.mesh, v.val[i], Jil.l);
            Jil.update(v_tot, en, core_state);
            Ji_tmp.insert(Jil);
#ifdef DEBUG
            out_file.open("check_Bessel.dat", std::fstream::out|std::fstream::app);
            out_file << "# r\tV(r)\trJ(r)" << std::endl;
            for(size_t j = 0; j < Jil.val.size(); j++){
                out_file << std::setprecision(8) << at.mesh.r[j] << " " << v_tot[j] << " " << Jil.val[j] << "\n";
            }
            out_file << "\n\n";
            out_file.close();
#endif
        }
        J[i].swap(Ji_tmp);
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
