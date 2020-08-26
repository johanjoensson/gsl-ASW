#include "augmented_spherical_wave.h"
#include "spherical_fun.h"
#include "numerov_solver.h"
#include "structure_const.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>

/*
Augmented_spherical_wave::Augmented_spherical_wave(std::vector<Hankel_container>& Hs_n, std::vector<Bessel_container>& Bs_n, Site_t<3> center_n, double kappa_n, lm l_n, spin s_n)
 : Hs_m(Hs_n), Bs_m(Bs_n), center(center_n), kappa(kappa_n), l(l_n), s(s_n)
{}
*/
/*
double Augmented_spherical_wave::operator()(const GSL::Vector &r) const
{
    double res = 0.;
    // Start with on-center value
    res += Hs_m[center.index()](l, kappa, s)(r)* cubic_harmonic(l, r).val;

    // Then do off-center contributions
    GSL::Vector R(3);
    for(auto Js : Bs_m){
        for(Augmented_Bessel Jj : Js){
            R =  - center.pos();
            Structure_constant B = Structure_constant(Js.back().l.n - 1, Jj.l, this->l);
            res += Jj(r)*cubic_harmonic(Jj.l, r).val*B(R);
        }
    }

    return res;
}
*/
