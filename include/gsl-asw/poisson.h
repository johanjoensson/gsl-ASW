#ifndef POISSON_EQ_H
#define POISSON_EQ_H
#include "numerov_solver.h"
#include "numerical-mesh/numerical-mesh.h"
#include "numerical-mesh/numerical-mesh-integration.h"
#include <vector>


class Poisson_equation{
public:
    // const Mesh_base<1, double>& mesh_m;
    std::vector<double> val_m;
    // int l_m;

    void solve(const Mesh_base<1, double>& mesh, const std::vector<double>& rho, const int l)
    {
        val_m.resize(mesh.size());
        std::vector<double> s(mesh.size(), 0), g(mesh.size(), 0);
        auto s_i = s.begin(), g_i = g.begin();
        auto rho_i = rho.begin();
        for(auto m_i = mesh.begin(); m_i != mesh.end(); m_i++){
            *(s_i++) = -4.*M_PI*2*m_i->r()*std::sqrt(GSL::pow_int(m_i->dr(), 3))*(*(rho_i++));
            *g_i =  1./m_i->dr()*(0.5*m_i->d3r() - 0.75*GSL::pow_int(m_i->d2r(), 2)/m_i->dr());
            if (l != 0){
                *g_i -= l*(l + 1)/m_i->r2()*GSL::pow_int(m_i->dr(), 2);
            }
            g_i++;
        }

        std::vector<double> init;
        /***********************************************************************
        * Shift of energy to match boundary condition VH(r(0)) = V0            *
        ***********************************************************************/
        double v0 = 0;
        if (l == 0){
            /*******************************************************************
            * Initial value, VH(0) = V0
            *******************************************************************/
            v0 = 8.*M_PI*simpson_integral<double>(mesh,
                [&](const Mesh_base<1, double>::mesh_point& p)
                {
                    return p.r()*rho[static_cast<size_t>(p.i())];
                }
            );
            auto r0 = mesh.r(0);
            auto dr0 = mesh.dr(0);
            auto F0p = 0/*v0*std::sqrt(dr0) - 0.5*r0*v0*GSL::pow_int(std::sqrt(dr0), -3)*/;
            init = std::vector<double> {
                /***************************************************************
                * Solve with boundary condition VH(r(0)) = VH'(r(0)) = 0       *
                ***************************************************************/
                0,
                /***************************************************************
                * O(h^5) estimate of r(x)V(r(x))|x=1, from
                * L. M. Quiroz González, and D. Thompson
                * Computers in Physics 11, 514 (1997); doi: 10.1063/1.168593
                ***************************************************************/
                (F0p*(1+g[2]/12) + 1./24*(7*s[0] + 6*s[1] - s[2]) + g[2]/36*(s[0] + 2*s[1]))/
                    (1 + g[1]/4 + g[1]*g[2]/18)
            };

        }else{
            init = std::vector<double> {
                GSL::pow_int(mesh.r(0), l + 1)/std::sqrt(mesh.dr(0)),
                GSL::pow_int(mesh.r(1), l + 1)/std::sqrt(mesh.dr(1)),
                GSL::pow_int(mesh.r(2), l + 1)/std::sqrt(mesh.dr(2))
            };
        }
        std::vector<double> empty(0);
        Numerov_solver sol;
        sol.solve(val_m.begin(), val_m.end(), g.begin(), g.end(), s.begin(),
            s.end(), init.begin(), init.end(), empty.begin(), empty.end(),
            val_m.end());

        auto S = mesh.r(mesh.size() - 1);
        double right_val = 8.*M_PI/(2*l + 1)*GSL::pow_int(1./S, l)*
            simpson_integral<double>(mesh,
                [&](const Mesh_base<1, double>::mesh_point& p)
                {
                    return GSL::pow_int(p.r(), l + 2)*rho[static_cast<size_t>(p.i())];
                }
            );

        double Al = (right_val - val_m.back()*std::sqrt(mesh.dr().back()) - S*v0)/GSL::pow_int(S, l + 1);
        auto v_i = val_m.begin();
        /***********************************************************************
        * Transformations : rVH(r)/sqrt(dr) -> rVH(r) -> rVH(r) + Alr^(l + 1) ->
        * -> VH(r) + Alr^l -> VH(r) + Alr^l + V0
        ***********************************************************************/
        for(auto m_i = mesh.begin(); m_i != mesh.end(); m_i++, v_i++){
            *v_i = (*v_i*std::sqrt(m_i->dr())
                    + Al*GSL::pow_int(m_i->r(), l + 1))/m_i->r()
                    + v0;
        }
        /***********************************************************************
        * In case r(0) = 0, set VH(0) = V0 explicitly                          *
        ***********************************************************************/
        val_m.front() = v0;
    }

    Poisson_equation()
     : val_m()
    {
        // solve();
    }



    std::vector<double> val() const noexcept { return this->val_m;}
    double val(const size_t i) const noexcept { return this->val_m[i];}
    // std::vector<double> rho() const noexcept { return this->rho_m;}
    // double rho(const size_t i) const noexcept { return this->rho_m[i];}

};

#endif // POISSON_EQ_H
