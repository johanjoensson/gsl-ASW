#ifndef POISSON_EQ_H
#define POISSON_EQ_H
#include "numerov_solver.h"
#include "numerical-mesh.h"
#include "numerical-mesh-integration.h"
#include <vector>


class Poisson_equation{
protected:
    const Mesh_base<1, double>& mesh_m;
    std::vector<double> rho_m, val_m;
    int l_m;

public:
    Poisson_equation(const Mesh_base<1, double>& mesh, const std::vector<double>& rho, const int l)
     : mesh_m(mesh), rho_m(rho), val_m(mesh.size()), l_m(l)
    {}

    void solve()
    {
        std::vector<double> s(mesh_m.size(), 0), g(mesh_m.size(), 0);
        auto s_i = s.begin(), g_i = g.begin(), rho_i = rho_m.begin();
        for(auto m_i = mesh_m.begin(); m_i != mesh_m.end(); m_i++){
            *(s_i++) = -8.*M_PI*m_i->r()*std::sqrt(GSL::pow_int(m_i->dr(), 3))*(*(rho_i++));
            *g_i =  1./m_i->dr()*(0.5*m_i->d3r() - 0.75*GSL::pow_int(m_i->d2r(), 2)/m_i->dr());
            if (l_m != 0){
                *g_i -= l_m*(l_m + 1)/m_i->r2()*GSL::pow_int(m_i->dr(), 2);
            }
            g_i++;
        }

        std::vector<double> init;
        double v0 = 0;
        if (l_m == 0){
            v0 = (8.*M_PI*simpson_integral<double>(mesh_m,
                [this](const double& i)
                {
                    return mesh_m.r(i)*rho_m[static_cast<size_t>(i)];
                }
            ));
            init = std::vector<double> {
                0,
                // v0*mesh_m.dr(0)/std::sqrt(mesh_m.dr(1))
                (v0*(1+g[2]/12) + 1./24*(6*s[1] - s[2]) + g[2]/36*2*s[1])/
                    (1 + g[1]/4 + g[1]*g[2]/18)
            };

        }else{
            init = std::vector<double> {
                GSL::pow_int(mesh_m.r(0), l_m + 1)/std::sqrt(mesh_m.dr(0)),
                GSL::pow_int(mesh_m.r(1), l_m + 1)/std::sqrt(mesh_m.dr(1)),
                GSL::pow_int(mesh_m.r(2), l_m + 1)/std::sqrt(mesh_m.dr(2))
            };
        }
        std::vector<double> empty(0);
        Numerov_solver sol;
        sol.solve(val_m.begin(), val_m.end(), g.begin(), g.end(), s.begin(),
            s.end(), init.begin(), init.end(), empty.begin(), empty.end(),
            val_m.end());

        auto v_i = val_m.begin();
        for(auto m_i = mesh_m.begin(); m_i != mesh_m.end(); m_i++, v_i++){
            *v_i *= std::sqrt(m_i->dr());
        }

        auto S = mesh_m.r(mesh_m.size() - 1);
        std::cout << "S = " << S << "\n";
        double right_val = 8.*M_PI/(2*l_m + 1)*GSL::pow_int(1./S, l_m)*
            simpson_integral<double>(mesh_m,
                [this](const double& i)
                {
                    return GSL::pow_int(mesh_m.r(i), l_m + 2)*rho_m[static_cast<size_t>(i)];
                }
            );

        std::cout << "Right boundary value = " << right_val << "\n";
        std::cout << "difference at right boundary " << std::abs(right_val - val_m.back()) << "\n";
        double Al = (right_val - val_m.back())/GSL::pow_int(S, l_m + 1);
        std::cout << "Al = " << Al << "\n";
        v_i = val_m.begin();
        for(auto m_i = mesh_m.begin(); m_i != mesh_m.end(); m_i++, v_i++){
             *v_i += Al*GSL::pow_int(m_i->r(), l_m + 1);
             *v_i /= m_i->r();
        }
        val_m.front() = v0;

    }

    std::vector<double> val() const { return this->val_m;}
    double val(const size_t i) const { return this->val_m[i];}
    std::vector<double> rho() const { return this->rho_m;}
    double rho(const size_t i) const { return this->rho_m[i];}

};

#endif // POISSON_EQ_H
