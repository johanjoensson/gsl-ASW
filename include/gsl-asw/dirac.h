#ifndef DIRAC_H
#define DIRAC_H
#include <vector>
#include <iomanip>
#include <GSLpp/special_functions.h>
#include <GSLpp/vector.h>
#include <GSLpp/matrix.h>
#include <GSLpp/ode.h>
#include <GSLpp/interp.h>
#include <algorithm>
#include <fstream>

#include <numerical-mesh.h>
#include <numerical-mesh-integration.h>

class Radial_Dirac_Equation{
    /***********************************************************************//**
    * Find the point where E - v(X_i) = 0 (classical inversion point).\n
    * \t__Iter_res__ - iterator to result container (e.g. std::vector)\n
    * \t__Iter_v__ - iterator to v container (e.g. std::vector)\n
    * Input:\n
    * \t__res_start__ - Iterator pointing to the first element of the result container.\n
    * \t__res_end__ - Iterator pointing to the element one past the last element of the result container.\n
    *\t__v_start__ - Iterator pointing to the first element of  the g container
    ***************************************************************************/
    template<class Iter_res, class Iter_v>
    std::array<Iter_res, 4> find_inversion_point(
        std::array<Iter_res, 4> res_start,
        std::array<Iter_res, 4> res_end,
        Iter_v v_start)
    {
        int init_sign = signum(energy_m - *v_start);
        std::array<Iter_res, 4> current = res_start;
        Iter_v current_v = v_start;
        while(current != res_end && signum(energy_m - (*current_v)) == init_sign){
            std::for_each (current.begin(), current.end(), [] (Iter_res it) { it++;});
            current_v++;
        }
        return current;
    }

protected:
    const Exponential_mesh<1, double>& mesh_m;
    std::vector<double> g_m, f_m;
    double energy_m;

    double F_norm() const
    {
        std::vector<double> integrand(mesh_p.size(), 0);
        for(size_t i = 0; i < mesh_m.size(); i++){
            integrand[i] = (GSL::pow_int(g_m[i], 2) +
                            GSL::pow_int(f_m[i], 2))*mesh_m.drx(i);
        }
        return simpson_integral(mesh_m, [&](const double& i)
            {
                return integrand[static_cast<size_t>(i)];
            });
    }

    template<class Iter, class T = double>
    double variational_de(Iter inv)
    {
        Iter P = g_m.begin();
        Iter Q = f_m.begin();
        for( ; P != inv ; P++, Q++){}

        return 12*(*P)*(*Q);
    }

public:
    Radial_Dirac_Equation(const double e_min, const double e_max,
      const std::vector<double>& v, const size_t l, const std::vector<double>& left_init_n,
      const std::vector<double>& right_init_n, const Logarithmic_mesh& mesh,
      const double tol = 1e-10)
      : Dirac_Equation(e_min, e_max, v, left_init_n, right_init_n, 1., tol),
        mesh_m(mesh)
    {
        if(mesh.size() != v.size()){
            throw std::runtime_error("Length of potential vector does not match length of radial mesh!");
        }
        auto r2 = mesh.r2_begin();
        auto v_c = v_m.begin();
        auto r = mesh.r_begin();
        for( ; v_c != v_m.end(); v_c++, r2++, r++){
            *v_c += static_cast<double>(l*(l + 1))/(*r2);
        }

        auto l_i = l_init_m.begin();
        for(auto  drx_i = mesh.drx_begin(); l_i != l_init_m.end(); l_i++, drx_i++){
            *l_i /= *drx_i;
        }
        auto r_i = r_init_m.rbegin();
        for(auto drx_i = mesh.drx_rbegin(); r_i != r_init_m.rend(); r_i++, drx_i++){
            *r_i /= *drx_i;
        }

    }

    void solve(const  size_t nodes) override
    {
        solve(nodes, 0.5*(e_min_m + e_max_m));
    }

    void solve(const  size_t nodes, const  double e_guess) override
    {    }

    double norm() const override
    {
        std::vector<double> integrand(mesh_p.size());
        for(size_t i = 0; i < mesh_p.size(); i++){
            integrand[i] = GSL::pow_int(psi_m[0][i], 2) +
                           GSL::pow_int(psi_m[2][i], 2);
        }
        return mesh_p.integrate_simpson(integrand);
    }

    void normalize() override
    {
        double n = this->norm();
        auto drx_i = mesh_m.drx_begin();
        auto g = g_m.begin(), f = f_m.begin();
        for(; g != g_m.end(), f != f_m.end(); g++, f++, drx_i++){
            *g /= std::sqrt(n);
            *f /= std::sqrt(n);
        }
    }

};

class Radial_Dirac_Equation_Central_Potential : public Radial_Dirac_Equation{
public:
    Radial_Dirac_Equation_Central_Potential(const std::vector<double>& v, const size_t l, const std::vector<double>& left_init,
      const std::vector<double>& right_init, const Logarithmic_mesh& mesh,
      const double tol = 1e-10)
    : Radial_Dirac_Equation(0, 0, v, l, left_init,
        right_init, mesh, tol)
    {
        e_min_m = *std::min_element(v_m.begin()+1, v_m.end());
        e_max_m = *(v_m.end());
    }
};
#endif // DIRAC_H
