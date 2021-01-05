#ifndef DIRAC_H
#define DIRAC_H
#include <vector>
#include <iomanip>
#include "numerov_solver.h"
#include "log_mesh.h"
#include <GSLpp/special_functions.h>
#include <GSLpp/vector.h>
#include <GSLpp/matrix.h>
#include <GSLpp/ode.h>
#include <GSLpp/interp.h>
#include <algorithm>
#include <fstream>

class Dirac_Equation{
protected:
    double energy_m, e_min_m, e_max_m;
    std::vector<double> v_m;
    std::array<std::vector<double>, 4> psi_m
    std::array<std::vector<double>, 4>l_init_m, r_init_m;
    double tol_m;

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
	std::array<Iter_res, 4> find_inversion_point(std::array<Iter_res, 4> res_start,
        std::array<Iter_res, 4> res_end,  Iter_v v_start)
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

private:
    const double h_m;

public:

    Dirac_Equation(const double e_min, const double e_max, const std::vector<double>& v, const std::array<std::vector<double>, 4>& left_init, const std::array<std::vector<double>, 4>& right_init, const double h = 1e-3, const double tol = 1e-10)
     : energy_m(0), e_min_m(e_min), e_max_m(e_max), v_m(v), psi_m(v.size(), 0),
       l_init_m(left_init), r_init_m(right_init), tol_m(tol), h_m(h)
    {}

    Dirac_Equation(const double e_min, const double e_max, const std::vector<double>& v,
    const std::array<std::vector<double>, 4>& left_init, const std::array<std::vector<double>, 4>& right_init, const Mesh&
    mesh, const double tol = 1e-10)
        : Dirac_Equation(e_min, e_max, v, left_init, right_init,
          mesh.dx(), tol)
    {};

    // Schroedinger_Equation() = default;
    Dirac_Equation(const Dirac_Equation &) = default;
    Dirac_Equation(Dirac_Equation &&) = default;
    virtual ~Dirac_Equation() = default;

    Dirac_Equation& operator=(const Dirac_Equation&) = default;
    Dirac_Equation& operator=(Dirac_Equation&&) = default;

    virtual double norm() const
    {
        double res = 0;
        auto psi1_it = psi_m[0].begin(), psi2_it = psi_m[1].begin(),
             psi3_it = psi_m[2].begin(), psi4_it  psi_m[3].begin();
        for(; psi1_it != psi_m[0].end(); psi1_it++, psi2_it++, psi3_it++, psi4_it++){
            res += (GSL::pow_int(*psi1_it, 2) + GSL::pow_int(*psi2_it, 2) +
                    GSL::pow_int(*psi3_it, 2) + GSL::pow_int(*psi4_it, 2))*h_m;
        }
        return res;
    }

    virtual void solve()
    {
/*
        std::vector<double> g(v_m.size(), 0), s(v_m.size(), 0);
        Numerov_solver sol;
        auto inv = psi_m.end();
        auto end_point = psi_m.rbegin();
        for( auto tmp = r_init_m.begin() ; tmp != r_init_m.end(); tmp++,
            end_point++ ){}
        int de;

        while(std::abs(e_max_m - e_min_m) > tol_m){

            energy_m = 0.5*(e_min_m + e_max_m);
            if(inv != find_inversion_point(end_point, psi_m.rend(),
                v_m.rbegin()).base()){
                    throw std::runtime_error("Unexpected change of matching point!");
                }

            for(auto it = v_m.begin(), g_c = g.begin(); it != v_m.end(); it++,
                g_c++){
                *g_c = (energy_m - *it)*GSL::pow_int(h_m, 2);
            }
            inv = sol.solve(psi_m.begin(), psi_m.end(), g.begin(), g.end(),
                  s.begin(), s.end(), l_init_m.begin(), l_init_m.end(),
                  r_init_m.begin(), r_init_m.end(), inv);

            if(inv != psi_m.end()){
                de = sol.derivative_diff(psi_m.begin(), psi_m.end(), inv,
                g.begin(), s.begin(), GSL::pow_int(h_m, 3));

            }else{
                de = 0;
            }

            if(de > 0 ) {
                e_max_m = energy_m;
            }else if(de < 0){
                e_min_m = energy_m;
            }else{
                e_max_m = e_min_m;
            }
        }
*/
    }

    virtual void solve(size_t nodes)
    {
        solve(nodes, 0.5*(e_min_m + e_max_m));
    }

    virtual void solve(const size_t nodes, const double e_guess)
    {
/*
        std::vector<double> g(v_m.size(), 0), s(v_m.size(), 0);
        Numerov_solver sol;
        auto inv = psi_m.end();
        auto end_point = psi_m.rbegin();
        for( auto tmp = r_init_m.begin() ; tmp != r_init_m.end(); tmp++, end_point++ ){}
        size_t n = nodes + 1;
        energy_m = e_guess;
        int de;

        while(n != nodes || std::abs(e_max_m - e_min_m) > tol_m){

            for(auto it = v_m.begin(), g_c = g.begin(); it != v_m.end(); it++,
                g_c++){

                *g_c = (energy_m - *it)*GSL::pow_int(h_m, 2);
            }
            inv = find_inversion_point(end_point, psi_m.rend(), v_m.rbegin()).base();
            if(inv != sol.solve(psi_m.begin(), psi_m.end(), g.begin(), g.end(),
                  s.begin(), s.end(), l_init_m.begin(), l_init_m.end(),
                  r_init_m.begin(), r_init_m.end(), inv)){
                      throw std::runtime_error("Unexpected change of matching point!");
                  }
            n = static_cast<size_t>(sol.count_nodes(psi_m.begin(),
                psi_m.end()));


            if(n > nodes){
                e_max_m = energy_m;
            }else if(n < nodes){
                e_min_m = energy_m;
            }

            if(inv != psi_m.end()){
                de = sol.derivative_diff(psi_m.begin(), psi_m.end(), inv,
                g.begin(), s.begin(), GSL::pow_int(h_m, 3));

            }else{
                de = 0;
            }

            if(n == nodes){
                if(de > 0 ) {
                    e_max_m = energy_m;
                }else if(de < 0){
                    e_min_m = energy_m;
                }else{
                    e_max_m = e_min_m;
                }
            }
            energy_m = 0.5*(e_min_m + e_max_m);
        }
*/
    }

    virtual void normalize()
    {
        double n = norm();
        auto normalize_component = [=](const std::vector<double>& component){
            for(auto it = component.begin(); it != component.end(); it++){
                *it /= n;
            }
        }
        std::for_each(psi_m.begin(), psi_m.end(), normalize_component);
    }

    std::array<std::vector<double>, 4>& psi()
    {
        return psi_m;
    }

    std::vector<double>& v()
    {
        return v_m;
    }

    double e()
    {
        return energy_m;
    }

};

class Radial_Dirac_Equation : public Dirac_Equation{
protected:
    const Logarithmic_mesh& mesh_p;

    double F_norm() const
    {
        std::vector<double> integrand(mesh_p.size(), 0);
        for(size_t i = 0; i < mesh_p.size(); i++){
            integrand[i] = (GSL::pow_int(psi_m[0][i], 2) +
                            GSL::pow_int(psi_m[2][i], 2))*mesh_p.drx(i);
        }
        return mesh_p.integrate_simpson(integrand);
    }

    template<class Iter, class T = double>
    double variational_de(Iter inv)
    {
        Iter P = psi_m[0].begin();
        Iter Q = psi_m[2].begin();
        for( ; != inv ; P++, Q++){}

        return 12*(*P)*(*Q);
    }

public:
    Radial_Dirac_Equation(const double e_min, const double e_max,
      const std::vector<double>& v, const size_t l, const std::vector<double>& left_init_n,
      const std::vector<double>& right_init_n, const Logarithmic_mesh& mesh,
      const double tol = 1e-10)
      : Dirac_Equation(e_min, e_max, v, left_init_n, right_init_n, 1., tol),
        mesh_p(mesh)
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

    using Dirac_Equation::solve;

    void solve(const  size_t nodes) override
    {
        solve(nodes, 0.5*(e_min_m + e_max_m));
    }

    void solve(const  size_t nodes, const  double e_guess) override
    {
        std::vector<double> g(v_m.size(), 0), s(v_m.size(), 0);
        Numerov_solver sol;
        size_t n = nodes + 1;
        double de = 10*tol_m;
        energy_m = e_guess;

        auto inv = psi_m.begin();

        auto end_point = psi_m.end(), v_end = v_m.end();
        for( auto tmp = r_init_m.begin(); tmp != r_init_m.end(); tmp++, end_point--, v_end-- ){}

        auto v_start = v_m.begin();
        for(auto tmp = l_init_m.begin(); tmp != l_init_m.end(); tmp++, v_start++){}
        v_start = std::min_element(v_start, v_end);

        auto start_point = psi_m.begin();
        for(auto tmp = v_m.begin() ; tmp != v_start; tmp++, start_point++){}

        while(n != nodes || (std::abs(de) > tol_m && std::abs(e_max_m - e_min_m) > tol_m)){
            inv = find_inversion_point(start_point, end_point, v_start);
            auto v_i = v_m.begin();
            auto g_i = g.begin();
            auto drx_i = mesh_p.drx_begin();
            for(; v_i != v_m.end(); v_i++, g_i++, drx_i++){
                *g_i = (energy_m - *v_i)*GSL::pow_int(*drx_i, 2) -
                    GSL::pow_int(mesh_p.A(), 2)/4;
            }

            sol.solve(psi_m.begin(), psi_m.end(), g.begin(), g.end(),
                  s.begin(), s.end(), l_init_m.begin(), l_init_m.end(),
                  r_init_m.rbegin(), r_init_m.rend(), inv);

            n = static_cast<size_t>(sol.count_nodes(++psi_m[0].begin(),--psi_m[0].end()));

            if(n > nodes){
                e_max_m = energy_m;
            }else if(n < nodes){
                e_min_m = energy_m;
            }

            if(n == nodes){
                if(inv != psi_m.end()){
                    de = variational_de(g.begin(), inv);
                }else{
                    de = 0;
                }

                if(de > 0){
                    e_min_m = energy_m;
                }else if(de < 0){
                    e_max_m = energy_m;
                }
                energy_m += de;
                energy_m = std::max(e_min_m, energy_m);
                energy_m = std::min(e_max_m, energy_m);

            }else{
                energy_m = 0.5*(e_min_m + e_max_m);
            }
        }
        auto it = psi_m.begin();
        auto drx_i = mesh_p.drx_begin();
        for(; it != psi_m.end(); it++, drx_i++){
            *it *= std::sqrt(*drx_i);
        }
    }

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
        auto drx_i = mesh_p.drx_begin();
        auto g = psi_m[0].begin(), f = psi_m[2].begin();
        for(; g != psi_m[0].end(), f != psi_m[2].end(); g++, f++, drx_i++){
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
