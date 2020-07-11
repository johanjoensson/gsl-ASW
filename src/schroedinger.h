#ifndef SCHOEDINGER_H
#define SCHOEDINGER_H
#include <vector>
#include <iomanip>
#include "numerov_solver.h"
#include "log_mesh.h"
#include <GSLpp/special_functions.h>
#include <algorithm>
#include <fstream>

class Schroedinger_Equation{
protected:
    double energy_m, e_min_m, e_max_m;
    std::vector<double> v_m;
    std::vector<double> psi_m, l_init_m, r_init_m;
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
	Iter_v find_inversion_point(Iter_res res_start, Iter_res res_end,  Iter_v
		v_start)
	{
		int init_sign = signum(energy_m - *v_start);
		Iter_res current = res_start;
		Iter_v current_v = v_start;
		while(current != res_end && signum(energy_m - (*current_v)) == init_sign){
			current++;
			current_v++;
		}
		return current;
	}

private:
    const double h_m;

public:

    Schroedinger_Equation(const double e_min, const double e_max, const std::vector<double>& v, const std::vector<double>& left_init, const std::vector<double>& right_init, const double h = 1e-3, const double tol = 1e-10)
     : energy_m(0), e_min_m(e_min), e_max_m(e_max), v_m(v), psi_m(v.size(), 0),
       l_init_m(left_init), r_init_m(right_init), tol_m(tol), h_m(h)
    {}

    Schroedinger_Equation(const double e_min, const double e_max, const std::vector<double>& v,
    const std::vector<double>& left_init, const std::vector<double>& right_init, const Mesh&
    mesh, const double tol = 1e-10)
        : Schroedinger_Equation(e_min, e_max, v, left_init, right_init,
          mesh.dx(), tol)
    {};

    // Schroedinger_Equation() = default;
    Schroedinger_Equation(const Schroedinger_Equation &) = default;
    Schroedinger_Equation(Schroedinger_Equation &&) = default;
    virtual ~Schroedinger_Equation() = default;

    Schroedinger_Equation& operator=(const Schroedinger_Equation&) = default;
    Schroedinger_Equation& operator=(Schroedinger_Equation&&) = default;

    virtual double norm() const
    {
        double res = 0;
        for(auto it = psi_m.begin(); it != psi_m.end(); it++){
            res += GSL::pow_int(*it, 2)*h_m;
        }
        return res;
    }

    virtual void solve()
    {
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

    }

    virtual void solve(size_t nodes)
    {
        solve(nodes, 0.5*(e_min_m + e_max_m));
    }

    virtual void solve(const size_t nodes, const double e_guess)
    {
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
    }

    virtual void normalize()
    {
        double n = norm();
        for(auto it = psi_m.begin(); it != psi_m.end(); it++){
            *it /= n;
        }
    }

    std::vector<double>& psi()
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

class Radial_Schroedinger_Equation : public Schroedinger_Equation{
protected:
    const Logarithmic_mesh& mesh_p;

    double F_norm() const
    {
        std::vector<double> integrand(mesh_p.size());
        for(size_t i = 0; i < mesh_p.size(); i++){
            integrand[i] = GSL::pow_int(psi_m[i], 2)*mesh_p.drx(i);
        }
        return mesh_p.integrate_simpson(integrand);
    }

    template<class Iter, class T = double>
    double variational_de(Iter g_start, Iter inv)
    {
        Iter g_i = g_start;
        auto F_i = psi_m.rbegin();
        auto drx_i = mesh_p.drx_rbegin();
        for( ; F_i != inv; F_i++, g_i++, drx_i++){}

        T n = F_norm();
        T Fi = *F_i;
        T Fm1 = *(F_i - 1);
        T Fp1 = *(F_i + 1);
        T drx = *drx_i;


        auto f = [](T g)->T{ return 1 + g/12;};


        T Fcusp = (Fm1*f(*(g_i - 1)) + Fp1*f(*(g_i + 1)))/(12 - 10*f(*g_i));
        // T Fcusp = (Fm1*f(*(g_i - 1)) + Fp1*f(*(g_i + 1)) + 10*Fi*f(*g_i))/12;
        T df = f(*g_i)*(Fi/Fcusp - 1.);
        return 12*Fcusp*Fcusp*df*drx/n;
    }

public:
    Radial_Schroedinger_Equation(const double e_min, const double e_max,
      const std::vector<double>& v, const size_t l, const std::vector<double>& left_init_n,
      const std::vector<double>& right_init_n, const Logarithmic_mesh& mesh,
      const double tol = 1e-10)
      : Schroedinger_Equation(e_min, e_max, v, left_init_n, right_init_n, 1., tol),
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

    using Schroedinger_Equation::solve;

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
        auto inv = psi_m.rend();
        auto end_point = psi_m.rbegin();
        for( auto tmp = r_init_m.begin() ; tmp != r_init_m.end(); tmp++, end_point++ ){}
        while(n != nodes || (std::abs(de) > tol_m && std::abs(e_max_m - e_min_m) > tol_m)){
            inv = find_inversion_point(end_point, psi_m.rend(), v_m.rbegin());

            auto v_i = v_m.begin();
            auto g_i = g.begin();
            auto drx_i = mesh_p.drx_begin();
            for(; v_i != v_m.end(); v_i++, g_i++, drx_i++){
                *g_i = (energy_m - *v_i)*GSL::pow_int(*drx_i, 2) -
                    GSL::pow_int(mesh_p.A(), 2)/4;
            }

            sol.solve(psi_m.rbegin(), psi_m.rend(), g.rbegin(), g.rend(),
                  s.rbegin(), s.rend(), r_init_m.rbegin(), r_init_m.rend(),
                  l_init_m.begin(), l_init_m.end(), inv);

            n = static_cast<size_t>(sol.count_nodes(++psi_m.rbegin(),--psi_m.rend()));

            if(n > nodes){
                e_max_m = energy_m;
            }else if(n < nodes){
                e_min_m = energy_m;
            }

            if(n == nodes){
                if(inv != psi_m.rend()){
                    de = variational_de(g.rbegin(), inv);
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
            integrand[i] = GSL::pow_int(psi_m[i], 2);
        }
        return mesh_p.integrate_simpson(integrand);
    }

    void normalize() override
    {
        double n = this->norm();
        auto drx_i = mesh_p.drx_begin();
        for(auto it = psi_m.begin(); it != psi_m.end(); it++, drx_i++){
            *it /= std::sqrt(n);
        }
    }

};

class Radial_Schroedinger_Equation_Central_Potential : public Radial_Schroedinger_Equation{
public:
    Radial_Schroedinger_Equation_Central_Potential(const std::vector<double>& v, const size_t l, const std::vector<double>& left_init,
      const std::vector<double>& right_init, const Logarithmic_mesh& mesh,
      const double tol = 1e-10)
    : Radial_Schroedinger_Equation(0, 0, v, l, left_init,
        right_init, mesh, tol)
    {
        e_min_m = *std::min_element(v_m.begin()+1, v_m.end());
        e_max_m = *(v_m.end());
    }
};
#endif // SCHOEDINGER_H
