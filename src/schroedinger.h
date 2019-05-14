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
		int init_sign = signum(*v_start);
		Iter_res current = res_start;
		Iter_v current_v = v_start;
		while(current != res_end && signum(energy_m - (*current_v)) == init_sign){
			current++;
			current_v++;
		}
		return current;
	}

    virtual void normalize()
    {
        double norm = 0;
        for(auto it = psi_m.begin(); it != psi_m.end(); it++){
            norm += GSL::pow_int(*it, 2)*h_m;
        }
        for(auto it = psi_m.begin(); it != psi_m.end(); it++){
            *it /= norm;
        }
    }

private:
    double h_m;

public:

    Schroedinger_Equation(double e_min, double e_max, std::vector<double>& v, std::vector<double>& left_init, std::vector<double>& right_init, double h = 1e-3, double tol = 1e-10)
     : energy_m(0), e_min_m(e_min), e_max_m(e_max), v_m(v), psi_m(v.size(), 0),
       l_init_m(left_init), r_init_m(right_init), tol_m(tol), h_m(h)
    {}

    Schroedinger_Equation(double e_min, double e_max, std::vector<double>& v,
    std::vector<double>& left_init, std::vector<double>& right_init, Mesh&
    mesh, double tol = 1e-10)
        : Schroedinger_Equation(e_min, e_max, v, left_init, right_init,
          mesh.dx(), tol)
    {};

    virtual ~Schroedinger_Equation() = default;

    virtual void solve()
    {
        std::vector<double> g(v_m.size(), 0), s(v_m.size(), 0);
        Numerov_solver sol;
        auto inv = psi_m.end();
        auto end_point = psi_m.rbegin();
        for( auto tmp = r_init_m.begin() ; tmp != r_init_m.end(); tmp++, end_point++ ){}
        int de;

        while(std::abs(e_max_m - e_min_m) > tol_m){
            energy_m = 0.5*(e_min_m + e_max_m);
            inv = find_inversion_point(end_point, psi_m.rend(), v_m.rbegin()).base();
            for(auto it = v_m.begin(), g_c = g.begin(); it != v_m.end(); it++, g_c++){
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

    void solve(size_t nodes)
    {
        std::vector<double> g(v_m.size(), 0), s(v_m.size(), 0);
        Numerov_solver sol;
        auto inv = psi_m.end();
        for( auto tmp = r_init_m.begin() ; tmp != r_init_m.end(); tmp++, inv-- ){}
        size_t n = nodes + 1;
        while(n != nodes){
            energy_m = 0.5*(e_min_m + e_max_m);

            for(auto it = v_m.begin(), g_c = g.begin(); it != v_m.end(); it++, g_c++){
                *g_c = (energy_m - *it)*GSL::pow_int(h_m, 2);
            }
            inv = sol.solve(psi_m.begin(), psi_m.end(), g.begin(), g.end(),
                  s.begin(), s.end(), l_init_m.begin(), l_init_m.end(),
                  r_init_m.begin(), r_init_m.end(), inv);
            n = static_cast<size_t>(sol.count_nodes(psi_m.begin(),
                psi_m.end()));
            if(n > nodes){
                e_max_m = energy_m;
            }else if(n < nodes){
                e_min_m = energy_m;
            }
            solve();
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
    Logarithmic_mesh mesh_p;

    double norm()
    {
        double res = 0;
        auto drx_i = mesh_p.drx_begin();
        auto it = psi_m.begin();
        for(; it != psi_m.end(); it++, drx_i++){
            res += GSL::pow_int(*it*(*drx_i), 2);
        }
        return res;
    }

    void normalize()
    {
        double n = norm();
        for(auto it = psi_m.begin(); it != psi_m.end(); it++){
            *it /= std::sqrt(n);
        }

    }
public:
    virtual ~Radial_Schroedinger_Equation() = default;
    Radial_Schroedinger_Equation(double e_min, double e_max,
      std::vector<double>& v, size_t l, std::vector<double>& left_init,
      std::vector<double>& right_init, Logarithmic_mesh& mesh,
      double tol = 1e-10)
      : Schroedinger_Equation(e_min, e_max, v, left_init, right_init, 1., tol),
        mesh_p(mesh)
    {
        if(mesh.size() != v.size()){
            throw std::runtime_error("Length of potential vector does not match length of radial mesh!");
        }
        auto r2 = mesh.r2_begin();
        auto v_c = v_m.begin();
        for( ; v_c != v_m.end(); v_c++, r2++){
            *v_c += static_cast<double>(l*(l + 1))/(*r2);
        }
    }

    using Schroedinger_Equation::solve;

    void solve(size_t nodes)
    {
        double e_max_old, e_min_old;
        std::vector<double> g(v_m.size(), 0), s(v_m.size(), 0);
        Numerov_solver sol;
        auto inv = psi_m.rend();
        auto end_point = psi_m.rend();
        for( auto tmp = r_init_m.begin() ; tmp != r_init_m.end(); tmp++, end_point++ ){}
        size_t n = nodes + 1;
        while(n != nodes){
            energy_m = 0.5*(e_max_m + e_min_m);
            inv = find_inversion_point(psi_m.rbegin(), end_point, v_m.rbegin());

            for(auto v_i = v_m.begin(), g_i = g.begin(), drx_i = mesh_p.drx_begin(); v_i != v_m.end(); v_i++, g_i++, drx_i++){
                *g_i = (energy_m - *v_i)*GSL::pow_int(*drx_i, 2) -
                    GSL::pow_int(mesh_p.A(), 2)/4;
            }

            inv = sol.solve(psi_m.rbegin(), psi_m.rend(), g.rbegin(), g.rend(),
                  s.rbegin(), s.rend(), r_init_m.rbegin(), r_init_m.rend(),
                  l_init_m.begin(), l_init_m.end(), inv);

            n = static_cast<size_t>(sol.count_nodes(++psi_m.begin(),
                --psi_m.end()));

            if(n > nodes){
                e_max_m = energy_m;
            }else if(n < nodes){
                e_min_m = energy_m;
            }
            e_max_old = e_max_m;
            e_min_old = e_min_m;
            std::cout << e_min_m << " < " << energy_m << " < " << e_max_m << "\n";
            solve();
            n = static_cast<size_t>(sol.count_nodes(++psi_m.begin(),
                --psi_m.end()));
            if(n > nodes){
                e_min_m = e_min_old;
                e_max_m = energy_m;
            }else if(n < nodes){
                e_max_m = e_max_old;
                e_min_m = energy_m;
            }
            std::cout << "number of nodes = " << nodes << "\n";
        }
    }

    void solve()
    {
        std::vector<double> g(v_m.size(), 0), s(v_m.size(), 0);
        Numerov_solver sol;
        auto inv = psi_m.rend();
        auto end_point = psi_m.rend();
        for( auto tmp = r_init_m.begin() ; tmp != r_init_m.end(); tmp++, end_point++ ){}
        double de = 10*tol_m;

        while(std::abs(de) > tol_m && std::abs(e_max_m - e_min_m) > tol_m){
            inv = find_inversion_point(psi_m.rbegin(), end_point,
                v_m.rbegin());

            auto g_i = g.begin();
            auto drx_i = mesh_p.drx_begin();
            for(auto v_i = v_m.begin(); v_i != v_m.end(); v_i++, g_i++, drx_i++){
                *g_i = (energy_m - *v_i)*GSL::pow_int(*drx_i, 2)
                    - GSL::pow_int(mesh_p.A(), 2)/4;
            }

            inv = sol.solve(psi_m.rbegin(), psi_m.rend(), g.rbegin(), g.rend(),
                  s.rbegin(), s.rend(), r_init_m.rbegin(), r_init_m.rend(),
                  l_init_m.begin(), l_init_m.end(), inv);

            if(inv != psi_m.rend()){
                de = variational_de(g.rbegin(), inv);
                std::cout << "Energy = " << energy_m << ", de = " << de << "\n";
            }else{
                de = 0;
            }
            if(de > 0){
                e_min_m = energy_m;
            }else if(de < 0){
                e_max_m = energy_m;
            }
            energy_m += de;
            energy_m = std::min(energy_m, e_max_m);
            energy_m = std::max(energy_m, e_min_m);
        }
        normalize();
    }

    template<class Iter, class T = double>
    double variational_de(Iter g_start, Iter inv)
    {
        Iter g_i = g_start;
        auto F_i = psi_m.rbegin();
        auto drx_i = mesh_p.drx_rbegin();
        for( ; F_i != inv; F_i++, g_i++, drx_i++){
        }
        T n = std::sqrt(norm());
        T Fi = *F_i/n;
        T Fm1 = *(F_i - 1)/n;
        T Fp1 = *(F_i + 1)/n;


        auto f = [](Iter g)->T{ return 1 + *g/12;};

        // T Fcusp = (Fm1*f(g_i - 1) + Fp1*f(g_i + 1))/(12 - 10*f(g_i));
        T Fcusp = (Fm1*f(g_i - 1) + Fp1*f(g_i + 1) + 10*Fi*f(g_i))/12;
        T df = f(g_i)*(Fi/Fcusp - 1.);
        return 12*Fcusp*Fcusp*df;
    }

};

class Radial_Schroedinger_Equation_Central_Potential : public Radial_Schroedinger_Equation{
public:
    ~Radial_Schroedinger_Equation_Central_Potential() = default;
    Radial_Schroedinger_Equation_Central_Potential(std::vector<double>& v, size_t l, std::vector<double>& left_init,
      std::vector<double>& right_init, Logarithmic_mesh& mesh,
      double tol = 1e-10)
    : Radial_Schroedinger_Equation(0, 0, v, l, left_init,
        right_init, mesh, tol)
    {
        if(l > 0){
            e_min_m = *std::min_element(++v_m.begin(), v_m.end());
        }else{
            e_min_m = v_m.back() - 1;
        }
        e_max_m = *--(--v_m.end());
        std::cout << "Energy inside range [" << e_min_m << ", " << e_max_m << "]\n";
    }
};
#endif // SCHOEDINGER_H
