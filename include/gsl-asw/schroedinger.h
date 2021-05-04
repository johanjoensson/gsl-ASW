#ifndef SCHOEDINGER_H
#define SCHOEDINGER_H
#include <vector>
#include <iomanip>
#include "numerov_solver.h"
#include "numerical-mesh.h"
#include "numerical-mesh-integration.h"
#include <GSLpp/special_functions.h>
#include <algorithm>
#include <fstream>

class Schroedinger_Equation{
    template<class T>
	static int signum(T val)
     noexcept
	{
		return (val > T(0)) - (val < T(0));
	}

protected:
    const Mesh_base<1, double>& mesh_m;
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
	Iter_res find_inversion_point(Iter_res res_start, Iter_res res_end,  Iter_v
		v_start)
     const noexcept
	{
		int init_sign = signum(energy_m - *v_start);
		Iter_res current = res_start;
		Iter_v current_v = v_start;
        while(current != res_end && signum(energy_m - (*current_v)) == init_sign){
			current++;
			current_v++;
		}
        return ++current;
	}

private:

public:

    Schroedinger_Equation(const double e_min, const double e_max, const std::vector<double>& v,
    const std::vector<double>& left_init, const std::vector<double>& right_init, const Mesh_base<1, double>&
    mesh, const double tol = 1e-10)
    : mesh_m(mesh), energy_m(0), e_min_m(e_min), e_max_m(e_max), v_m(v), psi_m(v.size(), 0),
      l_init_m(left_init), r_init_m(right_init), tol_m(tol)
    {};

    Schroedinger_Equation(const Schroedinger_Equation &) = default;
    Schroedinger_Equation(Schroedinger_Equation &&) = default;
    virtual ~Schroedinger_Equation() = default;

    Schroedinger_Equation& operator=(const Schroedinger_Equation&) = default;
    Schroedinger_Equation& operator=(Schroedinger_Equation&&) = default;


    virtual void solve()
    {
        std::vector<double> g(v_m.size(), 0), s(v_m.size(), 0);
        Numerov_solver sol;
        auto end_point = psi_m.rbegin();
        for( auto tmp = r_init_m.begin() ; tmp != r_init_m.end(); tmp++,
            end_point++ ){}
        auto inv = find_inversion_point(end_point, psi_m.rend(),
                    v_m.rbegin()).base();
        int de = 0;

        const double h = mesh_m.dr(0);

        while(std::abs(e_max_m - e_min_m) > tol_m){

            energy_m = 0.5*(e_min_m + e_max_m);
            for(auto it = v_m.begin(), g_c = g.begin(); it != v_m.end(); it++,
                g_c++){
                *g_c = (energy_m - *it)*GSL::pow_int(h, 2);
            }
            inv = find_inversion_point(end_point, psi_m.rend(), v_m.rbegin()).base();
            if(inv != sol.solve(psi_m.begin(), psi_m.end(), g.begin(), g.end(),
                  s.begin(), s.end(), l_init_m.begin(), l_init_m.end(),
                  r_init_m.begin(), r_init_m.end(), inv)){
                      throw std::runtime_error("Unexpected change of matching point!");
                  }

            if(inv != psi_m.end()){
                de = sol.derivative_diff(psi_m.begin(), psi_m.end(), inv,
                g.begin(), s.begin(), GSL::pow_int(h, 3));
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
        const double h = mesh_m.dr(0);

        while(n != nodes || std::abs(e_max_m - e_min_m) > tol_m){

            for(auto it = v_m.begin(), g_c = g.begin(); it != v_m.end(); it++,
                g_c++){

                *g_c = (energy_m - *it)*GSL::pow_int(h, 2);
            }
            inv = find_inversion_point(end_point, psi_m.rend(), v_m.rbegin()).base();
            if(inv != sol.solve(psi_m.begin(), psi_m.end(), g.begin(), g.end(),
                  s.begin(), s.end(), l_init_m.begin(), l_init_m.end(),
                  r_init_m.begin(), r_init_m.end(), inv)){
                      throw std::runtime_error("Unexpected change of matching point!");
                  }
            n = static_cast<size_t>(sol.count_nodes(psi_m.begin(), psi_m.end()));


            if(n > nodes){
                e_max_m = energy_m;
            }else if(n < nodes){
                e_min_m = energy_m;
            }

            if(inv != psi_m.end()){
                de = sol.derivative_diff(psi_m.begin(), psi_m.end(), inv,
                g.begin(), s.begin(), GSL::pow_int(h, 3));

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

    virtual double norm()
     const noexcept
    {
        double res = simpson_integral<double>(mesh_m,
            [this](const double& i)
            {
                return GSL::pow_int(psi_m[static_cast<size_t>(i)], 2);
            }
        );
        return res;
    }

    void normalize(const double N = 1)
     noexcept
    {
        double n = norm();
        std::transform(psi_m.begin(), psi_m.end(), psi_m.begin(),
        [=](const double val){
            return val*std::sqrt(N/n);
        });
    }

    std::vector<double>& psi()
     noexcept
    {
        return psi_m;
    }

    std::vector<double> psi()
     const noexcept
    {
        return psi_m;
    }

    double& psi(const size_t i)
     noexcept
    {
        return psi_m[i];
    }

    double psi(const size_t i)
     const noexcept
    {
        return psi_m[i];
    }

    std::vector<double>& v()
     noexcept
    {
        return v_m;
    }

    std::vector<double> v()
     const noexcept
    {
        return v_m;
    }

    double& v(const size_t i)
     noexcept
    {
        return v_m[i];
    }

    double v(const size_t i)
     const noexcept
    {
        return v_m[i];
    }

    double e()
     const noexcept
    {
        return energy_m;
    }

};

class Radial_Schroedinger_Equation : public Schroedinger_Equation{
protected:

    template<class Iter, class T = double>
    double variational_de(Iter g_start, Iter inv)
     const noexcept
    {
        Iter g_i = g_start;
        auto mesh_i = mesh_m.begin();
        for(auto F_i = psi_m.begin() ; F_i != inv; F_i++, g_i++, mesh_i++){}
        T n = simpson_integral<double>(mesh_m,
            [this](const double& i)
            {
                return GSL::pow_int(psi_m[static_cast<size_t>(i)], 2)*mesh_m.dr(i);
            }
        );
        T Fi = *(inv - 1);
        T Fm1 = *(inv - 2);
        T Fp1 = *inv;

        auto f = [](T g)->T{ return 1 + g/12.;};

        T Fcusp = (Fm1*f(*(g_i - 2)) + Fp1*f(*g_i) + 10*Fi*f(*(g_i - 1)))/12.;

        T df = f(*(g_i - 1))*(Fi/Fcusp - 1);

        return 12.*Fcusp*Fcusp*df/n;
    }

public:
    Radial_Schroedinger_Equation(const double e_min, const double e_max,
      const std::vector<double>& v, const size_t l, const std::vector<double>& left_init_n,
      const std::vector<double>& right_init_n, const Mesh_base<1, double>& mesh,
      const double tol = 1e-10)
      : Schroedinger_Equation(e_min, e_max, v, left_init_n, right_init_n, mesh, tol)
    {

        if(mesh_m.size() != v.size()){
            throw std::runtime_error("Length of potential vector does not match length of radial mesh!");
        }

        auto v_c = v_m.begin();
        for(auto mesh_i = mesh_m.begin() ; v_c != v_m.end(); v_c++, mesh_i++){
            *v_c += static_cast<double>(l*(l + 1))/mesh_i->r2();
        }
    }

    using Schroedinger_Equation::solve;

    void solve(const  size_t nodes)
     noexcept override
    {
        solve(nodes, 0.5*(e_min_m + e_max_m));
    }

    void solve(const  size_t nodes, const  double e_guess)
     noexcept override
    {
        std::vector<double> g(v_m.size(), 0), s(v_m.size(), 0);
        Numerov_solver sol;
        size_t n = nodes + 1;
        double de = 10*tol_m;
        energy_m = e_guess;

        auto inv = psi_m.begin();

        auto end_point = psi_m.end(), v_end = v_m.end();
        for( auto tmp = r_init_m.begin(); tmp != r_init_m.end(); tmp++, end_point--, v_end-- ){}


        std::vector<double> left_init;
        std::transform(l_init_m.begin(), l_init_m.end(), mesh_m.begin(), std::back_inserter(left_init),
            [=](const double psi, const Mesh_base<1, double>::mesh_point mi)
            {
                return psi*mi.r()/std::sqrt(mi.dr());
            }
        );
        std::vector<double> right_init;
        std::transform(r_init_m.begin(), r_init_m.end(), mesh_m.end() - static_cast<double>(r_init_m.size() - 1), std::back_inserter(right_init),
            [=](const double psi, const Mesh_base<1, double>::mesh_point mi)
            {
                return psi*mi.r()/std::sqrt(mi.dr());
            }
        );


        // Locate potential minima
        auto v_start = v_m.begin();
        for(auto tmp = left_init.begin(); tmp != left_init.end(); tmp++, v_start++){}
        v_start = std::min_element(v_start, v_end);

        // Inversion point must be to the right of the inversion minima, so we
        // can start the search for the inversion point, at the potential minima
        auto start_point = psi_m.begin();
        for(auto tmp = v_m.begin() ; tmp != v_start; tmp++, start_point++){}



        size_t variational_steps = 0;
        while(n != nodes || (std::abs(de) > tol_m && std::abs((e_max_m*e_max_m - e_min_m*e_min_m)/(e_max_m + e_min_m)) > tol_m)){
            inv = find_inversion_point(start_point, end_point, v_start);
            auto v_i = v_m.begin();
            auto g_i = g.begin();
            auto mesh_i = mesh_m.begin();
            for(; v_i != v_m.end(); v_i++, g_i++, mesh_i++){
                *g_i =  (energy_m - *v_i)*GSL::pow_int(mesh_i->dr(), 2)+
                        (
                            -3./4*GSL::pow_int(mesh_i->d2r()/mesh_i->dr(), 2)
                            +1./2*mesh_i->d3r()/mesh_i->dr()
                        );
            }

            sol.solve(psi_m.begin(), psi_m.end(), g.begin(), g.end(),
                s.begin(), s.end(), left_init.begin(), left_init.end(),
                right_init.rbegin(), right_init.rend(), inv);

            n = static_cast<size_t>(sol.count_nodes(++psi_m.begin(),--psi_m.end()));

            if(n > nodes){
                e_max_m = energy_m;
                energy_m = std::max(0.5*(e_min_m + e_max_m), e_max_m - 1);
            }else if(n < nodes){
                e_min_m = energy_m;
                energy_m = std::min(0.5*(e_min_m + e_max_m), e_min_m + 1);
            }else if(n == nodes){
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
                energy_m = std::min(std::max(energy_m, e_min_m), e_max_m);
                variational_steps++;
            }
        }

        auto it = psi_m.begin();
        auto mesh_i = mesh_m.begin();
        for(; it != psi_m.end(); it++, mesh_i++){
            *it *= std::sqrt(mesh_i->dr())/mesh_i->r();
        }

        it = psi_m.begin();
        for(auto li = l_init_m.begin(); li != l_init_m.end(); li++, it++){
            *it = *li;
        }

        auto scale = 1.;
        if(std::abs(r_init_m.back() - psi_m.back()) > tol_m){
            scale = r_init_m.back()/psi_m.back();
        }

        it = psi_m.begin();
        mesh_i = mesh_m.begin();
        for(; it != psi_m.end(); it++, mesh_i++){
            *it *= scale;
        }

    }

    virtual double norm()
     const noexcept override
    {
        double res = simpson_integral<double>(mesh_m,
            [this](const double& i)
            {
                return GSL::pow_int(mesh_m.r(i)*psi_m[static_cast<size_t>(i)], 2);
            }
        );
        return res;
    }

};

class Radial_Schroedinger_Equation_Central_Potential : public Radial_Schroedinger_Equation{
public:
    Radial_Schroedinger_Equation_Central_Potential(const std::vector<double>& v, const size_t l, const std::vector<double>& left_init,
      const std::vector<double>& right_init, const Mesh_base<1, double>& mesh,
      const double tol = 1e-10)
    : Radial_Schroedinger_Equation(0, 0, v, l, left_init,
        right_init, mesh, tol)
    {
        e_min_m = *std::min_element(v_m.begin()+1, v_m.end());
        e_max_m = v_m.back();
    }
};
#endif // SCHOEDINGER_H
