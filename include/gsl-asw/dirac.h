#ifndef DIRAC_H
#define DIRAC_H
#include <vector>
#include <iomanip>
#include <GSLpp/vector.h>
#include <GSLpp/matrix.h>
#include <GSLpp/ode.h>
#include <GSLpp/interp.h>
#include <GSLpp/special_functions.h>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <unistd.h>
#include <numerical-mesh.h>
#include <numerical-mesh-integration.h>



class Radial_Dirac_Equation{

	template<class T>
	static int signum(T val)
	 noexcept
	{
		return (val > T(0)) - (val < T(0));
	}

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
    Iter_res find_inversion_point(
        Iter_res res,
        Iter_res res_end,
        Iter_v v)
		 const noexcept
    {
        int init_sign = signum(energy_m - *v);
        while(res != res_end && signum(energy_m - (*v)) == init_sign){
            v++;
            res++;
        }
        return res;
    }

    size_t count_nodes(std::vector<double> fun, size_t i_inv)
	 const noexcept
    {
    	size_t n_nodes = 0;
		auto fi = fun.begin();
		for(;signum(*fi) == 0; fi++){
		}
		auto previous_sign = signum(*fi);
    	for(fi++; fi != fun.begin() + i_inv - 1; fi++){
            if (signum(*fi) != previous_sign && signum(*fi) != 0){
				n_nodes++;
				previous_sign = signum(*fi);
            }
    	}
    	return n_nodes;
    }



protected:
    const Exponential_mesh<1, double>& mesh_m;
    const std::vector<double>& v_m;
    const int kappa_m;
    std::vector<double> g_m, f_m;
	std::vector<double> left_init_g_m, left_init_f_m, right_init_g_m, right_init_f_m;
    double energy_m, abs_err_m;

    /*!*************************************************************************
    * Give a variational estimate of the energy correction required to get a
    * smooth solution, matching value and slope at the inversion point.
    ***************************************************************************/
    template<class T = double>
    T variational_de(T G_inv, T Fout_inv, T Fin_inv, T dr)
	 const noexcept
    {
		const double c = 2*137.035999174;
        T n = simpson_integral<double>(mesh_m,
            [this](const double& i)
            {
                return (GSL::pow_int(g_m[static_cast<size_t>(i)], 2)
					  + GSL::pow_int(f_m[static_cast<size_t>(i)], 2))
					  *mesh_m.dr(i);
            }
        );


        return c*G_inv*(Fout_inv - Fin_inv)/n*dr;
    }

public:
    Radial_Dirac_Equation(const Exponential_mesh<1, double>& mesh, const std::vector<double>& v,
        const int kappa, std::array<std::vector<double>, 2> left_init, std::array<std::vector<double>, 2> right_init, const double tol = 1e-10)
      : mesh_m(mesh), v_m(v), kappa_m(kappa), g_m(mesh.size()), f_m(mesh.size()),
	    left_init_g_m(left_init[0]), left_init_f_m(left_init[1]),
		right_init_g_m(right_init[0]), right_init_f_m(right_init[1]), energy_m(0), abs_err_m(tol)
    {}

	std::vector<double> g() const noexcept {return g_m;}
	std::vector<double> f() const noexcept {return f_m;}
	double g(const size_t i) const noexcept {return g_m[i];}
	double f(const size_t i) const noexcept {return f_m[i];}

	double energy() const noexcept {return energy_m;}

    /*!*************************************************************************
    * Solve the radial Dirac equations, starting with an energy guess equal to
    * e_guess, then refine the energy guess by making sure that the large
    * component of the solution has the required number of nodes. Then
    * use more refined methods, such as matching derivatives or variational
    * estimates of the energy, for solving the problem.
    ***************************************************************************/
    void solve(const  size_t nodes, const  double e_guess, const double e_min_guess, const double e_max_guess)
    {
        const double c = 2*137.035999174;

        double e_min = *std::min_element(std::begin(v_m) + 2, std::end(v_m));
		e_min = std::max(e_min_guess, e_min);
        double e_max = std::min(e_max_guess, v_m.back());
		std::vector<double> xs(v_m.size(), 0);
		std::iota(xs.begin(), xs.end(), 0);
		GSL::Interpolation Vlinear = GSL::Interpolation(xs, v_m, GSL::Interpolation_type::Linear);

        energy_m = e_guess;
        auto func = [=](double x, const std::vector<double>& y, std::vector<double>& f)
        {
            auto r = mesh_m.r(x);
			auto dr = mesh_m.dr(x);
			auto d2r = mesh_m.d2r(x);

			// Radial large component
			f[0] = dr*((-kappa_m/r - 0.5*d2r/(dr*dr))*y[0] + 1/c*(energy_m + c*c - Vlinear(x))*y[1]);
			// Radial small component
			f[1] = dr*(-1/c*(energy_m - Vlinear(x))*y[0] + (kappa_m/r - 0.5*d2r/(dr*dr))*y[1]);
        };

        GSL::ODE_solver solver(GSL::Runge_Kutta_Prince_Diamond, 2, abs_err_m, 100);

        solver.set_rhs(func);

		auto mesh_point = mesh_m.begin();
		for(auto g_l = left_init_g_m.begin(), g = g_m.begin();
			g_l != left_init_g_m.end(); g_l++, g++, mesh_point++){
			*g = *g_l*mesh_point->r()/std::sqrt(mesh_point->dr());
		}
		mesh_point = mesh_m.begin();
		for(auto f_l = left_init_f_m.begin(), f = f_m.begin();
			f_l != left_init_f_m.end(); f_l++, f++, mesh_point++){
			*f = *f_l*mesh_point->r()/std::sqrt(mesh_point->dr());
		}

        double tol = abs_err_m, de = 2*tol;
        size_t n = nodes - 1;
        while (n != nodes || (std::abs(de) > tol && std::abs(e_max - e_min) > tol)){

            size_t turning_index = g_m.size() - std::distance(std::rbegin(g_m),
                find_inversion_point(std::rbegin(g_m), std::rend(g_m), std::rbegin(v_m)));
			std::vector<double> y{	g_m[left_init_g_m.size() - 1],
									f_m[left_init_f_m.size() - 1] };
			for (size_t i = left_init_g_m.size(); i <= turning_index; i++){
                double x = i - 1;
                double h = 1;
				solver.apply_fixed_step(x, h, y);
                g_m[i] = y[0];
                f_m[i] = y[1];
            }
			std::vector<double> y_out{y[0], y[1]};

			n = count_nodes(g_m, turning_index);

			if (n == nodes) {
				for (size_t i = turning_index + 1; i < mesh_m.size(); i++){
					double x = i - 1;
					double h = 1;
					solver.apply_fixed_step(x, h, y);
					g_m[i] = y[0];
					f_m[i] = y[1];
				}
				n = count_nodes(g_m, turning_index);
			}
			n = nodes;

			if (n > nodes){
                e_max = energy_m;
				energy_m = 0.5*(e_max + e_min);
            }else if (n < nodes){
                e_min = energy_m;
				energy_m = 0.5*(e_max + e_min);
            } else if (n == nodes){
				solver.reset();
				mesh_point = mesh_m.end() - right_init_g_m.size();
				for(auto g_r = right_init_g_m.begin(), g = g_m.end() - right_init_g_m.size();
					g_r != right_init_g_m.end(); g_r++, g++, mesh_point++){
					*g = *g_r*mesh_point->r()/std::sqrt(mesh_point->dr());
				}
				mesh_point = mesh_m.end() - right_init_f_m.size();
				for(auto f_r = right_init_f_m.begin(), f = f_m.end() - right_init_f_m.size();
					f_r != right_init_f_m.end(); f_r++, f++, mesh_point++){
					*f = *f_r*mesh_point->r()/std::sqrt(mesh_point->dr());
				}
				y[0] = g_m[g_m.size() - right_init_g_m.size()];
				y[1] = f_m[f_m.size() - right_init_f_m.size()];
                for (size_t i = mesh_m.size() - right_init_g_m.size() - 1; i >= turning_index; i--){
                    double x = i + 1;
                    double h = -1;
					solver.apply_fixed_step(x, h, y);
                    g_m[i] = y[0];
                    f_m[i] = y[1];
                }
				std::vector<double> y_in{y[0], y[1]};
				solver.reset();

				auto scale = y_out[0]/y_in[0];
				for (size_t i = turning_index; i < mesh_m.size(); i++){
                    g_m[i] *= scale;

                }
				for (size_t i = turning_index; i < mesh_m.size(); i++){
					f_m[i] *= scale;
				}

			    de = variational_de(g_m[turning_index], y_out[1], y_in[1]*scale, mesh_m.dr(turning_index));
				if (de > 0){
                    e_min = energy_m;
                }else if (de < 0) {
                    e_max = energy_m;
                }
                energy_m += de;
                if(energy_m > e_max || energy_m < e_min){
                    energy_m = 0.5*(e_min + e_max);
                }
				std::cout << "de = " << de << "\n";
            }
        }
		auto g = g_m.begin(), f = f_m.begin();
		auto mesh_i = mesh_m.begin();
		double scale_g = 1, scale_f = 1;
		if(std::abs(right_init_g_m.back()) > 1e-16){
			scale_g = right_init_g_m.back()*mesh_m.r().back()/(g_m.back()*std::sqrt(mesh_m.dr().back()));
		}
		if(std::abs(right_init_f_m.back()) > 1e-16){
			scale_f = right_init_f_m.back()*mesh_m.r().back()/(f_m.back()*std::sqrt(mesh_m.dr().back()));
		}

		mesh_i = mesh_m.begin();
		g = g_m.begin(), f = f_m.begin();
		for(; mesh_i != mesh_m.end(); g++, f++, mesh_i++){
			*g *= std::sqrt(mesh_i->dr())/mesh_i->r();
			*f *= std::sqrt(mesh_i->dr())/mesh_i->r();
		}
		g = g_m.begin(), f = f_m.begin();
		for(auto gl = left_init_g_m.begin(); gl != left_init_g_m.end(); g++, gl++){
			*g = *gl;
		}
		for(auto fl = left_init_f_m.begin(); fl != left_init_f_m.end(); f++, fl++){
			*f = *fl;
		}
		mesh_i = mesh_m.begin();
		g = g_m.begin(), f = f_m.begin();
		for(; mesh_i != mesh_m.end(); g++, f++, mesh_i++){
			*g *= scale_g;
			*f *= scale_f;
		}

    }

	virtual double norm()
	 const noexcept
    {
        double res = simpson_integral<double>(mesh_m,
            [&](const double& i)
            {
                return GSL::pow_int(g_m[static_cast<size_t>(i)], 2) +
                       GSL::pow_int(f_m[static_cast<size_t>(i)], 2);
            }
        );
        return res;
    }

    /*!*************************************************************************
    * Normalize the solution.
    ***************************************************************************/
    void normalize(const double N = 1)
	 noexcept
    {
        double n = this->norm();

        std::transform(g_m.begin(), g_m.end(), g_m.begin(),
            [=](const double val)
            {
                return val*std::sqrt(N/n);
            });
            std::transform(f_m.begin(), f_m.end(), f_m.begin(),
                [=](const double val)
                {
                    return val*std::sqrt(N/n);
            });
    }

};

#endif // DIRAC_H
