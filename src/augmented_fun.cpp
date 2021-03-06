#include "augmented_fun.h"
#include "spherical_fun.h"
#include "schroedinger.h"
#include "envelope_fun.h"

#include <numeric>
#include <fstream>
#include <exception>


Augmented_function::Augmented_function( const lm l_n, const double kappa_n, const spin s_n,
    const Logarithmic_mesh& mesh_n)
 : En_m(0), l_m(l_n), kappa_m(kappa_n), s_m(s_n), mesh_m(mesh_n),
 val_m(mesh_n.size())
{}

double Augmented_function::operator()(const GSL::Vector& r) const
{
    double res = 0, radius = mesh_m.r_back();
    size_t t = 1;
    if(r.norm<double>() <= radius){
        while(mesh().r(t) < r.norm<double>() && t < mesh().size()){
            t++;
        }
        // r[t - 1] < |ri| and r[t] > |ri|
        if(t < mesh().size()){
            res = lerp(r.norm<double>(), mesh().r(t - 1), mesh().r(t),
            val_m[t - 1]*std::sqrt(mesh().drx(t-1))*mesh().drx(t-1), val_m[t]*std::sqrt(mesh().drx(t))*mesh().drx(t));
        }
    }
    return res;
}

double augmented_integral(const Augmented_function &a, const Augmented_function &b)
{
    if(std::abs(a.kappa() - b.kappa()) < 1e-10){
	    Logarithmic_mesh mesh = a.mesh();
	    std::vector<double> integrand(a.val().size(), 0);

	    for(size_t i = 0; i < integrand.size(); i++){
		    integrand[i] = a.val()[i]*b.val()[i];
	    }
	    return mesh.integrate_simpson(integrand);
    }else{
	    Envelope_Hankel h1(a.l(), a.kappa()), h2(b.l(), b.kappa());
	    return 1./(a.En() - b.En()) * wronskian(h1, h2, a.mesh().r_back());
    }
}

bool operator==(const Augmented_function &a, const Augmented_function &b)
{
    return a.l() == b.l();
}

bool operator!=(const Augmented_function &a, const Augmented_function &b)
{
    return !(a == b);
}

Augmented_Hankel::Augmented_Hankel(const lm l_n, const double kappa_n, const spin s_n,
    const Logarithmic_mesh& mesh_n)
 : Augmented_function(l_n, kappa_n, s_n, mesh_n)
{}

void Augmented_Hankel::update(std::vector<double>& v, const double en
    , const bool core)
{
    EH() = en;
    size_t nodes = static_cast<size_t>(std::max(0, l_m.n - l_m.l - 1));
    size_t last = mesh_m.size() - 1, lastbutone = mesh_m.size() - 2;

    Envelope_Hankel H(l_m, kappa_m);

    std::vector<double> l_init = {
        0,
        GSL::pow_int(mesh_m.r(1), l_m.l + 1),
        GSL::pow_int(mesh_m.r(2), l_m.l + 1)};
    std::vector<double> r_init;
    if(core){
        r_init = {1e-8, 0.};
    }else{
        r_init = {mesh_m.r(lastbutone)*H.barred_fun(mesh_m.r(lastbutone)),
                  mesh_m.r(last)*H.barred_fun(mesh_m.r(last))};
    }
    Radial_Schroedinger_Equation_Central_Potential se(v, static_cast<size_t>(l_m.l), l_init, r_init, mesh_m, 1e-10);
    se.solve(nodes, EH());
    if(core){
        se.normalize();
    }else{
        double scale = 1;
        if(std::abs(se.psi().back()) > 1e-12){
            scale = r_init.back()/se.psi().back();
        }
        for(size_t i = 0; i < mesh_m.size(); i++){
            se.psi()[i] *= scale;
        }
    }
    val_m = se.psi();
    EH() = se.e();
#ifdef DEBUG
    std::ofstream out_file("check_Hankel.dat", std::ios::app);
    for(size_t i = 0; i < mesh_m.size(); i++){
        out_file << mesh_m.r(i) << " " << val_m[i] << " " << v[i] << "\n";
    }
    out_file << "\n\n";
#endif //DEBUG
}

Augmented_Bessel::Augmented_Bessel(const lm l_n, const double kappa_n, const spin s_n,
    const Logarithmic_mesh& mesh_n)
 : Augmented_function(l_n, kappa_n, s_n, mesh_n)
{}

void Augmented_Bessel::update(std::vector<double>& v, const double en
    , const bool core)
{
    EJ() = en;
    size_t nodes = static_cast<size_t>(std::max(0, l_m.n - l_m.l - 1));
    size_t last = mesh_m.size() - 1, lastbutone = mesh_m.size() - 2;
    Envelope_Bessel I(l_m, kappa_m);

    if(!core){
        std::vector<double> l_init = {
            0,
            GSL::pow_int(mesh_m.r(1), l_m.l+1),
            GSL::pow_int(mesh_m.r(2), l_m.l+1)};

        std::vector<double> r_init = {
            I.barred_fun(mesh_m.r(lastbutone))
		        *mesh_m.r(lastbutone),
		    I.barred_fun(mesh_m.r(last))
		        *mesh_m.r(last)};

        Radial_Schroedinger_Equation_Central_Potential se(v, static_cast<size_t>(l_m.l), l_init, r_init, mesh_m, 1e-10);
        se.solve(nodes, EJ());
        double scale = r_init.back()/se.psi().back();
        for(size_t i = 0; i < mesh_m.size(); i++){
            se.psi()[i] *= scale;
        }
        val_m = se.psi();
        EJ() = se.e();
    }else{
        EJ() = 0.;
        val_m = std::vector<double>(val_m.size(), 0.);
    }
#ifdef DEBUG
    std::ofstream out_file("check_Bessel.dat", std::ios::app);
    for(size_t i = 0; i < mesh_m.size(); i++){
        out_file << mesh_m.r(i) << " " << val_m[i] << " " << v[i] << "\n";
    }
    out_file << "\n\n";
#endif //DEBUG
}

void Hankel_container::add_function(const Augmented_Hankel& H)
{
    auto cmp = [](const Augmented_Hankel& H1, const Augmented_Hankel& H2)
    {
        if(H1.s() == H2.s()){
            if(H1.kappa() == H2.kappa()){
                return H1.l() < H2.l();
            }else{
                return H1.kappa() < H2.kappa();
            }
        }else{
            return H1.s() < H2.s();
        }
        return false;
    };
    auto it = std::upper_bound(functions.begin(), functions.end(), H, cmp);
    functions.insert(it, H);
}

Augmented_Hankel& Hankel_container::get_function(const lm& l, const double& kappa, const spin& s)
{
    auto cmp = [](const Augmented_Hankel& H1, const Augmented_Hankel& H2)
    {
        if(H1.s() == H2.s()){
            if(H1.kappa() == H2.kappa()){
                return H1.l() < H2.l();
            }else{
                return H1.kappa() < H2.kappa();
            }
        }else{
            return H1.s() < H2.s();
        }
        return false;
    };
    auto it = std::lower_bound(functions.begin(), functions.end(),
        Augmented_Hankel ({l.n, l.l, -l.l}, kappa, s, Logarithmic_mesh()), cmp);
    if(it == functions.end()){
        throw std::runtime_error("Hankel function, " + l.to_string() + ", not found!\n");
    }
    return *it;
}

size_t Hankel_container::get_index(const lm& l, const double& kappa, const spin& s) const
{
    auto cmp = [](const Augmented_Hankel& H1, const Augmented_Hankel& H2)
    {
        if(H1.s() == H2.s()){
            if(H1.kappa() == H2.kappa()){
                return H1.l() < H2.l();
            }else{
                return H1.kappa() < H2.kappa();
            }
        }else{
            return H1.s() < H2.s();
        }
        return false;
    };
    auto it = std::lower_bound(functions.begin(), functions.end(),
            Augmented_Hankel ({l.n, l.l, -l.l}, kappa, s, Logarithmic_mesh()), cmp);
    if(it == functions.end()){
        throw std::runtime_error("Hankel function, " + l.to_string() + ", not found!\n");
    }
    return static_cast<size_t>(std::distance(functions.begin(), it));
}

void Bessel_container::add_function(const Augmented_Bessel& J)
{
    auto cmp = [](const Augmented_Bessel& J1, const Augmented_Bessel& J2)
    {
        if(J1.s() == J2.s()){
            if(J1.kappa() == J2.kappa()){
                return J1.l() < J2.l();
            }else{
                return J1.kappa() < J2.kappa();
            }
        }else{
            return J1.s() < J2.s();
        }
        return false;
    };
    auto it = std::upper_bound(functions.begin(), functions.end(), J, cmp);
    functions.insert(it, J);
}

Augmented_Bessel& Bessel_container::get_function(const lm& l, const double& kappa, const spin& s)
{
    auto cmp = [](const Augmented_Bessel& J1, const Augmented_Bessel& J2)
    {
        if(J1.s() == J2.s()){
            if(J1.kappa() == J2.kappa()){
                return J1.l() < J2.l();
            }else{
                return J1.kappa() < J2.kappa();
            }
        }else{
            return J1.s() < J2.s();
        }
        return false;
    };
    auto it = std::lower_bound(functions.begin(), functions.end(),
        Augmented_Bessel ({l.n, l.l, -l.l}, kappa, s, Logarithmic_mesh()), cmp);
    if(it == functions.end()){
        throw std::runtime_error("Bessel function, " + l.to_string() + ", not found!\n");
    }
    return *it;
}

size_t Bessel_container::get_index(const lm& l, const double& kappa, const spin& s) const
{
    auto cmp = [](const Augmented_Bessel& J1, const Augmented_Bessel& J2)
    {
        if(J1.s() == J2.s()){
            if(J1.kappa() == J2.kappa()){
                return J1.l() < J2.l();
            }else{
                return J1.kappa() < J2.kappa();
            }
        }else{
            return J1.s() < J2.s();
        }
        return false;
    };
    Augmented_Bessel J({l.n, l.l, -l.l}, kappa, s, Logarithmic_mesh());
    auto it = std::lower_bound(functions.begin(), functions.end(),
            Augmented_Bessel ({l.n, l.l, -l.l}, kappa, s, Logarithmic_mesh()), cmp);
    if(it == functions.end()){
        throw std::runtime_error("Bessel function, " + l.to_string() + ", not found!\n");
    }
    return static_cast<size_t>(std::distance(functions.begin(), it));
}

lm Bessel_container::min_lm() const
{
    auto cmp = [](const Augmented_Bessel& J1, const Augmented_Bessel& J2)
        {
            return J1.l() < J2.l();
        };
    auto it = std::min_element(functions.begin(), functions.end(), cmp);
    if(it == functions.end()){
        throw std::runtime_error("Minimum element not found!\n");
    }
    return it->l();

}


lm Bessel_container::max_lm() const
{
    auto cmp = [](const Augmented_Bessel& J1, const Augmented_Bessel& J2)
        {
            return J1.l() < J2.l();
        };
    auto it = std::max_element(functions.begin(), functions.end(), cmp);
    if(it == functions.end()){
        throw std::runtime_error("Maximum element not found!\n");
    }
    return it->l();

}
