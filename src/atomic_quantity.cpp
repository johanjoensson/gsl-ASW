#include <atomic_quantity.h>
#include <utils.h>
#include <GSLpp/basic_math.h>
#include <iomanip>
#include <numeric>
#include <poisson.h>
#if __cplusplus >= 201703L
#include <execution>
#endif

/*
double Atomic_quantity::operator()(GSL::Vector::Const_View r)
{
    double res = 0.;
    GSL::Vector ri(3);
    size_t t = 1;
    for(const auto& site : cr.sites()){
        size_t at_i = cr.atom_index(site);
        Atom at{cr.atom(site)};
        ri = r - site.pos();
        if(ri.norm() <= at.AS()){
            t = 1;

            // Linear interpolation to find value at point ri
            while(ri.norm() > at_meshes[at_i].r(t) && t < at_meshes[at_i].size()){
                t++;
            }
            if(t < at_meshes[at_i].size()){
                res += lerp(ri.norm(), at_meshes[at_i].r(t-1), at_meshes[at_i].r(t), val[at_i][t-1], val[at_i][t]);
            }
        }
    }

    return res;
}
*/

Potential::Potential(const std::vector<Exponential_mesh<1, double>>& at_meshes, const std::vector<size_t>& zs, const std::vector<size_t>& l_max,
         const Density& rho, Xc_func xcf, bool spinpol)
        : Atomic_quantity(
            ([&](const std::vector<Exponential_mesh<1, double>>& meshes)
            {
                std::vector<size_t> mesh_lengths(meshes.size());
                std::transform(meshes.begin(), meshes.end(), mesh_lengths.begin(),
                    [](const Exponential_mesh<1, double>& m){ return m.size();}
                );
                return mesh_lengths;
            }(at_meshes))
            , l_max, spinpol), atomic_m(0), hartree_m(0), exchange_correlation_pot_m(0),
          exchange_correlation_energy_m(0),
          xc_fun_m(xcf)
{
    size_t len = std::accumulate(std::cbegin(at_meshes), std::cend(at_meshes), 0,
        [](const size_t acc, const Exponential_mesh<1, double> m){ return acc + m.size();}
    );

    atomic_m = GSL::Block(len);
    GSL::Vector::View(atomic_m).set_zero();
/*    
    len = 0;
    for(const auto& at_mesh : at_meshes){
        for(size_t l : l_max){
            std::cout << " lmax = " << l << "\n";
            len += (l + 1)*at_mesh.size();
        }
    }
*/
    len =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::cbegin(at_meshes), std::cend(at_meshes), std::cbegin(l_max), 0,
        [](const size_t acc, const size_t v){ return acc + v;},
        [](const Exponential_mesh<1, double>& m, const size_t l){return m.size()*(l + 1);}
    );

    hartree_m = GSL::Block(len);
    GSL::Vector::View(hartree_m).set_zero();
    exchange_correlation_energy_m = GSL::Block(len);
    GSL::Vector::View(exchange_correlation_energy_m).set_zero();

    if(spinpol){
        std::cout << "SPINPOL!!!\n";
        len *= 2;
    }

    exchange_correlation_pot_m = GSL::Block(len);
    GSL::Vector::View(exchange_correlation_pot_m).set_zero();

    // std::transform(at_meshes.cbegin(), at_meshes.cend(), std::back_inserter(mesh_lengths_m),
    //     [](const Exponential_mesh<1, double>& m)
    //     {
    //         return m.size();
    //     }
    // );
    initial_pot(at_meshes, zs, rho);
}


GSL::Vector::View Potential::atomic(const size_t i)
{
    size_t offset = std::accumulate(std::cbegin(mesh_lengths_m), std::cbegin(mesh_lengths_m) + i, 0,
        [](const size_t acc, const size_t len){ return acc + len;}
    );

    return GSL::Vector::View(this->atomic_m, offset, mesh_lengths_m[i]);
}

GSL::Vector::Const_View Potential::atomic(const size_t i) const
{
    size_t offset = std::accumulate(std::cbegin(mesh_lengths_m), std::cbegin(mesh_lengths_m) + i, 0,
        [](const size_t acc, const size_t len){ return acc + len;}
    );

    return GSL::Vector::Const_View(atomic_m, offset, mesh_lengths_m[i]);
}

GSL::Matrix::View Potential::Hartree(const size_t i)
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::cbegin(mesh_lengths_m), std::cbegin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    return GSL::Matrix::View(this->hartree_m, offset, lmax_m[i] + 1, mesh_lengths_m[i]);
}

GSL::Matrix::Const_View Potential::Hartree(const size_t i) const
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::begin(mesh_lengths_m), std::begin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    return GSL::Matrix::Const_View(this->hartree_m, offset, lmax_m[i] + 1, mesh_lengths_m[i]);
}

GSL::Vector::View Potential::Hartree(const size_t i, const size_t l)
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::cbegin(mesh_lengths_m), std::cbegin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    return GSL::Vector::View(this->hartree_m, offset + l*mesh_lengths_m[i], mesh_lengths_m[i]);
}

GSL::Vector::Const_View Potential::Hartree(const size_t i, const size_t l) const
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::begin(mesh_lengths_m), std::begin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    return GSL::Vector::Const_View(this->hartree_m, offset + l*mesh_lengths_m[i], mesh_lengths_m[i]);
}

GSL::Matrix::View Potential::Vxc(const size_t i)
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::cbegin(mesh_lengths_m), std::cbegin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}
    return GSL::Matrix::View(this->exchange_correlation_pot_m, spin_scale*offset, lmax_m[i] + 1, spin_scale*mesh_lengths_m[i]);
}

GSL::Matrix::Const_View Potential::Vxc(const size_t i) const
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::begin(mesh_lengths_m), std::begin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}
    return GSL::Matrix::Const_View(this->exchange_correlation_pot_m, spin_scale*offset, lmax_m[i] + 1, spin_scale*mesh_lengths_m[i]);
}

GSL::Vector::View Potential::Vxc(const size_t i, const size_t l)
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::cbegin(mesh_lengths_m), std::cbegin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}
    return GSL::Vector::View(this->exchange_correlation_pot_m, spin_scale*(offset + l*mesh_lengths_m[i]), spin_scale*mesh_lengths_m[i]);
}

GSL::Vector::Const_View Potential::Vxc(const size_t i, const size_t l) const
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::begin(mesh_lengths_m), std::begin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}
    return GSL::Vector::Const_View(this->exchange_correlation_pot_m, spin_scale*(offset + l*mesh_lengths_m[i]), spin_scale*mesh_lengths_m[i]);
}

GSL::Matrix::View Potential::exc(const size_t i)
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::cbegin(mesh_lengths_m), std::cbegin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}
    return GSL::Matrix::View(this->exchange_correlation_energy_m, spin_scale*offset, lmax_m[i] + 1, spin_scale*mesh_lengths_m[i]);
}

GSL::Matrix::Const_View Potential::exc(const size_t i) const
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::begin(mesh_lengths_m), std::begin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}
    return GSL::Matrix::Const_View(this->exchange_correlation_energy_m, spin_scale*offset, lmax_m[i] + 1, spin_scale*mesh_lengths_m[i]);
}

GSL::Vector::View Potential::exc(const size_t i, const size_t l)
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::cbegin(mesh_lengths_m), std::cbegin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}
    return GSL::Vector::View(this->exchange_correlation_energy_m, spin_scale*(offset + l*mesh_lengths_m[i]), spin_scale*mesh_lengths_m[i]);
}

GSL::Vector::Const_View Potential::exc(const size_t i, const size_t l) const
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::begin(mesh_lengths_m), std::begin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}
    return GSL::Vector::Const_View(this->exchange_correlation_energy_m, spin_scale*(offset + l*mesh_lengths_m[i]), spin_scale*mesh_lengths_m[i]);

}
/*
Potential::Potential(const Crystal_t<3, Atom>& crystal, const std::vector<Exponential_mesh<1, double>>& at_meshes_n,
     int xc_func_id, std::function<double(const size_t, const double)> atomic_potential)
 : Atomic_quantity(crystal, at_meshes_n), electrostatic(crystal.atoms().size()), exchange_correlation(crystal.atoms().size()),
  at_pot(atomic_potential), MT_0(), xc_fun(xc_func_id)
{}
*/


double atomic_potential(const size_t Z, const double r)
{
    return -2.*static_cast<double>(Z)/r;
}

double mod_at_pot(const size_t Z, const double r)
{
    double alpha = 0.53625;
    double x = r*std::cbrt(Z)/0.88534;
    double phi = 1./((1 + alpha*x)*(1 + alpha*x));
    return -(1 + static_cast<double>(Z - 1)*phi)/r;
}

void Potential::initial_pot(const std::vector<Exponential_mesh<1, double>>& at_meshes, const std::vector<size_t>& zs, const Density& rho)
{
    calc_atomic(at_meshes, zs);
    calc_pot(at_meshes, rho);
}

void Potential::calc_atomic(const std::vector<Exponential_mesh<1, double>>& at_meshes, const std::vector<size_t>& zs)
{
    for(size_t at_i = 0; at_i < at_meshes.size(); at_i++){
        std::transform(std::cbegin(at_meshes[at_i]), std::cend(at_meshes[at_i]), this->atomic(at_i).begin(),
        [&](const Exponential_mesh<1, double>::mesh_point& p)
        {
            return atomic_potential(zs[at_i], p.r());
        });
        if(std::abs(at_meshes[at_i].r(0)) < 1e-16){
            this->atomic(at_i).front() = 0;
        }
    }
}

void Potential::calc_Hartree(const std::vector<Exponential_mesh<1, double>>& at_meshes, const Density& rho)
{
    Poisson_equation p;
    for(size_t at_i = 0; at_i < at_meshes.size(); at_i++){
        const GSL::Vector rho_tot_i(rho.tot(at_i));
        for( size_t l = 0; l < lmax_m[at_i]; l++){
            p.solve(at_meshes[at_i], {rho_tot_i.gsl_data()->data, rho_tot_i.gsl_data()->data + rho_tot_i.size()}, l);
            this->Hartree(at_i, l).copy(GSL::Vector::Const_View{p.val().data(), p.val().size()});
        }
    }
}

void Potential::calc_XC(const Density& rho)
{
    for(size_t at_i = 0; at_i < lmax_m.size(); at_i++){
        for( size_t l = 0; l < lmax_m[at_i]; l++){
            this->Vxc(at_i, l).swap(this->xc_fun_m.vxc(rho.valence(at_i, l)));
            this->exc(at_i, l).swap(this->xc_fun_m.exc(rho.valence(at_i, l)));
        }
    }
}

void Potential::calc_pot(const std::vector<Exponential_mesh<1, double>>& at_meshes, const Density& rho)
{
    calc_Hartree(at_meshes, rho);
    calc_XC(rho);
}

const GSL::Vector Potential::tot(size_t at, size_t l) const
{
    return this->atomic(at) + this->Hartree(at, l) + this->Vxc(at, l);
}

/*
double Potential::Xi0(const Site_t<3>& j, const double r)
{
    double res = 0, a = 0;
    for(const auto& site : cr.sites()){
        if(site != j){
            a = (j.pos() - site.pos()).norm();
            res += (a + r)/a*atomic_potential(cr.atom(site).Z(), a + r);
        }
    }
    return res;
}

void Potential::initial_pot(double vol)
{
    size_t nel = 0;
    for(const auto& at : cr.atoms()){
        nel += at.Z();
    }
    std::vector<double> rho;
    for(size_t i = 0; i < cr.atoms().size(); i++){
        electrostatic[i] = std::vector<double>(at_meshes[i].size());
        val[i] = std::vector<double>(at_meshes[i].size());
        // exchange_correlation[i] = std::vector<double>(at_meshes[i].size(), 0.);
        rho = std::vector<double>(at_meshes[i].size(), static_cast<double>(nel)/vol);
        exchange_correlation[i] = xc_fun.exc(rho);
        for(size_t j = 0; j < at_meshes[i].size(); j++){
            electrostatic[i][j] =
            at_pot(cr.atoms(i).Z(), at_meshes[i].r(j));
        }
    }
    std::vector<bool> done(cr.atoms().size(), false);
    double r = 0., r1 = 0., r2 = 0., drx = 0., drx1 = 0., drx2 = 0.;
    for(const auto& site : cr.sites()){
        size_t at_i = cr.atom_index(site);
        if(done[at_i]){
            continue;
        }
        std::vector<double> alpha(at_meshes[at_i].size());
        alpha[0] = 0.;
        alpha[1] = 1./6 * (4*Xi0(site, at_meshes[at_i].r(0))*at_meshes[at_i].dr(static_cast<size_t>(0))
    + (Xi0(site, at_meshes[at_i].r(1)) + Xi0(site, -at_meshes[at_i].r(1)))*at_meshes[at_i].dr(static_cast<size_t>(1)));
        electrostatic[at_i][0] += Xi0(site, 0);
        electrostatic[at_i][1] += alpha[1] * 1/at_meshes[at_i].r(1);
        for(size_t ri = 2; ri < at_meshes[at_i].size(); ri++){
            r = at_meshes[at_i].r(ri);
            r1 = at_meshes[at_i].r(ri-1);
            r2 = at_meshes[at_i].r(ri-2);
            drx = at_meshes[at_i].dr(ri);
            drx1 = at_meshes[at_i].dr(ri-1);
            drx2 = at_meshes[at_i].dr(ri-2);
            alpha[ri] = alpha[ri - 2];
            alpha[ri] += 1./6 * drx2 * (Xi0(site, r2) + Xi0(site, -r2));
            alpha[ri] += 4./6 * drx1 * (Xi0(site, r1) + Xi0(site, -r1));
            alpha[ri] += 1./6 * drx * (Xi0(site, r) + Xi0(site, -r));

            electrostatic[at_i][ri] += alpha[ri]/r;
        }
        done[at_i] = true;
        std::cout << "Done setting up electrostatic potential of atom " << at_i << "\n";
    }

    for(size_t i = 0; i < cr.atoms().size(); i++){
        for(size_t j = 0; j < at_meshes[i].size(); j++){
            val[i][j] = electrostatic[i][j]  + exchange_correlation[i][j] ;
        }
    }

    // Calculate MT_0 as average potential over all atomic spheres
    double MT_0_s = 0, areas = 0;
    for(size_t i = 0; i < cr.atoms().size(); i++){
        MT_0_s += val[i].back()*GSL::pow_int(cr.atoms(i).AS(), 2);
        areas += GSL::pow_int(cr.atoms(i).AS(), 2);
        std::cout << std::fixed << std::setprecision(12) << "MT_0 contribution = " << val[i].back()*GSL::pow_int(cr.atoms(i).AS(), 2) << "\n";
    }
    this->MT_0 = MT_0_s/areas;

	std::cout << "MT0 = " << this->MT_0 << std::endl;


    // Make potential relative to MT_0
    for(size_t i = 0; i < cr.atoms().size(); i++){
        for(size_t j = 0; j < at_meshes[i].size(); j++){
            val[i][j] -= this->MT_0;
        }
    }

}
*/


Density::Density(const std::vector<Exponential_mesh<1, double>>& at_meshes, const std::vector<size_t>& l_max,
                 bool spinpol)
        : Atomic_quantity(
            ([&](const std::vector<Exponential_mesh<1, double>>& meshes)
            {
                std::vector<size_t> mesh_lengths(meshes.size());
                std::transform(meshes.begin(), meshes.end(), mesh_lengths.begin(),
                    [](const Exponential_mesh<1, double>& m){ return m.size();}
                );
                return mesh_lengths;
            }(at_meshes))
            , l_max, spinpol), valence_m(0), core_m(0)
{
    size_t len =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::cbegin(at_meshes), std::cend(at_meshes), std::cbegin(l_max), 0,
        [](const size_t acc, const size_t v){ return acc + v;},
        [](const Exponential_mesh<1, double>& m, const size_t l){return m.size()*(l + 1);}
    );
    if(spinpol_m){
        len *= 2;
    }
    valence_m = GSL::Block(len);
    len = std::accumulate( std::cbegin(at_meshes), std::cend(at_meshes), 0,
        [](const size_t acc, const Exponential_mesh<1, double>& m)
        {
            return acc + m.size();
        }
    );
    if(spinpol_m){
        len *= 2;
    }
    core_m = GSL::Block(len);
}

Density::Density(const std::vector<size_t>& mesh_lengths, const std::vector<size_t> lmax,
    bool spinpol)
    : Atomic_quantity(mesh_lengths, lmax, spinpol), valence_m(0), core_m(0)
{
    size_t len =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::cbegin(mesh_lengths), std::cend(mesh_lengths), std::cbegin(lmax), 0,
        [](const size_t acc, const size_t v){ return acc + v;},
        [](const size_t m, const size_t l){return m*(l + 1);}
    );
    if(spinpol_m){
        len *= 2;
    }
    valence_m = GSL::Block(len);
    len = std::accumulate( std::cbegin(mesh_lengths), std::cend(mesh_lengths), 0,
        [](const size_t acc, const size_t m)
        {
            return acc + m;
        }
    );
    if(spinpol_m){
        len *= 2;
    }
    core_m = GSL::Block(len);
}

GSL::Matrix::View Density::valence(size_t i)
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::begin(mesh_lengths_m), std::begin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}

    return GSL::Matrix::View(this->valence_m, spin_scale*offset,  lmax_m[i] + 1, spin_scale*mesh_lengths_m[i]);
}

GSL::Matrix::Const_View Density::valence(size_t i) const
{
    size_t offset =
#if __cplusplus >= 201703L
    std::transform_reduce(std::execution::par_unseq,
#else
    std::inner_product(
#endif
        std::begin(mesh_lengths_m), std::begin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}

    return GSL::Matrix::Const_View(this->valence_m, spin_scale*offset,  lmax_m[i] + 1, spin_scale*mesh_lengths_m[i]);
}

GSL::Vector::View Density::valence(size_t i, size_t l)
{
    size_t offset =
    #if __cplusplus >= 201703L
        std::transform_reduce(std::execution::par_unseq,
    #else
        std::inner_product(
    #endif
            std::begin(mesh_lengths_m), std::begin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
            [](const size_t acc, const size_t len){ return acc + len;},
            [](const size_t len, const size_t l) {return len*(l + 1);}
        );

        size_t spin_scale = 1;
        if (spinpol_m){spin_scale = 2;}

        return GSL::Vector::View(this->valence_m, spin_scale*(offset + l*mesh_lengths_m[i]), spin_scale*mesh_lengths_m[i]);
}

GSL::Vector::Const_View Density::valence(size_t i, size_t l) const
{
    size_t offset =
    #if __cplusplus >= 201703L
        std::transform_reduce(std::execution::par_unseq,
    #else
        std::inner_product(
    #endif
    std::begin(mesh_lengths_m), std::begin(mesh_lengths_m)  + i, std::cbegin(lmax_m), 0,
        [](const size_t acc, const size_t len){ return acc + len;},
        [](const size_t len, const size_t l) {return len*(l + 1);}
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}

    return GSL::Vector::Const_View(this->valence_m, spin_scale*(offset + l*mesh_lengths_m[i]), spin_scale*mesh_lengths_m[i]);
}


GSL::Vector::View Density::core(size_t i)
{
    size_t offset = std::accumulate( std::cbegin(mesh_lengths_m), std::cbegin(mesh_lengths_m) + i, 0,
        [](const size_t acc, size_t len)
        {
            return acc + len;
        }
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}

    return GSL::Vector::View(this->core_m, spin_scale*offset, spin_scale*mesh_lengths_m[i]);
}

GSL::Vector::Const_View Density::core(size_t i) const
{
    size_t offset = std::accumulate( std::cbegin(mesh_lengths_m), std::begin(mesh_lengths_m) + i, 0,
        [](const size_t acc, size_t len)
        {
            return acc + len;
        }
    );

    size_t spin_scale = 1;
    if (spinpol_m){spin_scale = 2;}

    return GSL::Vector::Const_View(this->core_m, spin_scale*offset, spin_scale*mesh_lengths_m[i]);
}

GSL::Vector Density::tot(size_t i) const
{
    GSL::Vector v(this->mesh_lengths_m[i], 0);
    for(size_t l = 0; l < lmax_m[i]; l++){
        v += this->valence(i, l);
    }

    GSL::Vector::Const_View c(this->core(i));
    if(spinpol_m){
        GSL::Vector::Const_View vup = v.cview(0, v.size(), 2);
        GSL::Vector::Const_View vdown = v.cview(1, v.size(), 2);
        GSL::Vector::Const_View cup = c.cview(0, c.size(), 2);
        GSL::Vector::Const_View cdown = c.cview(1, c.size(), 2);

        return vup + vdown + cup + cdown;
    } else {
        return v + c;
    }
}
