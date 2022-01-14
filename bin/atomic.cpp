#include <dirac.h>
#include <schroedinger.h>
#include <numerical-mesh.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <poisson.h>
#include <atomic_quantity.h>
#include <GSLpp/error.h>
#include <utils.h>
#include <xc_func.h>
#include <mixer.h>

double etot(const Exponential_mesh<1, double>& m, const double e, const Density& rho,
    const Potential& pot)
{
    double Eks = e;
    double dExc = 4*M_PI*simpson_integral<double>(m,
        [&](const Mesh_base<1, double>::mesh_point& p)
        {
            size_t id = static_cast<size_t>(p.i());
            return pot.Vxc(0, 0)[id]*rho.tot(0)[id]*p.r2();
        }
    );

    double EH = 2*M_PI*simpson_integral<double>(m,
        [&](const Mesh_base<1, double>::mesh_point& p)
        {
            size_t id = static_cast<size_t>(p.i());
            return pot.Hartree(0,0)[id]*rho.tot(0)[id]*p.r2();
        }
    );

    double Exc = 4*M_PI*simpson_integral<double>(m,
        [&](const Mesh_base<1, double>::mesh_point& p)
        {
            size_t id = static_cast<size_t>(p.i());
            return pot.exc(0, 0)[id]*rho.tot(0)[id]*p.r2();
        }
    );

    std::cout << std::string(80, '=') << "\n";
    std::cout << "EKS = " << Eks << "\n";
    std::cout << "EH = " << EH << "\n";
    std::cout << "Exc = " << Exc << "\n";
    std::cout << "dExc/drho = " << dExc << "\n\n";

    std::cout << "Etot = " << Eks - EH + Exc - dExc  << "\n";
    std::cout << std::string(80, '=') << "\n";

    return Eks - EH + Exc - dExc;
}

int main()
{
    // GSL::Error_handler e_handler;
    // e_handler.off();



    /***************************************************************************
    * NIST reference value                                                     *
    ***************************************************************************/
    const double alpha = 7.2973525693e-3;
    const double c = 2./alpha;
    double mc2 = c*c/2;

    unsigned int z = 2;

    auto at_pot = [&] (const double r) {
        return -2.*z/r;
    };

    std::vector<std::tuple<ls, double, double>> states;
    ls nlm{1, 0, Spin::Up, -1};
    for(size_t ne = 0; ne < z; ne ++, nlm++){
        double e_SH = -static_cast<double>(z*z)/(nlm.n*nlm.n);
        states.push_back({nlm, 1, e_SH});
    }

    double r_min = 0, r_max = 8., a = 1e-2;
    size_t N = 1001;
    Exponential_mesh<1, double> m(r_min, r_max, a, N);


    // GSL::Vector Vtot(m.size()), Vat(m.size()), VH(m.size()),
    //             Vxc(m.size()), exc(m.size());
    // std::transform(m.begin(), m.end(), Vat.view().begin(),
    //     [=](const Exponential_mesh<1, double>::mesh_point& p)
    //     {
    //         return at_pot(p.r());
    //     });
    // Vat.view().front() = 0;

    Density rho_in({m}, {0}); //, rho_out({m}, {0});
    rho_in.core(0).set_zero();
    rho_in.valence(0,0).set_all(3.*z/(4*M_PI*GSL::pow_int(r_max, 3)));
    Potential pot({m}, {z}, {0}, rho_in, Xc_func(43), false);

    std::ofstream os;
    os.open("dirac_test.dat", std::ios::out);


    double alpha_mix = 0.8;
    Linear_mixer<Density> mixer(alpha_mix);
    double previous_etot = 0, e_new = 0;
    double e_diff = 1;
    for(auto it = 0; it < 1 && e_diff > 1e-10; it++){

        Density rho_out({m}, {0});
        double e_sum = 0;
        rho_out.valence(0, 0).set_zero();
        for(auto& state : states){
            ls nlm = std::get<0>(state);
            int kappa = nlm.kappa;
            double occ = std::get<1>(state);
            double& en = std::get<2>(state);
            auto beta = std::sqrt(GSL::pow_int(kappa, 2) - GSL::pow_int(z*alpha, 2));

            double g0, g1;
            try {
                g1 = 1;
                g0 = GSL::exp((beta - 1)*GSL::log(1 + (-m.dr(1) + m.d2r(1)/2)/m.r(1))).val;
            }catch (const std::runtime_error& e){
                std::cerr << e.what() << "\n";
                return 1;
            }
            double f0, f1;
            if(kappa > 0){
                f0 = g0*(kappa + beta)/(z*alpha);
                f1 = g1*(kappa + beta)/(z*alpha);
            }else{
                f0 = g0*(z*alpha)/(beta - kappa);
                f1 = g1*(z*alpha)/(beta - kappa);
            }

            std::array<std::vector<double>, 2> left_init {
                                                std::vector<double>{g0, g1},
                                                std::vector<double>{f0, f1}
                                            };
            std::array<std::vector<double>, 2> right_init {
                                                    std::vector<double>{std::numeric_limits<double>::epsilon()},
                                                    std::vector<double>{-std::numeric_limits<double>::epsilon()}
                                                };

            GSL::Vector vv(pot.tot(0,0).clone());
            std::vector<double> v(vv.gsl_data()->data, vv.gsl_data()->data + vv.size());
            Radial_Dirac_Equation dirac{m, v, kappa, left_init, right_init, 1e-12};

            dirac.solve(nlm.n - nlm.l - 1, en, en*1.5, v.back(), true);

            en = dirac.energy();
            dirac.normalize(1);
            auto g_v = dirac.g(), f_v = dirac.f();
            GSL::Vector::Const_View g(g_v.data(), g_v.size()), f(f_v.data(), f_v.size());
            rho_out.valence(0, 0) += occ/(4*M_PI)*(g*g + f*f);

            if (it == 0){
                double exact_energy = mc2/std::sqrt(
                    1 + GSL::pow_int(
                            z*alpha/(nlm.n - std::abs(kappa) + beta)
                        , 2)
                ) - mc2;
                std::cout << "\tE" << nlm << ":" << it << " = " << en << "\n";
                std::cout << "\t|e - e_exact| = " << std::abs(en - exact_energy) << "\n";
            }
            e_sum += occ*en;


            std::ofstream orbitals;
            orbitals.open("orbitals.dat", std::ios::app);
            for(size_t i = 0; i < m.size(); i++){
                orbitals << m.r(i) << " " << dirac.g(i) << " " << dirac.f(i) << "\n";
            }
            orbitals << "\n";
            orbitals.close();
        }

        /*
        for(size_t i = 0; i < m.size(); i++){
            auto dr = m.dr(i);
            os << m.r(i) << "  " << pot.atomic(0)[i] << "  " << pot.Hartree(0,0)[i] << "  " << pot.Vxc(0,0)[i]
               << "  " << rho_out.tot(0)[i] << "\n";
        }
        os << "\n\n";
        */
        e_new = etot(m, e_sum, rho_out, pot);
        e_diff = std::abs(e_new - previous_etot);
        previous_etot = e_new;
        std::cout << "\tEdiff = " << e_diff << "\n\n";

        rho_in = mixer.mix(rho_in, rho_out);
    }
    os.close();

    // std::cout << "exact single particle energy is " << exact_energy << "\n";
    std::cout << "Total energy = " << e_new << " Ry\n";
    // std::cout << "reference LDA value for total energy is -5.6696  Ry\n";


    return 0;
}
