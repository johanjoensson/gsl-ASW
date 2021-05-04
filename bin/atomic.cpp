#include <dirac.h>
#include <schroedinger.h>
#include <numerical-mesh.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <poisson.h>
#include <GSLpp/error.h>

int main()
{
    // GSL::Error_handler e_handler;
    // e_handler.off();

    const double c = 2*137.035999174;
    const double alpha = 2./c;

    unsigned int z = 92;

    auto V = [&] (const double r) {
        return -2.*z/r;
    };

    Exponential_mesh<1, double> m{0, 40, 1e-3, 50001};
    std::vector<double> pot;
    std::transform(m.begin(), m.end(), std::back_inserter(pot),
        [=](const Exponential_mesh<1, double>::mesh_point p)
        {
            return V(p.r());
        });

    std::cout << "epsilon = " << std::numeric_limits<double>::epsilon() << "\n";


    int n = 7, l = 0, s = 1;
    int kappa = -l - 1;
    if(s < 0){
        kappa = l;
    }
    auto beta = std::sqrt(GSL::pow_int(kappa, 2) - GSL::pow_int(z*alpha, 2));
    double rb;
    try {
        rb = GSL::exp((beta-1)*GSL::log(m.r(1))).val;
    }catch (const std::runtime_error& e){
        std::cerr << e.what() << "\n";
        return 1;
    }
    auto rb0 = rb*(1 + (beta - 1)/m.r(1)*m.dr(1));
    std::cout << "beta = " << beta << "\n";
    std::cout << "rb = " << rb << "\n";
    std::cout << "rb0 = " << rb0 << "\n";
    std::cout << "f/g = " << c*(beta + kappa)/(z*alpha) << "\n";
    std::array<std::vector<double>, 2> left_init {
                                        std::vector<double>{rb0, rb},
                                        std::vector<double>{rb0*c*(beta + kappa)/(z*alpha), rb*c*(beta + kappa)/(z*alpha)}
                                    };
    std::array<std::vector<double>, 2> right_init {
                                            std::vector<double>{1e-16},
                                            std::vector<double>{-1e-18}
                                        };

    Radial_Dirac_Equation dirac{m, pot, kappa, left_init, right_init, 1e-10};


    double mc2 = c*c/2;
    double exact_energy = mc2/std::sqrt(
        1 + GSL::pow_int(
                z*alpha/(n - std::abs(kappa) + std::sqrt(GSL::pow_int(kappa, 2)
                - GSL::pow_int(z*alpha, 2)))
            , 2)
    ) - mc2;
    std::cout << "Eguess = " << exact_energy << "\n";

    dirac.solve(n - l - 1, exact_energy, -200, -0);

    dirac.normalize();
    std::cout << "Energy = " << dirac.energy() << "(Ry), error ~ " <<
        std::abs(dirac.energy() - exact_energy)*1e6 << "(\u03BCRy) \n";


    std::vector<double> rho;
    auto g = dirac.g(), f = dirac.f();
    std::transform(g.begin(), g.end(), f.begin(), std::back_inserter(rho),
        [](const double g, const double f)
        {
            return g*g + f*f;
        }
    );

    Poisson_equation pos(m, rho, 0);
    pos.solve();
    std::ofstream os;
    os.open("dirac_test.dat", std::ios::out);
    for(size_t i = 0; i < m.size(); i++){
        auto dr = m.dr(i);
        os << m.r(i) << "  " << pot[i] << "  " << dirac.g(i) << "  " <<
              dirac.f(i) << "  " << rho[i] << "  " << pos.val(i) << "\n";
    }
    os.close();

    return 0;
}
