#include <dirac.h>
#include <numerical-mesh.h>
#include <vector>
#include <algorithm>
#include <limits>

int main()
{
    const double c = 2*137.035999174;
    const double alpha = 2./c;

    unsigned int z = 1;

    auto V = [&] (const double r) {
        return -2./r;
    };

    Exponential_mesh<1, double> m{0, 250, 0.001, 10001};
    std::vector<double> pot;
    std::transform(m.begin(), m.end(), std::back_inserter(pot),
    [=](const Exponential_mesh<1, double>::mesh_point p)
    {
        return V(p.r());
    });
    // pot.front() = 0;
    pot[0] = pot[1]*1000;

    std::cout << "epsilon = " << std::numeric_limits<double>::epsilon() << "\n";


    int n = 5, l = 2, s = -1;
    int kappa = -l - 1;
    if(s < 0){
        kappa = l;
    }
    auto beta = sqrt(GSL::pow_int(kappa, 2) - GSL::pow_int(2/c, 2));
    auto rb = GSL::exp(beta*GSL::log(m.r(1))).val;
    double lambda = sqrt(-0.5 - GSL::pow_int(-0.5/c, 2));
    std::array<std::vector<double>, 2> left_init {
                                        std::vector<double>{0, rb},
                                        std::vector<double>{0, rb*c*(beta + kappa)}
                                    };
    std::array<std::vector<double>, 2> right_init {
                                            std::vector<double>{2*std::numeric_limits<double>::epsilon()},
                                            std::vector<double>{-2*std::numeric_limits<double>::epsilon()}
                                        };

    Radial_Dirac_Equation dirac{m, pot, kappa, left_init, right_init, std::numeric_limits<double>::epsilon()};

    dirac.solve(n - l - 1, -1./GSL::pow_int(n, 2), -1.1, 0);

    double exact_energy = -1./GSL::pow_int(n, 2)*(
         1 + GSL::pow_int(z*alpha, 2)/n*(1./std::abs(kappa) - 3./(4*n)));
    std::cout << "Energy = " << dirac.energy() << "(Ry), error ~ " <<
        std::abs(dirac.energy() - exact_energy)*1e6 << "(\u03BCRy) \n";
    std::cout << "N = " << dirac.norm() << "\n";
    dirac.normalize();
    std::ofstream os;
    os.open("dirac_test.dat", std::ios::out);
    for(size_t i = 0; i < m.size(); i++){
        auto dr = m.dr(i);
        os << m.r(i) << "  " << pot[i] << "  " << dirac.g(i)*sqrt(dr) << "  " << dirac.f(i)*sqrt(dr) << "\n";
    }
    os.close();
    return 0;
}