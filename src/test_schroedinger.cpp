#include <schroedinger.h>
#include <vector>
#include <fstream>
#include <GSLpp/special_functions.h>
#include <GSLpp/error.h>
#include <cmath>
#include <algorithm>

std::fstream outfile;

void harmonic_oscillator()
{
    unsigned int len = 2000;
	double l = 12;
	double j = -l/2;
	std::vector<double> v(len, 0);
	for(auto it = v.begin(); it != v.end(); it++, j += l/static_cast<double>(len - 1)){
		*it = GSL::pow_int(j, 2);
	}
	int n = 8;
    auto exact_energy = 2*(n + 0.5);

	auto fac = [](int n_){
		double res = GSL::fact(static_cast<unsigned int>(n_)).val;
		return res;
    };
	std::vector<double> left = {
		GSL::pow_int(-1., n)/std::sqrt(GSL::pow_int(2, n)*fac(n))*1./std::sqrt(std::sqrt(M_PI))*
			GSL::exp(-0.5*GSL::pow_int(-l/2, 2)).val*GSL::hermite_phys(n, -l/2).val,
		GSL::pow_int(-1., n)/std::sqrt(GSL::pow_int(2, n)*fac(n))*1./std::sqrt(std::sqrt(M_PI))*
			GSL::exp(-0.5*GSL::pow_int(-l/2 + 1./static_cast<double>(len-1), 2)).val*GSL::hermite_phys(n, -l/2 + 1./static_cast<double>(len - 1)).val};
	std::vector<double> right = {
		1./std::sqrt(GSL::pow_int(2, n)*fac(n))*1./std::sqrt(std::sqrt(M_PI))*
			GSL::exp(-0.5*GSL::pow_int(l/2 - 1./static_cast<double>(len - 1), 2)).val* GSL::hermite_phys(n, l/2 - 1./(len - 1)).val,
		1./std::sqrt(GSL::pow_int(2, n)*fac(n))*1./std::sqrt(std::sqrt(M_PI))*
			GSL::exp(-0.5*GSL::pow_int(l/2, 2)).val*GSL::hermite_phys(n, l/2).val
	 };
     Linear_mesh<1, double> m(-l/2, l/2, len);
	Schroedinger_Equation se(0, 20, v, left, right, m, 1e-6);
	// Schroedinger_Equation se(0, 50, v, left, right, l/static_cast<double>(len - 1), 1e-6);
	se.solve(static_cast<size_t>(n));
    se.normalize();
    std::cout << "Found energy : " << se.e() << "\n";
    std::cout << "Exact energy = " << exact_energy << "\n";
    std::cout << "\tError = " << std::abs(se.e() - exact_energy) <<
              " (Ry), relative error = " <<
              std::abs((se.e() - exact_energy)/exact_energy)*100 << " (%)\n";
	for(auto val : se.psi()){
		outfile << val << "\n";
	}
    outfile << "\n" << std::endl;

}

void coulomb_potential()
{
    unsigned int len = 50001;
    double rmax = 52;
    int z = 92;
    int n = 1;
    int l = 0;
    auto exact_energy = -GSL::pow_int(static_cast<double>(z)/n, 2);


    std::vector<double> v(len);
    Exponential_mesh<1, double> mesh(0, rmax, 1e-3, len);
    auto pot = [z](double r){return -2.*z/r;};
    auto mesh_it = mesh.begin();
    for(auto it = v.begin(); it != v.end(); it++, mesh_it++){
	    *it = pot(mesh_it->r());
    }
    std::vector<double> left;
    if (l == 0){
        double A = 1, C = 2*z;
        left = std::vector<double>{
            GSL::pow_int(-1, (n - l - 1))*(A - C/2*mesh.r(0)),
            GSL::pow_int(-1, (n - l - 1))*(A - C/2*mesh.r(1)),
            GSL::pow_int(-1, (n - l - 1))*(A - C/2*mesh.r(2))
        };
    } else {
        left = std::vector<double>{
            GSL::pow_int(-1, (n - l - 1))*GSL::pow_int(mesh.r(0), l),
            GSL::pow_int(-1, (n - l - 1))*GSL::pow_int(mesh.r(1), l),
            GSL::pow_int(-1, (n - l - 1))*GSL::pow_int(mesh.r(2), l)
        };
    }

    std::vector<double> right = {1e-15, 0};

    unsigned int nodes = static_cast<unsigned int>(std::max(n - l - 1, 0));
    Radial_Schroedinger_Equation_Central_Potential se(v, static_cast<unsigned int>(l), left, right, mesh, 1e-15);
    se.solve(nodes, -0.1, true);
    // se.normalize();

    std::cout << "Found energy : " << se.e() << "\n";
    std::cout << "Exact energy = " << exact_energy << "\n";
    std::cout << "\tError = " << std::abs(se.e() - exact_energy) <<
              " (Ry), relative error = " <<
              std::abs((se.e() - exact_energy)/exact_energy)*100 << " (%)\n";
    auto psi_i = se.psi().begin();
    auto v_i = se.v().begin();
    mesh_it = mesh.begin();
    while(psi_i != se.psi().end()){
	    outfile << mesh_it->r() << " " << *psi_i << " " << *v_i << "\n";
	    psi_i++;
	    v_i++;
	    mesh_it++;
    }
    outfile << "\n\n"<< std::flush;
}

int main()
{
    outfile.open("schroedinger_test.dat", std::ios::out);
    std::cout << "Harmonic Oscillator" << "\n";
    std::cout << std::string(80,'=') << "\n";
    harmonic_oscillator();
    std::cout << "\n";
    std::cout << "Coloumb potential" << "\n";
    std::cout << std::string(80,'=') << "\n";
    coulomb_potential();
    std::cout << "\n";
    outfile.close();
	return 0;
}
