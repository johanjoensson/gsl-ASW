#include <log_mesh.h>
//#include "numerov_solver.h"
#include <GSLpp/special_functions.h>
#include <iostream>
#include <fstream>

int main()
{
    std::fstream out_file("mesh_test.dat", std::ios::out);
    Logarithmic_mesh mesh(2*M_PI, 501);
    auto sinc = [](double x){return GSL::sinc(10*x).val;};
    std::vector<double> f(mesh.size());
    for(size_t i = 0; i < mesh.size(); i++){
        f[i] = sinc(mesh.r(i));
        out_file << mesh.r(i) << " " << mesh.r(i+1) - mesh.r(i) << " " << mesh.drx(i) << " " << f[i] << "\n";
    }
    out_file.close();
    std::cout << mesh.integrate(f) << "\n";
    std::cout << mesh.integrate_simpson(f) << "\n";
    return 0;
}
