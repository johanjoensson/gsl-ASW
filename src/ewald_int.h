#ifndef EWALD_INT_H
#define EWALD_INT_H

#include <vector>
#include "log_mesh.h"

class Integral_hankel{
        Logarithmic_mesh &mesh;
        double ewald_param;
	double kappa;
    public:
        std::vector<double> h;

        void set_ewald_param(double eta);
	void set_kappa(double kappa);

        Integral_hankel(Logarithmic_mesh &mesh);
};

#endif //EWALD_INT_H
