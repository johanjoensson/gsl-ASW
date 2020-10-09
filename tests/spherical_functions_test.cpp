#include <gtest/gtest.h>
#include "../src/spherical_fun.h"
#include "GSLpp/vector.h"

#define TOL 1e-12


TEST(CubicHarmonic, Y00)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5), static_cast<double>(j - 5), static_cast<double>(k - 5)});
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({0,0}, v).val, 1./(std::sqrt(4*M_PI)), TOL) << v;
    }

}

TEST(CubicHarmonic, Y1m1)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({1, -1}, v).val, std::sqrt(3./(4*M_PI))*v[1], TOL) << v;
    }
}

TEST(CubicHarmonic, Y10)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({1, 0}, v).val, std::sqrt(3./(4*M_PI))*v[2], TOL) << v;
    }
}

TEST(CubicHarmonic, Y11)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({1, 1}, v).val, std::sqrt(3./(4*M_PI))*v[0], TOL) << v;
    }
}

TEST(CubicHarmonic, Y2m2)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({2, -2}, v).val, std::sqrt(15./(4*M_PI))*v[0]*v[1], TOL) << v;
    }
}

TEST(CubicHarmonic, Y2m1)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({2, -1}, v).val, std::sqrt(15./(4*M_PI))*v[1]*v[2], TOL) << v;
    }
}

TEST(CubicHarmonic, Y20)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({2, 0}, v).val, std::sqrt(5./(16*M_PI))*(3*GSL::pow_int(v[2], 2) - GSL::pow_int(v.norm<double>(), 2)), TOL) << v;
    }
}

TEST(CubicHarmonic, Y21)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({2, 1}, v).val, std::sqrt(15./(4*M_PI))*v[0]*v[2], TOL) << v;
    }
}

TEST(CubicHarmonic, Y22)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({2, 2}, v).val, std::sqrt(15./(16*M_PI))*(GSL::pow_int(v[0], 2) - GSL::pow_int(v[1], 2)), TOL) << v;
    }
}

TEST(CubicHarmonic, Y3m3)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
	double x = v[0], y = v[1];
	double x2 = x*x, y3 = y*y*y;
        ASSERT_NEAR(cubic_harmonic({3, -3}, v).val, std::sqrt(35/(32*M_PI))*(3*x2*y - y3), TOL) << v;
    }
}

TEST(CubicHarmonic, Y3m2)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
	double x = v[0], y = v[1], z = v[2];
        ASSERT_NEAR(cubic_harmonic({3, -2}, v).val, std::sqrt(105/(4*M_PI))*x*y*z, TOL) << v;
    }
}

TEST(CubicHarmonic, Y3m1)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
	double x = v[0], y = v[1], z = v[2], r2 = v.norm2();
	double z2 = z*z;
        ASSERT_NEAR(cubic_harmonic({3, -1}, v).val, std::sqrt(21/(32*M_PI))*(5*y*z2 - y*r2), TOL) << v;
    }
}

TEST(CubicHarmonic, Y30)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
	double z = v[2], r2 = v.norm2();
	double z3 = z*z*z;
        ASSERT_NEAR(cubic_harmonic({3, 0}, v).val, std::sqrt(7/(16*M_PI))*(5*z3 - 3*z*r2), TOL) << v;
    }
}

TEST(CubicHarmonic, Y31)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
	double x = v[0], y = v[1], z = v[2], r2 = v.norm2();
	double z2 = z*z;
        ASSERT_NEAR(cubic_harmonic({3, 1}, v).val, std::sqrt(21/(32*M_PI))*(5*x*z2 - x*r2), TOL) << v;
    }
}

TEST(CubicHarmonic, Y32)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
	double x = v[0], y = v[1], z = v[2];
	double x2 = x*x, y2 = y*y;
        ASSERT_NEAR(cubic_harmonic({3, 2}, v).val, std::sqrt(105/(16*M_PI))*(x2 - y2)*z, TOL) << v;
    }
}

TEST(CubicHarmonic, Y33)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({static_cast<double>(i - 5) , static_cast<double>(j - 5) , static_cast<double>(k - 5) });
            }
        }
    }

    for(auto v : vecs){
	double x = v[0], y = v[1], z = v[2];
	double x3 = x*x*x, y2 = y*y;
        ASSERT_NEAR(cubic_harmonic({3, 3}, v).val, std::sqrt(35/(32*M_PI))*(x3 - 3*x*y2), TOL) << v;
    }
}
