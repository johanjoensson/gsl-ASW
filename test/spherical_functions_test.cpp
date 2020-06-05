#include <gtest/gtest.h>
#include "../src/spherical_fun.h"
#include "GSLpp/vector.h"


TEST(CubicHarmonic, Y00)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({i - 5, j - 5, k - 5});
            }
        }
    }

    for(auto v : vecs){
        ASSERT_DOUBLE_EQ(cubic_harmonic({0,0}, v).val, 1./(std::sqrt(4*M_PI)));
    }

}

TEST(CubicHarmonic, Y1m1)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({i - 5 , j - 5 , k - 5 });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({1, -1}, v).val, std::sqrt(3./(4*M_PI))*v[1], 5e-12);
    }
}

TEST(CubicHarmonic, Y10)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({i - 5 , j - 5 , k - 5 });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({1, 0}, v).val, std::sqrt(3./(4*M_PI))*v[2], 5e-12);
    }
}

TEST(CubicHarmonic, Y11)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({i - 5 , j - 5 , k - 5 });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({1, 1}, v).val, std::sqrt(3./(4*M_PI))*v[0], 5e-12);
    }
}

TEST(CubicHarmonic, Y2m2)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({i - 5 , j - 5 , k - 5 });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({2, -2}, v).val, std::sqrt(15./(4*M_PI))*v[0]*v[1], 5e-12);
    }
}

TEST(CubicHarmonic, Y2m1)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({i - 5 , j - 5 , k - 5 });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({2, -1}, v).val, std::sqrt(15./(4*M_PI))*v[1]*v[2], 5e-12);
    }
}

TEST(CubicHarmonic, Y20)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({i - 5 , j - 5 , k - 5 });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({2, 0}, v).val, std::sqrt(5./(16*M_PI))*(3*GSL::pow_int(v[2], 2) - GSL::pow_int(v.norm<double>(), 2)), 5e-12);
    }
}

TEST(CubicHarmonic, Y21)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({i - 5 , j - 5 , k - 5 });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({2, 1}, v).val, std::sqrt(15./(4*M_PI))*v[0]*v[2], 5e-12);
    }
}

TEST(CubicHarmonic, Y22)
{
    std::array<GSL::Vector, 1000> vecs;
    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            for(int k = 0; k < 10; k++){
                vecs[100*i + 10*j + k] = GSL::Vector({i - 5 , j - 5 , k - 5 });
            }
        }
    }

    for(auto v : vecs){
        ASSERT_NEAR(cubic_harmonic({2, 2}, v).val, std::sqrt(15./(16*M_PI))*(GSL::pow_int(v[0], 2) - GSL::pow_int(v[1], 2)), 5e-12);
    }
}
