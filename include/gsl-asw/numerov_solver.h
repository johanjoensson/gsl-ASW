#ifndef NUMEROV_SOLVER_H
#define NUMEROV_SOLVER_H

#include <stdexcept>
#include <iostream>

/***************************************************************************//**
* A class for solving differential equations using the Numerov method
*******************************************************************************/
class Numerov_solver{
	template<class T>
	static int signum(T val)
	 noexcept
	{
		return (val > T(0)) - (val < T(0));
	}

	/***********************************************************************//**
	* Find the point where g(x_i) = 0 (classical inversion point).\n
	* \t__Iter_res__ - iterator to result container (e.g. std::vector)\n
	* \t__Iter_g__ - iterator to g container (e.g. std::vector)\n
	* Input:\n
	* \t__res_start__ - Iterator pointing to the first element of the result container.\n
	* \t__res_end__ - Iterator pointing to the element one past the last element of the result container.\n
	*\t__g_start__ - Iterator pointing to the first element of  the g container
	***************************************************************************/

	template<class Iter_res, class Iter_g, class Iter_s, class Iter_init>
	Iter_res solve_direction(Iter_res res_start, Iter_res res_end, Iter_init
		init_start, Iter_init init_end, Iter_g g_start, Iter_s s_start)
	 const noexcept
	{
		Iter_res current = res_start;
		Iter_g current_g = g_start;
		Iter_s current_s = s_start;

		// Set initial conditions
		for(Iter_init& init = init_start; init != init_end && current != res_end;
			init++){
			*current = *init;
			current++;

			current_g++;

			current_s++;
		}

		auto f = [](Iter_g g){return 1 + *g/12;};
		// Solve equation
		while(current != res_end){
			*current = 	*(current - 1)*(12 - 10*f(current_g - 1)) -
						*(current - 2)*f(current_g - 2) +
						1./12*(*current_s + 10*(*(current_s -1)) + *(current_s - 2));

			*current /= f(current_g);
			current++;
			current_g++;
			current_s++;
		}
		return current;
	}
	public:
		Numerov_solver() = default;

	/***********************************************************************//**
	* Use Numerov's algorithm to solve the differential equation.
	* \t__Iter_res__ - iterator to result container (e.g. std::vector)\n
	* \t__Iter_g__ - iterator to g container (e.g. std::vector)\n
	* Input:\n
	* \t__res_start__ - iterator to the first element in the result container\n
	* \t__res_end__ - iterator one step part the last element in the result container\n
	* \t__g_start__ - iterator to the start of the g container\n
	* \t__g_end__ - iterator to one step past the last element of the g container\n
	* \t__s_start__ - iterator to the start of the s container\n
	* \t__s_end__ - iterator to one step past the last element of the s container\n
	* \t__left_start__ - iterator to the start of the container containing the left hand initial conditions\n
	* \t__left_end__ - iterator to one step past the last element of the left boundary conditions container\n
	* \t__right_start__ - iterator to the start of the container containing the right hand initial conditions\n
	* \t__right_end__ - iterator to one step past the last element of the right boundary conditions container\n
	Returns:\n
	\t Iterator to the matching point.\n
	NOTE : The elements of g must be multiplied with the step size squared (h^2).
	***************************************************************************/
	template<class Iter_res, class Iter_g, class Iter_s, class Iter_left_init , class Iter_right_init, class T = double>
	Iter_res solve(Iter_res res_start, Iter_res res_end,
			Iter_g g_start, Iter_g g_end, Iter_s s_start, Iter_s s_end,
			Iter_left_init left_start, Iter_left_init left_end,
			Iter_right_init right_start, Iter_right_init right_end, Iter_res inv)
	 const noexcept
	{
		// If no initial conditions at the right hand side, only solve to the right
		if(right_start == right_end){
			inv = res_end;
		}
		// Solve | -> .
		inv = solve_direction(res_start, inv, left_start, left_end,
			g_start, s_start);
		// If we found an inversion point, solve in the other direction as well
		if(inv != res_end){
			// inv now points to the element after the matching point
			// Extract value at matching point
			T scale = (*(inv-1));
			// Solve . <- |
			auto rinv = solve_direction(std::reverse_iterator<Iter_res>(res_end),
				std::reverse_iterator<Iter_res>(inv - 1), right_start, right_end,
				std::reverse_iterator<Iter_g>(g_end),
				std::reverse_iterator<Iter_s>(s_end));
			// Match values at matching point
			scale /= *(rinv - 1);
			for(auto iter = inv - 1; iter != res_end; iter++){
				*iter *= scale;
			}

		}
		return inv;
	}


	/***********************************************************************//**
	* Solve in one direction only.
	* Overloaded version of solve.
	**************************************************************************/
	template<class Iter_res, class Iter_g, class Iter_s, class Iter_left_init , class Iter_right_init, class T = double>
	Iter_res solve(Iter_res res_start, Iter_res res_end,
			Iter_g g_start, Iter_g g_end, Iter_s s_start, Iter_s s_end,
			Iter_left_init left_start, Iter_left_init left_end)
	 const noexcept
	{
		solve(res_start, res_end, g_start, g_end, s_start, s_end, left_start,
			left_end, res_end, res_end);
	}

	/***********************************************************************//**
	* Use Numerov's algorithm to estimate the difference in derivatives at the matching point.
	* \t__Iter_res__ - iterator to result container (e.g. std::vector)\n
	* \t__Iter_g__ - iterator to g container (e.g. std::vector)\n
	* \t__Iter_s__ - iterator to s container (e.g. std::vector)\n
	* \t__T__ - Class representing numerical values (e.g. double)\n
	* Input:\n
	* \t__res_start__ - iterator to the first element in the result container\n
	* \t__res_end__ - iterator one step past the last element in the result container\n
	* \t__inv__ - iterator to result at inversion point\n
	* \t__g_start__ - iterator to the start of the g container\n
	* \t__s_start__ - iterator to the start of the s container\n
	Returns:\n
	\t Approximation of sign(y'_R(x_inv) - y'_L(x_inv)), sign of the difference in derivatives.\n
	NOTE : The elements of g must be premultiplied with the step size squared (h^2).
	***************************************************************************/
	template<class Iter_res, class Iter_g, class Iter_s, class T = double>
	int derivative_diff(Iter_res res_start, Iter_g res_end, Iter_res inv,
		Iter_g g_start, Iter_s s_start, double tol = 1e-6)
	 const noexcept
	{
		T val{0};
		if(inv == res_end){
			return signum(val);
		}
		// Iter_res res_c = res_start;
		Iter_g g_c = g_start;
		Iter_s s_c = s_start;
		// while(res_c != inv && res_c != res_end){
		for(auto res_c = res_start; res_c != inv; res_c++, g_c++, s_c++){}
		// val = *(res_c - 1) + *(res_c + 1) - (2 - *g_c)*(*res_c) - *s_c;
		val = *(inv - 2) + *inv - (2 - *(g_c - 1))*(*(inv - 1)) - *(s_c - 1);
		if(std::abs(val) < tol){
			val = 0;
		}
		return signum(val);
	}


	/***************************************************************************
	* Derivative expressions found using the same approach as                  *
	* V.I. Tselyaev,                                                           *
	* A generalized Numerov method for linear second-order differential        *
	* equations involving a first derivative term,                             *
	* Journal of Computational and Applied Mathematics,                        *
	* Volume 170, Issue 1,                                                     *
	* 2004,                                                                    *
	* https://doi.org/10.1016/j.cam.2003.12.042.                               *
	* (https://www.sciencedirect.com/science/article/pii/S0377042704000159)    *
	***************************************************************************/
	template<class Iter_res, class Iter_g, class Iter_s, class Iter_left_init , class Iter_right_init, class T = double>
	Iter_res derivative(Iter_res f_start, Iter_res f_end,
			Iter_g g_start, Iter_g g_end, Iter_s s_start, Iter_s s_end,
			Iter_left_init left_start, Iter_left_init left_end)
	 const noexcept
	{
		return f_start;
	}

};

#endif //NUMEROV_SOLVER_H
