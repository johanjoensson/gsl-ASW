#ifndef NUMEROV_SOLVER_H
#define NUMEROV_SOLVER_H

#include <stdexcept>
#include <iostream>

template<class T>
int signum(T val)
{
	return (val > T(0)) - (val < T(0));
}

/***************************************************************************//**
* A class for solving differential equations using the Numerov method
*******************************************************************************/
class Numerov_solver{
	/***********************************************************************//**
	* Find the point where g(x_i) = 0 (classical inversion point).\n
	* \t__Iter_res__ - iterator to result container (e.g. std::vector)\n
	* \t__Iter_g__ - iterator to g container (e.g. std::vector)\n
	* Input:\n
	* \t__res_start__ - Iterator pointing to the first element of the result container.\n
	* \t__res_end__ - Iterator pointing to the element one past the last element of the result container.\n
	*\t__g_start__ - Iterator pointing to the first element of  the g container
	***************************************************************************/
	template<class Iter_res, class Iter_g>
	Iter_g find_inversion_point(Iter_res res_start, Iter_res res_end,  Iter_g
		g_start)
	{
		int init_sign = signum(*g_start);
		Iter_res current = res_start;
		Iter_g current_g = g_start;
		while(current != res_end && signum(*current_g) == init_sign){
			current++;
			current_g++;
		}
		return current;
	}

	template<class Iter_res, class Iter_g, class Iter_s, class Iter_init>
	Iter_res solve_direction(Iter_res res_start, Iter_res res_end, Iter_init
		init_start, Iter_init init_end, Iter_g g_start, Iter_s s_start)
	{
		Iter_res current = res_start, previous = res_end, pprevious = res_end;
		Iter_g current_g = g_start, previous_g, pprevious_g;
		Iter_s current_s = s_start, previous_s, pprevious_s;

		// Set initial conditions
		for(Iter_init& init = init_start; init != init_end && current !=
				res_end; init++){
			*current = *init;
			pprevious = previous;
			previous = current;
			current++;

			pprevious_g = previous_g;
			previous_g = current_g;
			current_g++;

			pprevious_s = previous_s;
			previous_s = current_s;
			current_s++;
		}

		// Solve equation
		while(current != res_end){
			*current = 	2*(*previous)*(1 - 5./12*(*previous_g)) -
						(*pprevious)*(1 + 1./12*(*pprevious_g)) +
						1./12*((*current_s) + 10*(*previous_s) +(*pprevious_s));
			*current /= 1 + 1./12*(*current_g);

			pprevious = previous;
			pprevious_g = previous_g;
			pprevious_s = previous_s;
			previous = current;
			previous_g = current_g;
			previous_s = current_s;
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
			Iter_right_init right_start, Iter_right_init right_end)
		{
			Iter_res range_start = res_start, range_end = res_end;
			// Find matching point, inv points to the matching point
			// or res_end if no matching point was found
			Iter_res inv = find_inversion_point(range_start, range_end,
				g_start);
			// If we found an inversion point, we want to include it in our
			// numerical solution
			if(inv != res_end){
				inv++;
			}
			std::reverse_iterator<Iter_res> rinv;
			// Solve | -> .
			inv = solve_direction(res_start, inv, left_start, left_end,
				g_start, s_start);
			// inv now points to the element after the matching point
			// Extract value at matching point
			T scale = *(--inv);
			// Solve . <- |
			rinv = solve_direction(std::reverse_iterator<Iter_res>(res_end--),
				std::reverse_iterator<Iter_res>(inv)++, right_start, right_end,
				std::reverse_iterator<Iter_g>(g_end--),
				std::reverse_iterator<Iter_s>(s_end--));
			// Match values at matching point
			scale /= *(--rinv);
			for(auto& iter = inv; iter != res_end; iter++){
			 	*iter *= scale;
			}
			return inv;
		}
};

#endif //NUMEROV_SOLVER_H
