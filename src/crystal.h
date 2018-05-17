#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <vector>
#include "atom.h"

class Crystal {
	std::vector<Atom> atoms;
	public:
		Crystal();
};
#endif //CRYSTAL_H
