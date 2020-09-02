#include <gtest/gtest.h>
#include "GSLpp/error.h"

int main(int argc, char **argv)
{
	GSL::Error_handler e_handler;
	e_handler.off();
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
