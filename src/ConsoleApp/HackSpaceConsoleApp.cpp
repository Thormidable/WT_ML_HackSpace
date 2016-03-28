
#include "HackSpaceConsoleApp.h"

using namespace GAI;

int main()
{
	FeedForwardDense<Float64>  lFFDense(2, 3, 3, 2);

	MatrixDynamic<Float64> lInput;
	lInput.resize(3, 2);
	lInput << 3, 5, 5, 1, 10, 2;

	auto lResult = lFFDense.Process(lInput);

	return 0;
}

