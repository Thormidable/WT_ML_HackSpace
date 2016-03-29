
#include "HackSpaceConsoleApp.h"

using namespace GAI;

int main()
{

	time_t lTime;
	time(&lTime);
	RandomSeed(lTime);

	FeedForwardDense<Float64, NeuronFunctions::Sigmoid<Float64>>  lFFDense(2, 3, 1, 1);

	MatrixDynamic<Float64> lInput;
	//lInput.resize(3, 2);
	//lInput << 0.3, 1.0, 0.5, 0.2, 1.0, 0.4;
	//lInput << 0.3, 0.5, 1.0, 1.0, 0.2, 0.4;
	//lInput << 0.3, 0.3, 0.3, 0.3, 0.3, 0.3;
	//lInput << 0.3, 0.4, 0.5, 0.75, 0.82, 0.93;

	lInput.resize(4, 2);
	lInput << 3, 5, 10, 6,
		5, 1, 2, 1.5;

	MatrixDynamic<Float64> lExpected;
	/*lExpected.resize(3, 1);
	lExpected << 0.75, 0.82, 0.93;*/
	lExpected.resize(4, 1);
	lExpected << 0.75, 0.82, 0.93, 0.7; 

	DisplayMatrix(lInput, "Input Matrix");
	DisplayMatrix(lExpected, "Expected Results Matrix");

	if (!lFFDense.TestCostFunction(lInput, lExpected)) return 1;

	lFFDense.Train(lInput, lExpected);

	system("pause");
	return 0;
}

