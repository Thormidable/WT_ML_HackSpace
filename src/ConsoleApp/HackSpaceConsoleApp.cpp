
#include "HackSpaceConsoleApp.h"

using namespace GAI;

int main()
{

	time_t lTime;
	time(&lTime);
	RandomSeed(UInt32(lTime));

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

	MatrixDynamic<Float64> lTestingInput;
	lTestingInput.resize(4, 2);
	lTestingInput << 4, 4, 9, 6,
					5.5, 1, 2.5, 2;
	MatrixDynamic<Float64> lTestingExpected;
	lTestingExpected.resize(4, 1);
	lTestingExpected << 0.7, 0.89, 0.85, 0.75;

	if (!lFFDense.TestCostFunction(lInput, lExpected)) return 1;

	lFFDense.Train(lInput, lExpected);

	DisplayMatrix(lFFDense.GetResultMatrix(), "Trained Results Matrix");

	Float64 lfTrainedError = lFFDense.CalculateCostValues(lInput, lExpected);
	Float64 lfTestingError = lFFDense.CalculateCostValues(lTestingInput, lTestingExpected);

	if (lfTestingError> lfTrainedError*Float64(10))
	{
		printf("Testing does not match Training results : %f is much greater than %f\n", lfTestingError, lfTrainedError);
	}
	

	system("pause");
	return 0;
}

