
#include "HackSpaceConsoleApp.h"

using namespace GAI;


template<class T, class TData> void DisplayInputOutput(const T &lFFDense, const MatrixDynamic<TData> &lInput, const MatrixDynamic<TData> &lExpected)
{
	DisplayMatrix(lInput, "Input Matrix");
	DisplayMatrix(lExpected, "Expected Results Matrix");
	DisplayMatrix(lFFDense.GetResultMatrix(), "Trained Results Matrix");
}

int main()
{

	time_t lTime;
	time(&lTime);
	RandomSeed(UInt32(lTime));

	//FeedForwardDense<Float64, NeuronFunctions::Sigmoid<Float64>>  lFFDense(2, 3, 1, 1);

	//// Smooth data
	//MatrixDynamic<Float64> lInput;
	//lInput.resize(3, 2);	
	//lInput <<3, 5, 10, 5, 1, 2;

	//MatrixDynamic<Float64> lExpected;
	//lExpected.resize(3, 1);
	//lExpected << 0.75, 0.82, 0.93;

	//if (!lFFDense.TestCostFunction(lInput, lExpected)) return 1;

	//DisplayMatrix(lInput, "Input Matrix");
	//DisplayMatrix(lExpected, "Expected Results Matrix");
	//lFFDense.Randomise();
	//lFFDense.Train(lInput, lExpected);
	//DisplayMatrix(lFFDense.GetResultMatrix(), "Trained Results Matrix");

	//// Heavy Noisy Data
	//lInput.resize(4, 2);
	//lInput << 3, 5, 10, 6,
	//	5, 1, 2, 1.5;

	//lExpected.resize(4, 1);
	//lExpected << 0.75, 0.82, 0.93, 0.7; 
	//
	//if (!lFFDense.TestCostFunction(lInput, lExpected)) return 1;
	//DisplayMatrix(lInput, "Input Matrix");
	//DisplayMatrix(lExpected, "Expected Results Matrix");
	//lFFDense.Randomise();
	//lFFDense.Train(lInput, lExpected);
	//DisplayMatrix(lFFDense.GetResultMatrix(), "Trained Results Matrix");

	////Different Data Set to test validity of fit
	//MatrixDynamic<Float64> lTestingInput;
	//lTestingInput.resize(4, 2);
	//lTestingInput << 4, 4, 9, 6,
	//				5.5, 1, 2.5, 2;
	//MatrixDynamic<Float64> lTestingExpected;
	//lTestingExpected.resize(4, 1);
	//lTestingExpected << 0.7, 0.89, 0.85, 0.75;

	//Float64 lfTrainedError = lFFDense.CalculateCostValues(lInput, lExpected);
	//Float64 lfTestingError = lFFDense.CalculateCostValues(lTestingInput, lTestingExpected);

	//if (lfTestingError > lfTrainedError*Float64(10))
	//{
	//	printf("Testing does not match Training results : %f is much greater than %f\n", lfTestingError, lfTrainedError);
	//}

	//StochasticTraining Test
	FeedForwardDense<Float64, NeuronFunctions::Sigmoid<Float64>>  lFFDense2(10, 15, 1, 1);

	MatrixDynamic<Float64> lSchotasticTraining;
	//lSchotasticTraining.resize(8, 2);
	//lSchotasticTraining << 4, 4, 9, 6, 3, 5, 10, 6,		
	//	5.5, 1, 2.5, 2, 5, 1, 2, 1.5;

	/*lSchotasticTraining.resize(4, 2);
	lSchotasticTraining << 3, 5, 10, 6, 5, 1, 2, 1.5;*/
	lSchotasticTraining = TestData::FullBreastCancerInputs();

	MatrixDynamic<Float64> lSchotasticExpected;
	//lSchotasticExpected.resize(8, 1);
	//lSchotasticExpected << 0.7, 0.89, 0.85, 0.75, 0.75, 0.82, 0.93, 0.7;
	/*lSchotasticExpected.resize(4, 1);
	lSchotasticExpected << 0.75, 0.82, 0.93, 0.7;*/
	lSchotasticExpected = TestData::FullBreastCancerResults();	

	lSchotasticTraining.conservativeResize(lSchotasticTraining.rows(),lFFDense2.GetInputDims());
	lSchotasticExpected *= 0.5;
	lSchotasticExpected  = lSchotasticExpected.array() + 0.25;
	NormaliseInputDataSet(lSchotasticTraining);

	if (!lFFDense2.TestCostFunction(lSchotasticTraining, lSchotasticExpected)) return 1;
	lFFDense2.Train(lSchotasticTraining, lSchotasticExpected);
	lFFDense2.Save(L"./MalignancyTest.NN");

	//lFFDense2.Randomise();
	//lFFDense2.StochasticTraining(lSchotasticTraining, lSchotasticExpected,10,1,0.0000001);
	
	system("pause");
	return 0;
}
