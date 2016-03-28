
#include "GAILibrary.h"

using namespace GAI;

template<class T, class NeuronFunction> GAI::FeedForwardDense<T,NeuronFunction>::FeedForwardDense(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims)
{
	Resize(lInputDims, lLayerSize, lLayers, lOutputDims);
}

template<class T, class NeuronFunction> void GAI::FeedForwardDense<T,NeuronFunction>::Resize(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims)
{
	if (lLayerSize == 0) printf("ERROR. Numbers of Layers too low. for Forward Dense");

	if (mLayerSize != lLayerSize || lLayers != mLayers || (mLayerWeights.size() == 0 || mLayerWeights[0].rows() != lInputDims) || mLayerWeights.back().cols() != lOutputDims)
	{
		mLayerSize = lLayerSize;
		mLayers = lLayers;
		mLayerWeights.resize(lLayers+1);
		mLayerWeights[0].resize(lInputDims, lLayerSize);

		for (Size lLayer = 1; lLayer < lLayers; ++lLayer){mLayerWeights[lLayer].resize(lLayerSize, lLayerSize);}

		mLayerWeights.back().resize(lLayerSize, lOutputDims);		
	}

	Randomise();

}

template<class T, class NeuronFunction> void GAI::FeedForwardDense<T,NeuronFunction>::Randomise()
{
	if (mLayerWeights.size() == 0) throw;

	for (auto &i : mLayerWeights)
	{		
		//MatrixIterator(i, [](T &lData){lData = Random<T>(0, 1.0); });
		MatrixIterator(i, [](T &lData){lData = 0.5f; });
	}
}

template<class T, class NeuronFunction> void GAI::FeedForwardDense<T, NeuronFunction>::Train(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected)
{
	auto lDJDW = CalculateCostGradient(lInputValues, lExpected);	

	T lfScale = 0.1f;
	for (Size lLayer = 0; lLayer < mLayerWeights.size(); ++lLayer)
	{
		mLayerWeights[lLayer] -= lDJDW[lLayer] * lfScale;
	}
}

template<class T, class NeuronFunction> T GAI::FeedForwardDense<T, NeuronFunction>::CalculateCostValues(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected)
{
	auto lResults = ProcessFull(lInputValues);
	auto lEval = (lResults - lExpected);	
	return (lEval.cwiseProduct(lEval).array()).sum()*T(0.5);
}

template<class T, class NeuronFunction> std::vector<MatrixDynamic<T>> GAI::FeedForwardDense<T, NeuronFunction>::CalculateCostGradient(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected)
{
	auto lYHat = ProcessFull(lInputValues);

	MatrixDynamic<T> lD3;
	std::vector<MatrixDynamic<T>> lDJdW;
	lDJdW.resize(mLayerWeights.size());
	lDJdW.back() = -(lExpected-lYHat);

	for (Size lLayer = lDJdW.size()-1; lLayer > 0; --lLayer)
	{
		//Calculate sigmoidPrime of Input Energys
		auto lFPrimeZ3 = mZ[lLayer];
//		auto lfTestDiff = (lFPrimeZ3 - (mA[lLayer - 1] * mLayerWeights[lLayer])).eval();
		MatrixIterator<T>(lFPrimeZ3, &NeuronFunction::Prime);

		//Calculate ErrorFactor * sigmoidPrime(Input energies)
		lD3 = lDJdW[lLayer].cwiseProduct(lFPrimeZ3);
		lDJdW[lLayer] = mA[lLayer - 1].transpose() * lD3;
		if(lLayer > 1) lDJdW[lLayer - 1] = (lD3*mLayerWeights[lLayer].transpose()).eval();
	}

	auto lFPrimeZ3 = (mZ[0]).eval();
	MatrixIterator<T>(lFPrimeZ3, &NeuronFunction::Prime);

	lD3 = (lD3*mLayerWeights[1].transpose())*lFPrimeZ3;
	lDJdW[0] = lInputValues.transpose()*lD3;

	return lDJdW;
}

template<class T, class NeuronFunction> Bool GAI::FeedForwardDense<T, NeuronFunction>::TestCostFunction(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected)
{
	auto lfCalculatedGrads = CalculateCostGradient(lInputValues, lExpected);

	MatrixDynamic<T> lPerturbation;
	lPerturbation.resize(1, 3);

	T lEps = 0.0001f;

	for (Size lRow = 0; lRow < mLayerWeights.back().rows(); ++lRow)
	{
		for (Size lCol = 0; lCol < mLayerWeights.back().cols(); ++lCol)
		{
			auto lWeighted = mLayerWeights;
			mLayerWeights.back()(lRow, lCol) += lEps;
			T lLoss1 = CalculateCostValues(lInputValues, lExpected);

			mLayerWeights = lWeighted;
			mLayerWeights.back()(lRow, lCol) -= lEps;
			T lLoss2 = CalculateCostValues(lInputValues, lExpected);

			mLayerWeights = lWeighted;

			auto lfGrads = (lLoss1 - lLoss2) / (2 * lEps);

			auto lDiffResult = (lfCalculatedGrads.back()(lRow, lCol) - lfGrads);
			if (abs(lDiffResult) > 0.0000001f) return false;
		}
	}
	return true;
}

template<class T, class NeuronFunction> MatrixDynamic<T> GAI::FeedForwardDense<T,NeuronFunction>::Process(const MatrixDynamic<T> &lInputValues)
{
	MatrixDynamic<T> lResult = lInputValues;
	lResult = ProcessLayer(lResult, mLayerWeights[0]);

	for (Size lLayer = 1; lLayer < mLayerWeights.size(); ++lLayer)
	{
		lResult = ProcessLayer(lResult, mLayerWeights[lLayer]);
	}

	return lResult;
}

template<class T, class NeuronFunction> MatrixDynamic<T> GAI::FeedForwardDense<T, NeuronFunction>::ProcessFull(const MatrixDynamic<T> &lInputValues)
{
	mA.resize(mLayerWeights.size() - 1);
	mZ.resize(mLayerWeights.size());

	mA[0] = lInputValues;
	for (Size lLayer = 0; lLayer < mA.size(); ++lLayer)
	{
		mZ[lLayer] = mA[lLayer] * mLayerWeights[lLayer];
		mA[lLayer] = ProcessLayer(mA[lLayer], mLayerWeights[lLayer]);
		if (lLayer < mA.size() - 1) mA[lLayer + 1] = mA[lLayer];
	}
	mZ.back() = mA.back()*mLayerWeights.back();
	return ProcessLayer(mA.back(), mLayerWeights.back());
}

template<class T, class NeuronFunction> MatrixDynamic<T> GAI::FeedForwardDense<T, NeuronFunction>::ProcessLayer(const MatrixDynamic<T> &lInputs, const MatrixDynamic<T> &lWeights)
{		
	auto lResult = (lInputs*lWeights).eval();
	MatrixIterator(lResult, &NeuronFunction::Function);
	return lResult;
}

#define INSTANCE_NEURON_FUNCTIONS(T) \
template class GAI::FeedForwardDense<T,NeuronFunctions::Sigmoid<T> >;

GAI_EXPAND(GAI_FOR_EACH_FLOAT_TYPE(INSTANCE_NEURON_FUNCTIONS))
