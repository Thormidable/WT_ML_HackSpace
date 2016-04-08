
#pragma once

namespace GAI
{
/*
	template<class T, class NeuronFunction = NeuronFunctions::Sigmoid<T> > class FeedForwardEnergised : public FeedForwardDense<T,NeuronFunction>
	{
	public:
		FeedForwardEnergised(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims);		

		void Randomise();
		void Train(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);

		MatrixDynamic<T> Process(const MatrixDynamic<T> &lInputValues);
		MatrixDynamic<T> ProcessFull(const MatrixDynamic<T> &lInputValues);
		
		std::vector<MatrixDynamic<T>> CalculateCostGradient(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);
		
	protected:
		std::vector<MatrixDynamic<T>> mEnergiserWeights;
		std::vector<MatrixDynamic<T>> mEnergiserGradients;
	};
*/
}
