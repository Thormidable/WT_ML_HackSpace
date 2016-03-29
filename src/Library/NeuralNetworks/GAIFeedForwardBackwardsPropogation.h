
#pragma once

namespace GAI
{

	template<class T, class NeuronFunction = NeuronFunctions::Sigmoid<T> > class FeedForwardDense
	{
	public:
		FeedForwardDense(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims);
		void Resize(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims);
		void Randomise();

		void SetRegularisation(T lReg){ mRegularisationFactor = lReg; };
		T GetRegularisation(){ return mRegularisationFactor; }

		void Train(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);

		MatrixDynamic<T> Process(const MatrixDynamic<T> &lInputValues);
		MatrixDynamic<T> ProcessFull(const MatrixDynamic<T> &lInputValues);
		inline MatrixDynamic<T> GetResultMatrix(){ return mOutputValues; }

		T CalculateCostValues(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);

		Bool TestCostFunction(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);

		void MoveWeightsToArray(lbfgsfloatval_t *);
		void MoveArrayToWeights(const lbfgsfloatval_t *);

		std::vector<MatrixDynamic<T>> CalculateCostGradient(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);
		T CalculateCostValuesFromResults(const MatrixDynamic<T> &lExpected);

	protected:		
		
		MatrixDynamic<T> mInputWeights;
		std::vector<MatrixDynamic<T>> mLayerWeights;		
		std::vector<MatrixDynamic<T>> mA;
		std::vector<MatrixDynamic<T>> mZ;
		MatrixDynamic<T> mOutputValues;

		Size mInputDims;
		Size mLayers;
		Size mLayerSize;

		T mRegularisationFactor;
	};

}
