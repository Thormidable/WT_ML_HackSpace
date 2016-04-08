
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
		T GetRegularisation()const{ return mRegularisationFactor; }

		void Train(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);
		void StochasticTraining(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected, Size lBatchSize, T lfScale = T(1.0), T lfTolerance = T(0.00001));
		void TrainBatch(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected, T lfScale);

		MatrixDynamic<T> Process(const MatrixDynamic<T> &lInputValues);
		MatrixDynamic<T> ProcessFull(const MatrixDynamic<T> &lInputValues);
		inline MatrixDynamic<T> GetResultMatrix()const{ return mOutputValues; }

		T CalculateCostValues(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);

		Bool TestCostFunction(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);

		void MoveWeightsToArray(lbfgsfloatval_t *);
		void MoveArrayToWeights(const lbfgsfloatval_t *);
		void ReweightGradients(std::vector<MatrixDynamic<T>> &lGradients);

		std::vector<MatrixDynamic<T>> CalculateCostGradient(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);
		T CalculateCostValuesFromResults(const MatrixDynamic<T> &lExpected)const;

		Size GetInputDims()const{ return mLayerWeights.size() > 0 ? mLayerWeights.front().rows() : 0; }
		Size GetOutputDims()const{ return mLayerWeights.size() > 0 ? mLayerWeights.back().cols() : 0; }
		Size GetLayerSize()const{ return mLayerWeights.size() > 0 ? mLayerWeights.front().cols() : 0; }
		Size GetNumLayers()const{ return mLayerWeights.size() > 0 ? mLayerWeights.size() - 1 : 0; }

	protected:		
		
		MatrixDynamic<T> mInputWeights;
		std::vector<MatrixDynamic<T>> mLayerWeights;		
		std::vector<MatrixDynamic<T>> mA;
		std::vector<MatrixDynamic<T>> mZ;
		MatrixDynamic<T> mOutputValues;

		
		T mRegularisationFactor;
	};

}
