
#pragma once

namespace GAI
{
	namespace NeuronFunctions
	{
		template<class T> class Sigmoid
		{
		public:
			static inline void Function(T &lValue)
			{
				lValue = 1 / (1 + std::exp(-lValue));
			}
			static inline void Prime(T &lValue)
			{
				lValue = exp(-lValue) / std::pow(1 + exp(-lValue),2);
			}
		};

		template<class T> class TanH
		{
		public:
			static void Function(T& lValue)
			{				
				lValue = std::tanh(lValue);
			}
			static void Prime(T &lValue)
			{
				return std::atanh(lValue);
			}
		};
	}

	template<class T> using MatrixDynamic = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	template<class T, class Lambda> void MatrixIterator(MatrixDynamic<T> &i, Lambda lFunc)
	{
		T *lData = i.data();
		Size lElements = i.rows()*i.cols();
		for (Size lE = 0; lE < lElements; ++lE){ lFunc(*lData++); }
	}

	template<class T, class NeuronFunction = NeuronFunctions::Sigmoid<T> > class FeedForwardDense
	{
	public:
		FeedForwardDense(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims);
		void Resize(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims);
		void Randomise();

		void Train(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);
		T CalculateCostValues(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);
		std::vector<MatrixDynamic<T>> CalculateCostGradient(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);
		MatrixDynamic<T> Process(const MatrixDynamic<T> &lInputValues);
		MatrixDynamic<T> ProcessFull(const MatrixDynamic<T> &lInputValues);

		Bool TestCostFunction(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);
	protected:		
		
		MatrixDynamic<T> ProcessLayer(const MatrixDynamic<T> &lInputs, const MatrixDynamic<T> &lWeights);

		MatrixDynamic<T> mInputWeights;
		std::vector<MatrixDynamic<T>> mLayerWeights;		
		std::vector<MatrixDynamic<T>> mA;
		std::vector<MatrixDynamic<T>> mZ;
		MatrixDynamic<T> mOutputWeights;

		Size mInputDims;
		Size mLayers;
		Size mLayerSize;
	};

}
