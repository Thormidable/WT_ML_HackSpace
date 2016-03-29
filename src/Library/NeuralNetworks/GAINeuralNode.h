
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

		template<class T> class Linear
		{
		public:
			static inline void Function(T &lValue)
			{
				(void)lValue;
			}
			static inline void Prime(T &lValue)
			{
				lValue = 1;
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

	template<class T, class NeuronFunction = NeuronFunctions::Sigmoid<T> > class FeedForwardDense
	{
	public:
		FeedForwardDense(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims);
		void Resize(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims);
		void Randomise();

		void SetRegularisation(T lReg){ mRegularisationFactor = lReg; };
		T GetRegularisation(){ return mRegularisationFactor; }

		void Train(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);
		T CalculateCostValues(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);
		std::vector<MatrixDynamic<T>> CalculateCostGradient(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);
		MatrixDynamic<T> Process(const MatrixDynamic<T> &lInputValues);
		MatrixDynamic<T> ProcessFull(const MatrixDynamic<T> &lInputValues);

		T UpdateNeuronWeights(const std::vector<MatrixDynamic<T>> &lGradients, T lfScale);
		Bool TestCostFunction(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected);

		void MoveWeightsToArray(lbfgsfloatval_t *);
		void MoveArrayToWeights(const lbfgsfloatval_t *);

		inline MatrixDynamic<T> GetResultMatrix(){ return mOutputValues; }

		T CalculateCostValuesFromResults(const MatrixDynamic<T> &lExpected);

		template<class U,class Lambda> static inline void IterateMatrixList(std::vector<MatrixDynamic<T>> &lList, U *g, Lambda lFunc)
		{
			for (auto &i : lList)
			{
				for (Size lRow = 0; lRow < Size(i.rows()); ++lRow)
				{
					for (Size lCol = 0; lCol < Size(i.cols()); ++lCol)
					{
						lFunc(*g++, i(lRow, lCol));
					}
				}
			};
		}
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
