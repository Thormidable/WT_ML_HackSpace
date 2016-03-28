
#pragma once

namespace GAI
{

	template<class T> using NeuronFunction = T(*)(T);

	template<class T> T Sigmoid(T);
	template<class T> T TanH(T);

	template<class T> using MatrixDynamic = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

	template<class T> class FeedForwardDense
	{
	public:
		FeedForwardDense(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims);
		void Resize(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims);
		void Randomise();

		MatrixDynamic<T> Process(const MatrixDynamic<T> &lInputValues);
	protected:		
		
		void ProcessLayer(MatrixDynamic<T> &lInputs, const MatrixDynamic<T> &lWeights);

		MatrixDynamic<T> mInputWeights;
		std::vector<MatrixDynamic<T>> mLayerWeights;		
		MatrixDynamic<T> mOutputWeights;

		Size mInputDims;
		Size mLayers;
		Size mLayerSize;
	};

}
