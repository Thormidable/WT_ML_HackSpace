
#include "GAILibrary.h"

using namespace GAI;

template<class T> GAI::FeedForwardDense<T>::FeedForwardDense(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims)
{
	Resize(lInputDims, lLayerSize, lLayers, lOutputDims);
}

template<class T> void GAI::FeedForwardDense<T>::Resize(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims)
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

template<class T> void GAI::FeedForwardDense<T>::Randomise()
{
	if (mLayerWeights.size() == 0) throw;

	for (auto &i : mLayerWeights)
	{
		T *lData = i.data();
		Size lElements = i.rows()*i.cols();
		for (Size lE = 0; lE < lElements; ++lE){*lData++ = Random<T>(0, 1.0);}
	}
}

template<class T> MatrixDynamic<T> GAI::FeedForwardDense<T>::Process(const MatrixDynamic<T> &lInputValues)
{
	MatrixDynamic<T> lResult = lInputValues;
	ProcessLayer(lResult, mLayerWeights[0]);

	for (Size lLayer = 1; lLayer < mLayerWeights.size(); ++lLayer)
	{
		ProcessLayer(lResult, mLayerWeights[lLayer]);
	}

	return lResult.colwise().sum();
}

template<class T> void GAI::FeedForwardDense<T>::ProcessLayer(MatrixDynamic<T> &lInputs, const MatrixDynamic<T> &lWeights)
{	
	Size lElements = lWeights.cols()*lInputs.rows();
	lInputs = lInputs*lWeights;
	T *lPtr = lInputs.data();
	for (Size r = 0; r < lElements; ++r){ *lPtr++ = Sigmoid(*lPtr); }	
}

template<class T> T GAI::Sigmoid(T lFunc)
{
	return 1 / (1 + std::exp(-lFunc));
}
template<class T> T GAI::TanH(T lFunc)
{
	return std::tanh(lFunc);
}

#define INSTANCE_NEURON_FUNCTIONS(T) \
template T GAI::Sigmoid(T lFunc);\
template T GAI::TanH(T lFunc);\
template class GAI::FeedForwardDense<T>;

GAI_EXPAND(GAI_FOR_EACH_FLOAT_TYPE(INSTANCE_NEURON_FUNCTIONS))
