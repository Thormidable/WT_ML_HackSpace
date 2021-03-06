
#include "GAILibrary.h"
//
//using namespace GAI;
//
//template<class T, class NeuronFunction> GAI::FeedForwardEnergised<T, NeuronFunction>::FeedForwardEnergised(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims) : FeedForwardDense<T, NeuronFunction>(lInputDims,lLayerSize,lLayers,lOutputDims)
//{
//	
//}
//
//template<class T, class NeuronFunction> void GAI::FeedForwardDense<T,NeuronFunction>::Resize(Size lInputDims, Size lLayerSize, Size lLayers, Size lOutputDims)
//{
//	if (lLayerSize == 0) printf("ERROR. Numbers of Layers too low. for Forward Dense");
//
//	if (mLayerSize != lLayerSize || lLayers != mLayers || (mLayerWeights.size() == 0 || mLayerWeights[0].rows() != lInputDims) || mLayerWeights.back().cols() != lOutputDims)
//	{
//		mLayerSize = lLayerSize;
//		mLayers = lLayers;
//		mLayerWeights.resize(lLayers+1);
//		mLayerWeights[0].resize(lInputDims, lLayerSize);
//
//		for (Size lLayer = 1; lLayer < lLayers; ++lLayer){mLayerWeights[lLayer].resize(lLayerSize, lLayerSize);}
//
//		mLayerWeights.back().resize(lLayerSize, lOutputDims);		
//	}
//
//	Randomise();
//
//}
//
//template<class T, class NeuronFunction> void GAI::FeedForwardDense<T,NeuronFunction>::Randomise()
//{
//	if (mLayerWeights.size() == 0) throw;
//
//	for (auto &i : mEnergiserWeights)
//	{		
//		MatrixIterator(i, [](T &lData){lData = Random<T>(0.0, 1.0); });
//	}
//}
//
//template<class T, class NeuronFunction> void GAI::FeedForwardDense<T, NeuronFunction>::MoveArrayToWeights(const lbfgsfloatval_t *xCursor)
//{
//	IterateMatrixList<T,const lbfgsfloatval_t>(mLayerWeights, xCursor, [](const lbfgsfloatval_t &lOne, T &lMat){lMat = T(lOne); });
//}
//
//template<class T, class NeuronFunction> void GAI::FeedForwardDense<T, NeuronFunction>::MoveWeightsToArray(lbfgsfloatval_t *xCursor)
//{	
//	IterateMatrixList<T,lbfgsfloatval_t>(mLayerWeights, xCursor, [](lbfgsfloatval_t &lOne, T &lMat){lOne = lbfgsfloatval_t(lMat); });
//}
//
//template<class T, class NeuronFunction> void GAI::FeedForwardDense<T, NeuronFunction>::Train(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected)
//{
//	T lfLastMove = 10.0f;
//	
//	printf("Starting with a Cost Value of : %f\n", CalculateCostValues(lInputValues, lExpected));
//
//	Size lParameters = 0;
//	for (auto &i : mLayerWeights){lParameters += i.rows()*i.cols();}
//	
//	lbfgsfloatval_t fx;
//	lbfgsfloatval_t *x = lbfgs_malloc(lParameters);
//	lbfgs_parameter_t param;
//	if (x == NULL) {printf("ERROR: Failed to allocate a memory block for variables.\n");return ;}
//
//	MoveWeightsToArray(x);
//
//	lbfgs_parameter_init(&param);
//
//	NeuronTraining<T, NeuronFunction> lInstance(this,lInputValues,lExpected);
//
//	param.max_iterations = 1000;
//	param.m = 6;
//	param.past = 0;
//	param.orthantwise_c = 0;
//	param.ftol = 0.0000001;
//	param.gtol = 0.9;
//
//	lbfgs(lParameters, x, &fx, NeuronTrainingEvaluation<T, NeuronFunction>, NeuronTrainingProgress<T, NeuronFunction>, &lInstance, &param);
//	
//	MoveArrayToWeights(x);
//	lbfgs_free(x);
//	
//	printf("Optimised To a Cost Value of : %f in %u iterations\n",CalculateCostValues(lInputValues, lExpected),lInstance.mIterations);
//}
//
//template<class T, class NeuronFunction> T GAI::FeedForwardDense<T, NeuronFunction>::CalculateCostValues(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected)
//{
//	ProcessFull(lInputValues);
//	return CalculateCostValuesFromResults(lExpected);
//}
//
//template<class T, class NeuronFunction> T GAI::FeedForwardDense<T, NeuronFunction>::CalculateCostValuesFromResults(const MatrixDynamic<T> &lExpected)
//{
//	auto lEval = (mOutputValues - lExpected);	
//
//	T lfRegulisation = T(0);
//	for (auto &i : mLayerWeights){ lfRegulisation += i.cwiseProduct(i).sum(); }
//	lfRegulisation *= (mRegularisationFactor / 2)*T(lExpected.rows());
//
//	return (lEval.cwiseProduct(lEval).array()).sum()*T(0.5) + lfRegulisation;
//}
//
//template<class T, class NeuronFunction> std::vector<MatrixDynamic<T>> GAI::FeedForwardDense<T, NeuronFunction>::CalculateCostGradient(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected)
//{
//	auto lYHat = ProcessFull(lInputValues);
//
//	MatrixDynamic<T> lD3;
//	std::vector<MatrixDynamic<T>> lDJdW;
//	lDJdW.resize(mLayerWeights.size());
//	lDJdW.back() = -(lExpected-lYHat);
//
//	for (Size lLayer = lDJdW.size()-1; lLayer > 0; --lLayer)
//	{
//		//Calculate sigmoidPrime of Input Energys
//		auto lFPrimeZ3 = mZ[lLayer];
//		MatrixIterator<T>(lFPrimeZ3, &NeuronFunction::Prime);
//
//		//Calculate ErrorFactor * sigmoidPrime(Input energies)
//		lD3 = lDJdW[lLayer].cwiseProduct(lFPrimeZ3);
//		lDJdW[lLayer] = (mA[lLayer - 1].transpose() * lD3) + (mLayerWeights[lLayer] * mRegularisationFactor*T(lInputValues.rows()));
//		if(lLayer > 1) lDJdW[lLayer - 1] = (lD3*mLayerWeights[lLayer].transpose()).eval();
//	}
//
//	auto lFPrimeZ3 = mZ[0];
//	MatrixIterator<T>(lFPrimeZ3, &NeuronFunction::Prime);
//
//	lD3 = (lD3*(mLayerWeights[1].transpose())).cwiseProduct(lFPrimeZ3);
//	lDJdW[0] = (lInputValues.transpose()*lD3) + (mLayerWeights[0] * mRegularisationFactor*T(lInputValues.rows()));
//
//	return lDJdW;
//}
//
//template<class T, class NeuronFunction> Bool GAI::FeedForwardDense<T, NeuronFunction>::TestCostFunction(const MatrixDynamic<T> &lInputValues, const MatrixDynamic<T> &lExpected)
//{
//	auto lfCalculatedGrads = CalculateCostGradient(lInputValues, lExpected);
//
//	T lEps = 0.0001f;
//
//	for (Size lLayer = mLayerWeights.size(); lLayer-- >0; )
//	{
//		auto lWeighted = mLayerWeights[lLayer];
//
//		for (Size lRow = 0; lRow < Size(mLayerWeights[lLayer].rows()); ++lRow)
//		{
//			for (Size lCol = 0; lCol < Size(mLayerWeights[lLayer].cols()); ++lCol)
//			{
//				mLayerWeights[lLayer](lRow, lCol) += lEps;
//				T lLoss1 = CalculateCostValues(lInputValues, lExpected);
//
//				mLayerWeights[lLayer] = lWeighted;
//				mLayerWeights[lLayer](lRow, lCol) -= lEps;
//				T lLoss2 = CalculateCostValues(lInputValues, lExpected);
//
//				mLayerWeights[lLayer] = lWeighted;
//
//				auto lfGrads = (lLoss1 - lLoss2) / (2 * lEps);
//
//				auto lDiffResult = (lfCalculatedGrads[lLayer](lRow, lCol) - lfGrads);
//				if (abs(lDiffResult) > 0.0000001f)
//					return false;
//			}
//		}
//	}
//	return true;
//}
//
//template<class T, class NeuronFunction> MatrixDynamic<T> GAI::FeedForwardDense<T,NeuronFunction>::Process(const MatrixDynamic<T> &lInputValues)
//{
//	mOutputValues = lInputValues;
//	mOutputValues *= mLayerWeights[0];
//	MatrixIterator(mOutputValues, &NeuronFunction::Function);	
//
//	for (Size lLayer = 1; lLayer < mLayerWeights.size(); ++lLayer)
//	{
//		mOutputValues*= mLayerWeights[lLayer];
//		MatrixIterator(mOutputValues, &NeuronFunction::Function);
//	}
//
//	return mOutputValues;
//}
//
//template<class T, class NeuronFunction> MatrixDynamic<T> GAI::FeedForwardDense<T, NeuronFunction>::ProcessFull(const MatrixDynamic<T> &lInputValues)
//{
//	mA.resize(mLayerWeights.size() - 1);
//	mZ.resize(mLayerWeights.size());
//
//	mA[0] = lInputValues;
//	for (Size lLayer = 0; lLayer < mA.size(); ++lLayer)
//	{
//		mA[lLayer] *= mLayerWeights[lLayer];
//		mZ[lLayer] = mA[lLayer];
//		MatrixIterator(mA[lLayer], &NeuronFunction::Function);		
//		if (lLayer < mA.size() - 1) mA[lLayer + 1] = mA[lLayer];
//	}
//
//	mZ.back() = mA.back()*mLayerWeights.back();
//	mOutputValues = mZ.back();
//	MatrixIterator(mOutputValues, &NeuronFunction::Function);
//	return mOutputValues;
//}
//
//#define INSTANCE_NEURON_FUNCTIONS(T) \
//template class GAI::FeedForwardDense<T,NeuronFunctions::Sigmoid<T> >;\
//template class GAI::FeedForwardDense<T,NeuronFunctions::Linear<T> >;
//
//GAI_EXPAND(GAI_FOR_EACH_FLOAT_TYPE(INSTANCE_NEURON_FUNCTIONS))
