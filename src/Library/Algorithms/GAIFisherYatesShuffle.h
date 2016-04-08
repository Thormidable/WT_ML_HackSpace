
#include "GAILibrary.h"

namespace GAI
{
	template<class T> void FisherYatesShuffle(T *lData, Size lElements)
	{
		for (Size i = lElements; --i > 1;)
		{			
			std::Swap(lData[Random<Size>(0, i)], lData[i]);
		}
	}

	template<class T> void FisherYatesShuffleLinked(MatrixDynamic<T> &lData, MatrixDynamic<T> &lExpected)
	{
		Eigen::Matrix<T, 1, Eigen::Dynamic> lRow;
		Eigen::Matrix<T, 1, Eigen::Dynamic> lExpectedRow;
		for (Size i = lData.rows(); --i > 1;)
		{
			Size lRowID = Random<Size>(0, i);

			lRow = lData.middleRows<1>(i);
			lData.middleRows<1>(i) = lData.middleRows<1>(lRowID);
			lData.middleRows<1>(lRowID) = lRow;

			lExpectedRow = lExpected.middleRows<1>(i);
			lExpected.middleRows<1>(i) = lExpected.middleRows<1>(lRowID);
			lExpected.middleRows<1>(lRowID) = lExpectedRow;
		}
	}


	template<class T> void FisherYatesShuffleRows(MatrixDynamic<T> &lData)
	{
		Eigen::Matrix<T,1,Eigen::Dynamic> lRow;
		for (Size i = lData.rows(); --i > 1;)
		{
			Size lRowID = Random<Size>(0, i);
			lRow = lData.middleRows<1>(i);
			lData.middleRows<1>(i) = lData.middleRows<1>(lRowID);
			lData.middleRows<1>(lRowID) = lRow;
		}
	}

	template<class T> void FisherYatesShuffleCols(MatrixDynamic<T> &lData)
	{
		Eigen::Matrix<T, Eigen::Dynamic,1> lCol;
		for (Size i = lData.cols(); --i > 1;)
		{
			Size lColID = Random<Size>(0, i);
			lCol = lData.middleCols<1>(i);
			lData.middleCols<1>(i) = lData.middleCols<1>(lColID);
			lData.middleCols<1>(lColID) = lCol;
		}
	}
}