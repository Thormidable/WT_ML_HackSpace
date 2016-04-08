
#pragma once

namespace GAI
{
	template<class T> using MatrixDynamic = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	template<class T, class Lambda> void MatrixIterator(MatrixDynamic<T> &i, Lambda lFunc)
	{
		T *lData = i.data();
		Size lElements = i.rows()*i.cols();
		for (Size lE = 0; lE < lElements; ++lE){ lFunc(*lData++); }
	}
	template<class T> void DisplayMatrix(const MatrixDynamic<T> &lMatrix, std::string lsName)
	{
		printf("Matrix : %s\n", lsName.c_str());
		for (Size lRows = 0; lRows < Size(lMatrix.rows()); ++lRows)
		{
			printf("\t");
			for (Size lColumns = 0; lColumns < Size(lMatrix.cols()); ++lColumns)
			{
				printf("%f", lMatrix(lRows, lColumns));
				if (lColumns < Size(lMatrix.cols() - 1)) printf(" , ");
				else printf("\n");
			}
		}
	}


	template<class T,class U, class Lambda> static inline void IterateMatrixList(std::vector<MatrixDynamic<T>> &lList, U *g, Lambda lFunc)
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

	template<class T> void NormaliseInputDataSet(MatrixDynamic<T> &lMatrix)
	{
		for (Size lCol = 0; lCol < Size(lMatrix.cols()); ++lCol)
		{
			lMatrix.middleCols<1>(lCol) /= lMatrix.middleCols<1>(lCol).maxCoeff();
		}
	}
}
