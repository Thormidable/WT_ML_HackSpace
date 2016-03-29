
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
}
