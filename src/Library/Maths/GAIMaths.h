
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

	template<class T> void Write(const MatrixDynamic<T> &lMatrix, std::ostream &lStream)
	{
		UInt64 lRows = sizeof(T);
		lStream.write((char*)&lRows, sizeof(UInt64));
		lRows = lMatrix.rows();
		lStream.write((char*)&lRows, sizeof(UInt64));
		lRows = lMatrix.cols();
		lStream.write((char*)&lRows, sizeof(UInt64));
		lStream.write((char*)lMatrix.data(), sizeof(T)*lMatrix.rows()*lMatrix.cols());
	}

	template<class T> void Read(MatrixDynamic<T> &lMatrix, std::istream &lStream)
	{
		UInt64 lRows;
		UInt64 lCols;
		lStream.read((char*)&lRows, sizeof(UInt64));
		if (lRows != sizeof(T)) { printf("ERROR : Matrix Type is wrong size\n"); }
		
		lStream.read((char*)&lRows, sizeof(UInt64));
		lStream.read((char*)&lCols, sizeof(UInt64));

		lMatrix.resize(Int32(lRows), Int32(lCols));
		lStream.read((char*)lMatrix.data(), sizeof(T)*lRows*lCols);
	}

	template<class T> void Write(const std::vector<T> &lVector, std::ostream &lStream)
	{
		UInt64 lRows = sizeof(T);
		lStream.write((char*)&lRows, sizeof(UInt64));
		lRows = lVector.size();
		lStream.write((char*)&lRows, sizeof(UInt64));
		for (auto &i : lVector)
		{
			Write(i, lStream);
		}
	}

	template<class T> void Read(std::vector<T> &lVector, std::istream &lStream)
	{
		UInt64 lRows;
		lStream.read((char*)&lRows, sizeof(UInt64));
		if (lRows != sizeof(T)) { printf("ERROR : Matrix Type is wrong size\n"); }
		
		lStream.read((char*)&lRows, sizeof(UInt64));
		lVector.resize(UInt32(lRows));
		for (auto &i : lVector)
		{
			Read(i, lStream);
		}
	}

	template<class T> void Write(const T &lVector, std::ostream &lStream)
	{
		lStream.write((char*)&lVector, sizeof(T));
	}

	template<class T> void Read(T &lVector, std::istream &lStream)
	{
		lStream.read((char*)&lVector, sizeof(T));
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
