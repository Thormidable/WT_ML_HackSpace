
#pragma once

namespace GAI
{
	/*!
	\brief Namespace containing functions for fitting functions etc.
	*/
	namespace Base
	{

		/*! 
		\brief Class for storing a linear Polynomial Function
		*/
		template<UInt8 liTerms,class tD=Float32> class tPolynomialFunction
		{
			tD mpTerms[liTerms];
		public:
			tPolynomialFunction(){};
			tPolynomialFunction(const tPolynomialFunction &lOther)
			{
				std::copy(lOther.mpTerms, lOther.mpTerms + liTerms, mpTerms);
			};
			tPolynomialFunction(tPolynomialFunction &&lOther)
			{
				std::copy(lOther.mpTerms, lOther.mpTerms + liTerms, mpTerms);
			};
			tD FindFx(tD lfXValue)
			{
				tD lfPoly=1.0;
				tD lfAns=0.0f;		
				for(UInt8 liTerm=0;liTerm<liTerms;++liTerm)
				{
					lfAns+=lfPoly*mpTerms[liTerm];
					lfPoly*=lfXValue;
				}
				return lfAns;
			};

			tD FindFxShape(tD lfXValue)
			{
				tD lfPoly=lfXValue;
				tD lfAns=0.0f;		
				for(UInt8 liTerm=1;liTerm<liTerms;++liTerm)
				{
					lfAns+=lfPoly*mpTerms[liTerm];
					lfPoly*=lfXValue;
				}
				return lfAns;
			};

			tD GetTerm(UInt8 liTerm){return mpTerms[liTerm];};
			void SetTerm(UInt8 liTerm,tD lfTerm){mpTerms[liTerm]=lfTerm;};
		};


		/*!
		\brief Will Perform Least Square Residual Curve Fitting using a First Order Polynomial
		\param liPoints is the Number of Data Points.
		\param lpX is the X Value of Each Point.
		\param lpFx is the F(x) Value for each X Point.
		\return Will return a pointer (if solution can be found ) to a tPolynomialFunction<2,tD>.
		The Polynomial returned will take X Values and return F(x) Values.
		*/
		template<class tD> tPolynomialFunction<2, tD> LeastSquaresRegresion1(Size liPoints, tD *lpFx, tD *lpX);

		/*!
		\brief Will Perform Least Square Residual Curve Fitting using a First Order Polynomial assuming that the x axis is the location of the values in the array
		\param liPoints is the Number of Data Points.
		\param lpFx is the F(x) Value for each X Point.
		\return Will return a pointer (if solution can be found ) to a tPolynomialFunction<2,tD>.
		The Polynomial returned will take X Values and return F(x) Values.
		*/
		template<class tD> tPolynomialFunction<2, tD> LeastSquaresRegresion1(Size liPoints, tD *lpFx);

		/*!
		\brief Will Perform Least Square Residual Curve Fitting using a Second Order Polynomial
		\param liPoints is the Number of Data Points.
		\param lpX is the X Value of Each Point.
		\param lpFx is the F(x) Value for each X Point.
		\return Will return a pointer (if solution can be found ) to a tPolynomialFunction<2,tD>.
		The Polynomial returned will take X Values and return F(x) Values.
		*/
		template<class tD> tPolynomialFunction<3, tD> LinearLeastSquaresRegresion2(Size liPoints, tD *lpX, tD *lpFx);
		/*!
		\brief Will Perform Least Square Residual Curve Fitting using a Second Order Polynomial
		\param liPoints is the Number of Data Points.
		\param lpFx is the F(x) Value for each X Point.
		\return Will return a pointer (if solution can be found ) to a tPolynomialFunction<2,tD>.
		The Polynomial returned will take X Values and return F(x) Values. Assumes that value of each point is it's index
		*/
		template<class tD> tPolynomialFunction<3, tD> LinearLeastSquaresRegresion2(Size liPoints, tD *lpFx);


	};
};