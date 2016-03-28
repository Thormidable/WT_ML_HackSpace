
#include "GAILibrary.h"

template<class tD> GAI::Base::tPolynomialFunction<2, tD> GAI::Base::LeastSquaresRegresion1(Size liPoints, tD *lpFx, tD *lpX)
		{
			tPolynomialFunction<2, tD> lpPoly;
			if(liPoints>2)
			{
				tD lfPoints = static_cast<tD>(liPoints);
				tD lfSumX=0.0f;
				tD lfSumY=0.0f;
				tD lfSumXY=0.0f;
				tD lfSumX2=0.0f;
				for (Size i = 0; i<liPoints; i++)
				{			
					lfSumX+=lpX[i];
					lfSumY+=lpFx[i];
					lfSumXY+=lpX[i]*lpFx[i];
					lfSumX2+=lpX[i]*lpX[i];			
				}
				lpPoly.SetTerm(0,(lfSumY*lfSumX2-lfSumX*lfSumXY)/(lfPoints*lfSumX2-lfSumX*lfSumX));
				lpPoly.SetTerm(1,(lfSumXY-lfSumX*lfSumY)/(lfSumX2-lfSumX*lfSumX));
			}
			return lpPoly;
		}

		/*!
		\brief Will Perform Least Square Residual Curve Fitting using a First Order Polynomial assuming that the x axis is the location of the values in the array
		\param liPoints is the Number of Data Points.
		\param lpFx is the F(x) Value for each X Point.
		\return Will return a pointer (if solution can be found ) to a tPolynomialFunction<2,tD>.
		The Polynomial returned will take X Values and return F(x) Values.
		*/
template<class tD> GAI::Base::tPolynomialFunction<2, tD> GAI::Base::LeastSquaresRegresion1(Size liPoints, tD *lpFx)
		{
			tPolynomialFunction<2, tD> lpPoly;
			if (liPoints>2)
			{
				tD lfPoints = static_cast<tD>(liPoints);
				tD lfSumX = 0.0f;
				tD lfSumY = 0.0f;
				tD lfSumXY = 0.0f;
				tD lfSumX2 = 0.0f;
				for (Size i = 0; i<liPoints; i++)
				{
					lfSumX += i;
					lfSumY += lpFx[i];
					lfSumXY += i * lpFx[i];
					lfSumX2 += i * i;
				}

				lpPoly.SetTerm(0, (lfSumY*lfSumX2 - lfSumX*lfSumXY) / (lfPoints*lfSumX2 - lfSumX*lfSumX));
				lpPoly.SetTerm(1, (lfSumXY - lfSumX*lfSumY) / (lfSumX2 - lfSumX*lfSumX));
			}
			return lpPoly;
		}

		/*!
		\brief Will Perform Least Square Residual Curve Fitting using a Second Order Polynomial
		\param liPoints is the Number of Data Points.
		\param lpX is the X Value of Each Point.
		\param lpFx is the F(x) Value for each X Point.
		\return Will return a pointer (if solution can be found ) to a tPolynomialFunction<2,tD>.
		The Polynomial returned will take X Values and return F(x) Values.
		*/
template<class tD> GAI::Base::tPolynomialFunction<3, tD> GAI::Base::LinearLeastSquaresRegresion2(Size liPoints, tD *lpX, tD *lpFx)
		{
			if(liPoints>3)
			{	
				tD lfPoints=static_cast<tD>(liPoints);
				tD lfSumX=0.0f;
				tD lfSumY=0.0f;
				tD lfSumXY=0.0f;
				tD lfSumX2=0.0f;
				tD lfSumX3=0.0f;
				tD lfSumX2Y=0.0f;
				tD lfSumX4=0.0f;
				for (Size i = 0; i<liPoints; i++)
				{			
					lfSumX+=lpX[i];
					lfSumY+=lpFx[i];
					lfSumXY+=lpX[i]*lpFx[i];
					lfSumX2+=lpX[i]*lpX[i];
					lfSumX3+=lpX[i]*lpX[i]*lpX[i];
					lfSumX4+=lpX[i]*lpX[i]*lpX[i]*lpX[i];
					lfSumX2Y+=lpX[i]*lpX[i]*lpFx[i];
				}

				tD SXX=lfSumX2-(lfSumX*lfSumX/lfPoints);
				tD SXY=lfSumXY-(lfSumX*lfSumY/lfPoints);
				tD SXX2=lfSumX3-(lfSumX3)-(lfSumX*lfSumX2/lfPoints);
				tD SX2Y=lfSumX2Y-(lfSumX2*lfSumY/lfPoints);
				tD SX2X2=lfSumX4-(lfSumX2*lfSumX2/lfPoints);

				tPolynomialFunction<3,tD> lpPoly;
				tD a=(SX2Y*SXX-SXY*SXX2)/(SXX*SX2X2-SXX2*SXX2);
				tD b=(SXY*SX2X2-SX2Y*SXX2)/(SXX*SX2X2-SXX2*SXX2);
				lpPoly.SetTerm(2,a);
				lpPoly.SetTerm(1,b);
				lpPoly.SetTerm(0,(lfSumY/lfPoints) - (b*lfSumX/lfPoints)-(a*lfSumX2/lfPoints));

				return lpPoly;
			}

			return tPolynomialFunction<3, tD>();
		}
		/*!
		\brief Will Perform Least Square Residual Curve Fitting using a Second Order Polynomial
		\param liPoints is the Number of Data Points.
		\param lpFx is the F(x) Value for each X Point.
		\return Will return a pointer (if solution can be found ) to a tPolynomialFunction<2,tD>.
		The Polynomial returned will take X Values and return F(x) Values. Assumes that value of each point is it's index
		*/
template<class tD> GAI::Base::tPolynomialFunction<3, tD> GAI::Base::LinearLeastSquaresRegresion2(Size liPoints, tD *lpFx)
{
	if (liPoints > 3)
	{
		tD lfPoints = static_cast<tD>(liPoints);
		tD lfSumX = 0.0f;
		tD lfSumY = 0.0f;
		tD lfSumXY = 0.0f;
		tD lfSumX2 = 0.0f;
		tD lfSumX3 = 0.0f;
		tD lfSumX2Y = 0.0f;
		tD lfSumX4 = 0.0f;
		for (Size i = 0; i < liPoints; i++)
		{
			lfSumX += i;
			lfSumY += lpFx[i];
			lfSumXY += i * lpFx[i];
			lfSumX2 += i * i;
			lfSumX3 += i * i * i;
			lfSumX4 += i * i * i * i;
			lfSumX2Y += i * i * lpFx[i];
		}

		tD SXX = lfSumX2 - (lfSumX*lfSumX / lfPoints);
		tD SXY = lfSumXY - (lfSumX*lfSumY / lfPoints);
		tD SXX2 = lfSumX3 - (lfSumX3)-(lfSumX*lfSumX2 / lfPoints);
		tD SX2Y = lfSumX2Y - (lfSumX2*lfSumY / lfPoints);
		tD SX2X2 = lfSumX4 - (lfSumX2*lfSumX2 / lfPoints);

		tPolynomialFunction<3, tD> lpPoly;
		tD a = (SX2Y*SXX - SXY*SXX2) / (SXX*SX2X2 - SXX2*SXX2);
		tD b = (SXY*SX2X2 - SX2Y*SXX2) / (SXX*SX2X2 - SXX2*SXX2);
		lpPoly.SetTerm(2, a);
		lpPoly.SetTerm(1, b);
		lpPoly.SetTerm(0, (lfSumY / lfPoints) - (b*lfSumX / lfPoints) - (a*lfSumX2 / lfPoints));

		return lpPoly;
	}

	return tPolynomialFunction<3, tD>();
}

#define GAI_REGRESION_INSTANCING(tD) \
template GAI::Base::tPolynomialFunction<2, tD> GAI::Base::LeastSquaresRegresion1(Size liPoints, tD *lpFx, tD *lpX);\
template GAI::Base::tPolynomialFunction<2, tD> GAI::Base::LeastSquaresRegresion1(Size liPoints, tD *lpFx);\
template GAI::Base::tPolynomialFunction<3, tD> GAI::Base::LinearLeastSquaresRegresion2(Size liPoints, tD *lpX, tD *lpFx);\
template GAI::Base::tPolynomialFunction<3, tD> GAI::Base::LinearLeastSquaresRegresion2(Size liPoints, tD *lpFx);

GAI_EXPAND(GAI_FOR_EACH_FLOAT_TYPE(GAI_REGRESION_INSTANCING));
