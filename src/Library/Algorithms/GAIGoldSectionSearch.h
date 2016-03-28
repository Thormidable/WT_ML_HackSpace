
#pragma once

namespace GAI
{
	
		/*! 
		\brief Collection of Optimisation Functions
		*/
		namespace Optimisation
		{

			/*!
			\brief This will find the minimum for the single dimensional function lpFunction between the bounds lfLowBound and lfUpperBound with a maximum error of lfAccuracy Required.
			\param lpFunction This is the function to evaluate ( F(x) ).
			\param lfLowBound This is the bound on the lowest x value to consider feasible.
			\param lfUpperBound This is the bound on the Highest x value to consider feasible.
			\param lfAccuracyRequired sets the exit condition. Once the gap in bounds is reduced to less than lfAccuracyRequired the algorithm will return.
			\return This will return the x position of the minimum value of y for y = lpFunction(x).

			*/
			template<class tData,typename tFunction> tData GoldSectionSearch(tFunction lpFunction, tData lfLowBound, tData lfUpperBound, tData lfAccuracyRequired)
			{
				const tData lfcPhi = static_cast<tData>((1.0 + sqrt(5.0)) / 2.0);
				const tData lfcResPhi = 2 - lfcPhi;
				tData lfCentrePoint = lfLowBound + lfcResPhi*(lfUpperBound-lfLowBound);

				tData x;
				tData Fx;
				tData Fb;

				Fb = lpFunction(lfCentrePoint);
				do
				{
					Bool lbGreaterHalf = (lfUpperBound - lfCentrePoint > lfCentrePoint - lfLowBound);
					if (lbGreaterHalf) x = lfCentrePoint + lfcResPhi * (lfUpperBound - lfCentrePoint);
					else x = lfCentrePoint - lfcResPhi * (lfCentrePoint - lfLowBound);

					if (abs(lfUpperBound - lfLowBound) < lfAccuracyRequired * (abs(lfCentrePoint) + abs(x))) 
						return (lfUpperBound + lfLowBound) / 2; 

					Fx = lpFunction(x);

					if(abs(Fx) == std::numeric_limits<tData>::infinity() || Fx!=Fx )
					{
//						GAI::LogError("Function Evaluation Returned a Non Real Value");
						return std::numeric_limits<tData>::infinity();
					}

					if (Fx < Fb) 
					{
						if (lbGreaterHalf) lfLowBound=lfCentrePoint;
						else lfUpperBound=lfCentrePoint;
						lfCentrePoint=x;
						//b is now x so Fb is Fx
						Fb=Fx;
					}
					else 
					{
						if(Fx == Fb)
						{
							//Inverted Section Comparison.
							if (lbGreaterHalf) x = lfCentrePoint - lfcResPhi * (lfCentrePoint - lfLowBound);
							else x = lfCentrePoint + lfcResPhi * (lfUpperBound - lfCentrePoint);

							if (abs(lfUpperBound - lfLowBound) < lfAccuracyRequired * (abs(lfCentrePoint) + abs(x))) return (lfUpperBound + lfLowBound) / 2; 

							Fx = lpFunction(x);

							if (Fx < Fb) 
							{
								if (lbGreaterHalf) lfUpperBound=lfCentrePoint;
								else lfLowBound=lfCentrePoint;
								lfCentrePoint=x;
								//b is now x so Fb is Fx
								Fb=Fx;
							}
							else
							{
								if(Fx==Fb)
								{
									return Fb;
								}

								if (lbGreaterHalf) lfLowBound=x;
								else lfUpperBound=x;
							}
						}
						else
						{
							//Fb Is still current.
							if (lbGreaterHalf) lfUpperBound=x;
							else lfLowBound=x; 
						}
					}
				}while(1);
			};
		}

}
