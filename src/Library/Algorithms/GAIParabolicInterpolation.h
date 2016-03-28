
#pragma once

namespace GAI
{
	namespace Base
	{
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
			template<typename tFunction,class tData> tData SuccessiveParabolicInterpolation(tFunction lpFunction,tData lfLowBound, tData lfUpperBound, tData lfAccuracyRequired) 
			{
				tData lfMinSpread = 0.01;
				tData lpFx[3];
				tData lpX[3];
				lpX[0] = lfLowBound;
				lpX[1] = (lfUpperBound+lfLowBound)/2;
				lpX[2] = lfUpperBound;



				//Find Three Points
				tData lpFL =lpFx[0] = lpFunction(lfLowBound);
				lpFx[1] = lpFunction(lpX[1]);
				tData lpFU = lpFx[2] = lpFunction(lfUpperBound);

				Bool lbError=false;
				Int32 liCount = 100;
				do
				{
					lfMinSpread*=88.8;
					//Generate The parabola ( Ax^2+BX+C = Y)
					// C is not necessary for finding the minimum


					tData lfA = ((lpX[1]*lpX[1]-lpX[0]*lpX[0])*(lpX[2]-lpX[0])-(lpX[2]*lpX[2]-lpX[0]*lpX[0])*(lpX[1]-lpX[0]));

					//if true the parabola is now a straight line. The function has failed.
					if(lfA==0.0) 
						return lpX[1];

					lfA = ((lpFx[1] - lpFx[0])*(lpX[2]-lpX[0])-(lpFx[2] - lpFx[0])*(lpX[1]-lpX[0]))/lfA;

					if (lfA>0.0)
					{
						tData lfB = ( (lpFx[2]-lfA*lpX[2]*lpX[2])-(lpFx[1]-lfA*lpX[1]*lpX[1]) ) / (lpX[2] - lpX[1]);
						tData lfNew = -lfB/(2*lfA);

						//Find the half to insert the new value into.
						if(lfNew > lpX[1])
						{
							if(lfNew>lfUpperBound)
							{ 
								if(lbError)
									return lpX[1];
								lbError=true;
								lpX[1]=lpX[2];
								lpFx[1]=lpFx[2];

								lpX[2] = lfNew;
								lpFx[2] = lpFunction(lpX[2]);
							}
							else 
							{
								lbError=false;
								lpX[0]=lpX[1];
								lpFx[0]=lpFx[1];

								if(abs(lfNew-lpX[1]) <= lfAccuracyRequired*abs(lfNew)+(std::numeric_limits<tData>::epsilon)()*0.0001) 
									return lfNew;

								lpX[1]=lfNew;
								lpFx[1]=lpFunction(lpX[1]);

								if(abs(lpX[0]-lpX[1])<lfMinSpread) lpX[0] = lpX[1]-lfMinSpread;
							}
						}
						else
						{
							if(lfNew<lfLowBound)
							{ 
								if(lbError) 
									return lpX[1];
								lbError=true;
								lpX[1]=lpX[0];
								lpFx[1]=lpFx[0];

								lpX[0] = lfNew;
								lpFx[0] = lpFunction(lpX[0]);
							}
							else 
							{
								lbError=false;
								lpX[2]=lpX[1];
								lpFx[2] = lpFx[1];

								if(abs(lfNew-lpX[1]) <= lfAccuracyRequired*abs(lfNew)+(std::numeric_limits<tData>::epsilon)()*0.0001) 
									return lfNew;

								lpX[1]=lfNew;
								lpFx[1]=lpFunction(lpX[1]);

								if(abs(lpX[2]-lpX[1])<lfMinSpread) lpX[2] = lpX[1]-lfMinSpread;
							}
						}
					}
					else
					{
						// Have found a maxima.
						return GoldSectionSearch(lpFunction,lfLowBound,lfUpperBound,lfAccuracyRequired);
					}



				}while(--liCount);
				//	}while(1);

				if(!liCount) return GoldSectionSearch(lpFunction,lfLowBound,lfUpperBound,lfAccuracyRequired);
				else return lpX[1];
			};
		}
	}
}
