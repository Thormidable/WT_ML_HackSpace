
#pragma once

namespace GAI
{
	/*!
	\brief Function for searching a list for a value in a binary fashion (half and check repeat).	
	\param lpData List of data to be sorted. 
	\param lfValue Value to be found in the list.
	\param liMin Lowest valid Index of the Array.
	\param liMax Highest valid Index of the Array
	\return Index of Value or Index of closest Value in the Array.
	This function requires a sorted ascending list - It will not sort the values itself. If unsure of whether sorted, use a sorting algorithm to sort before passing data to this function. 
	Data type of list also requires the operators < and > to work for the values.
	*/
	template<class tT, typename tLambda> Int32 BinarySearch(tT lfValue, Int32 liMin, Int32 liMax, tLambda ExtractValue)
	{
		if (liMin > liMax)
		{
			GAI::Maths::Swap(liMin, liMax);
		}

		Int32 liMid = liMin;
		//printf("lfValue : %f    ",lfValue);
		// continue searching while [liMin,liMax] is not empty
		while (liMax > liMin)
		{

			liMid = (liMin+liMax)/2;

			// determine which subarray to search
			if (ExtractValue(liMid) <  lfValue)
			{  // change min index to search upper subarray
				liMin = liMid + 1;
				//	printf("<");
			}
			else
			{
				if (ExtractValue(liMid) > lfValue)
				{// change max index to search lower subarray
					liMax = liMid - 1;
					//printf(">");
				}
				else
				{
					//printf("Found %u\n",liMid);
					// lfValue found at index liMid
					return liMid;
				}
			}
		}
		// lfValue not found
		//printf("Found %u\n",liMid);
		if (liMin>=0 && ExtractValue(liMin) == lfValue) return liMin;
		else return -1;
	}

	template<class tT> Int32 BinarySearch(tT lfValue,tT *lpData, Int32 liMin, Int32 liMax)
	{
		auto lLambda = [&](UInt32 liPos)
		{
			return lpData[liPos]; 
		};
		return BinarySearch(lfValue, liMin, liMax, lLambda);
	}


	/*!
	\brief Function for searching a list for a value in a binary fashion (half and check repeat).
	\param lpData List of data to be sorted.
	\param lfValue Value to be found in the list.
	\param liMin Lowest valid Index of the Array.
	\param liMax Highest valid Index of the Array
	\return Index of Value or Index of closest Value in the Array.
	This function requires a sorted descending list - It will not sort the values itself. If unsure of whether sorted, use a sorting algorithm to sort before passing data to this function.
	Data type of list also requires the operators < and > to work for the values.
	*/
	template<class tT, typename tLambda> Int32 BinarySearchInv(tT lfValue, Int32 liMin, Int32 liMax, tLambda ExtractValue)
	{
		if (liMin > liMax)
		{
			GAI::Maths::Swap(liMin, liMax);
		}

		Int32 liMid;
		// continue searching while [liMin,liMax] is not empty
		while (liMax > liMin)
		{
			liMid = (liMin + liMax) / 2;
			// determine which subarray to search
			if (ExtractValue(liMid) >  lfValue)
			{  
				liMin = liMid + 1;
			}
			else
			{
				if (ExtractValue(liMid) < lfValue) liMax = liMid - 1;
				else return liMid;
			}
		}

		if (liMin >= 0 && ExtractValue(liMin) == lfValue) return liMin;
		else return -1;
	}

	template<class tT> Int32 BinarySearchInv(tT lfValue, tT *lpData, Int32 liMin, Int32 liMax)
	{
		return BinarySearchInv(lfValue, liMin, liMax, [&](UInt32 liPos){return lpData[liPos]; });
	}
}
