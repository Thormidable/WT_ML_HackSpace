
#pragma once

namespace GAI
{
	/*!
	\brief Bubble Sort - Only efficient for mostly sorted data.
	*/
	template<class tType,typename lfLamda> void BubbleSort(tType *arr, UInt32 max,lfLamda Comparison)
	{
		tType tmp;
		for(UInt32 i=0;i<max;i++)
		{
			for(UInt32 x=0; x<max-1-i; x++)
			{
				if(Comparison(arr[x],arr[x+1]))
				{
					//r.push_back(rnd);
					tmp = arr[x];
					arr[x] = arr[x+1];
					arr[x+1] = tmp;
				}
			}
		}
	}

	// BUBBLE SORT
	template<class tType> void BubbleSort(tType *arr, UInt32 max)
	{
		tType tmp;
		for(UInt32 i=0;i<max;i++)
		{
			for(UInt32 x=0; x<max-1-i; x++)
			{
				if(arr[x] > arr[x+1])
				{
					//r.push_back(rnd);
					tmp = arr[x];
					arr[x] = arr[x+1];
					arr[x+1] = tmp;
				}
			}
		}
	}
}
