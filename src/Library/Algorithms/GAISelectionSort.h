
#pragma once

namespace GAI
{
	/*!
	\brief Selection Sort 
	*/
	template<class tType,typename lfLamda> void SelectionSort(tType *arr, UInt32 max,lfLamda Comparison)
	{
		tType tmp;
		UInt32 min;
		for (UInt32 i = 0; i<max; i++)
		{
			min = i;
			for(UInt32 x=i; x<max; x++)
			{
				if(Comparison(arr[x],arr[min]))
				{
					min = x;
				}
			}
			tmp = arr[i];
			arr[i] = arr[min];
			arr[min] = tmp;
		}
	};

	// SELECTION SORT
	template<class tType> void SelectionSort(tType *arr, UInt32 max)
	{
		tType tmp;
		UInt32 min;
		for(UInt32 i=0;i<max;i++)
		{
			min = i;
			for(UInt32 x=i; x<max; x++)
			{
				if(arr[x] < arr[min])
				{
					min = x;
				}
			}
			tmp = arr[i];
			arr[i] = arr[min];
			arr[min] = tmp;
		}
	}
};
