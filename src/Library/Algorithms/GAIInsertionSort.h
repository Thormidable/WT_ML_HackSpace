
#pragma once

namespace GAI
{
	/*!
	\brief Insertion Sort - Only efficient for mostly sorted data.
	*/
	template<class tType, typename lfLamda> void InsertionSort(tType *arr, UInt32 max, lfLamda Comparison)
	{
		UInt32 i, j;
		tType index;
		for (i = 1; i < max; i++)
		{
			index = arr[i];
			j = i;
			while ((j > 0) && (Comparison(arr[j - 1], index)))
			{
				arr[j] = arr[j - 1];
				j = j - 1;
			}
			arr[j] = index;
		}

	};

	// INSERTION SORT
	template<class tType> void InsertionSort(tType *arr, UInt32 max)
	{
		UInt32 i, j;
		tType index;
		for(i=1; i < max;i++)
		{
			index  = arr[i];
			j = i;
			while((j > 0) && (arr[j-1] > index))
			{
				arr[j] = arr[j-1];
				j = j-1;
			}
			arr[j] = index;
		}

	}
}
