
#pragma once

namespace GAI
{
	/*!
	*\brief This is my Optimised Radix sort for sorting arrays of data.

	*/
	template<class tScale> void RadixSort(UInt8 *lpList,UInt8 *lpBucket,tScale liSize)	
	{
		tScale liPos;
		tScale liSum,liOld;
		UInt8 *lpSwitch;
		tScale lpStarts[0x100];

		//Reset the bucket sizes
		memset(lpStarts,0,sizeof(tScale)*0x100);
		//Find the Size of each bucket;
		for(liPos=0;liPos<liSize;++liPos) {++lpStarts[lpList[liPos] ];}

		//Convert bucket sizes to bucket starts.
		liSum=liOld=0;
		for(liPos=0;liPos<0x100;++liPos)
		{
			liSum+=lpStarts[liPos];
			lpStarts[liPos]=liOld;
			liOld=liSum;
		}

		//Fill the Buckets.
		for(liPos=0;liPos<liSize;++liPos){lpBucket[lpStarts[lpList[liPos] ]++ ]=lpList[liPos];}

		lpSwitch=lpBucket;
		lpBucket=lpList;
		lpList=lpSwitch;

		memcpy(lpList,lpBucket,sizeof(UInt8)*liSize);
	}

	/// Radix Sort. See RadixStort(int*,int*,unsigned int).
	template<class tScale> void RadixSort(UInt16 *lpList,UInt16 *lpBucket,tScale liSize)	
	{
		UInt8 liPass;
		tScale liPos;
		tScale liSum,liOld;
		UInt16 *lpSwitch;
		tScale lpStarts[0x100];

		for(liPass=0;liPass<sizeof(UInt16)*8;liPass+=8)
		{
			//Reset the bucket sizes
			memset(lpStarts,0,sizeof(tScale)*0x100);
			//Find the Size of each bucket;
			for(liPos=0;liPos<liSize;++liPos) {++lpStarts[(lpList[liPos]>>liPass)&0xff];}

			//Convert bucket sizes to bucket starts.
			liSum=liOld=0;
			for(liPos=0;liPos<0x100;++liPos)
			{
				liSum+=lpStarts[liPos];
				lpStarts[liPos]=liOld;
				liOld=liSum;
			}

			//Fill the Buckets.
			for(liPos=0;liPos<liSize;++liPos){lpBucket[lpStarts[(lpList[liPos]>>liPass)&0xff]++]=lpList[liPos];}

			lpSwitch=lpBucket;
			lpBucket=lpList;
			lpList=lpSwitch;
		}
	}

	/// Thisis a radix sort for unsigned integers. See RadixStort(UInt8*,UInt8*,unsigned int).
	template<class tScale> void RadixSort(UInt32 *lpList,UInt32 *lpBucket,tScale liSize)
	{
		UInt8 liPass;
		tScale liPos;
		tScale liSum,liOld;
		UInt32 *lpSwitch;
		tScale lpStarts[0x100];

		for(liPass=0;liPass<sizeof(UInt32)*8;liPass+=8)
		{
			//Reset the bucket sizes
			memset(lpStarts,0,sizeof(tScale)*0x100);
			//Find the Size of each bucket;
			for(liPos=0;liPos<liSize;++liPos) {++lpStarts[(lpList[liPos]>>liPass)&0xff];}

			//Convert bucket sizes to bucket starts.
			liSum=liOld=0;
			for(liPos=0;liPos<0x100;++liPos)
			{
				liSum+=lpStarts[liPos];
				lpStarts[liPos]=liOld;
				liOld=liSum;
			}

			//Fill the Buckets.
			for(liPos=0;liPos<liSize;++liPos){lpBucket[lpStarts[(lpList[liPos]>>liPass)&0xff]++]=lpList[liPos];}

			lpSwitch=lpBucket;
			lpBucket=lpList;
			lpList=lpSwitch;
		}
	}

	/// Thisis a radix sort for unsigned integers. See RadixStort(UInt8*,UInt8*,unsigned int).
	template<class tScale> void RadixSort(UInt64 *lpList,UInt64 *lpBucket,tScale liSize)
	{
		UInt8 liPass;
		tScale liPos;
		tScale liSum,liOld;
		UInt64 *lpSwitch;
		tScale lpStarts[0x100];

		for(liPass=0;liPass<sizeof(UInt64)*8;liPass+=8)
		{
			//Reset the bucket sizes
			memset(lpStarts,0,sizeof(tScale)*0x100);
			//Find the Size of each bucket;
			for(liPos=0;liPos<liSize;++liPos) {++lpStarts[(lpList[liPos]>>liPass)&0xff];}

			//Convert bucket sizes to bucket starts.
			liSum=liOld=0;
			for(liPos=0;liPos<0x100;++liPos)
			{
				liSum+=lpStarts[liPos];
				lpStarts[liPos]=liOld;
				liOld=liSum;
			}

			//Fill the Buckets.
			for(liPos=0;liPos<liSize;++liPos){lpBucket[lpStarts[(lpList[liPos]>>liPass)&0xff]++]=lpList[liPos];}

			lpSwitch=lpBucket;
			lpBucket=lpList;
			lpList=lpSwitch;
		}
	}

	/*!
	\brief Wrapper function for RadixSorting ANY Data type. Assumes that each ascending bit is more significant than the previous.
	* \param lpList This points to the list of data to be sorted.
	* \param lpBucket This points to the array the function can use as a bucket array. Should be the same length as lpList.
	* \param liSize This is the Size of the data array the function will sort.
	* Radix sort works by dividing the values into 256 'buckets' which then allows the function to
	* literally insert the objects into the correct place in the list. This is repeated for each
	* byte from least significant byte to most (I think).
	*/
	template<class tType,class tScale> void RadixSort(tType *lpList,tType *lpBucket,tScale liSize)
	{
		if(sizeof(tType)==sizeof(UInt8)) RadixSort(reinterpret_cast<UInt8*>(lpList),reinterpret_cast<UInt8*>(lpBucket),liSize);
		if(sizeof(tType)==sizeof(UInt16)) RadixSort(reinterpret_cast<UInt16*>(lpList),reinterpret_cast<UInt16*>(lpBucket),liSize);
		if(sizeof(tType)==sizeof(UInt32)) RadixSort(reinterpret_cast<UInt32*>(lpList),reinterpret_cast<UInt32*>(lpBucket),liSize);
		if(sizeof(tType)==sizeof(UInt64)) RadixSort(reinterpret_cast<UInt64*>(lpList),reinterpret_cast<UInt64*>(lpBucket),liSize);
	}
}
