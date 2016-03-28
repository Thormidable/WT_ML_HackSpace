
#pragma once

namespace GAI
{

	template < typename T > inline T Random(T lMinVal, T lMaxVal) 
	{
		Int64 lA = rand();
		Int64 lB = rand();
		return static_cast<T>((lA | (lB << 32)) % (lMaxVal - lMinVal + 1) + lMinVal);
	}

	template < > inline Float32 Random(Float32 lMinVal, Float32 lMaxVal) { return (static_cast<Float32>(rand()) / RAND_MAX)*(lMaxVal - lMinVal) + lMinVal; }
	template < > inline Float64 Random(Float64 lMinVal, Float64 lMaxVal) { return (static_cast<Float64>(rand()) / RAND_MAX)*(lMaxVal - lMinVal) + lMinVal; }
	
	template < typename T > inline T Random() { return Random(0, (std::numeric_limits<T>::max)()); }
	template < > inline Float32 Random() { return Random(0.0f, 1.0f); }
	template < > inline Float64 Random() { return Random(0.0, 1.0); }
	
	template < typename T > inline T ZeroedRandom();
	template < > inline Int8 ZeroedRandom() { return Random(std::numeric_limits<Int8>::lowest(), std::numeric_limits<Int8>::max()); }
	template < > inline Int16 ZeroedRandom() { return Random(std::numeric_limits<Int16>::lowest(), std::numeric_limits<Int16>::max()); }
	template < > inline Int32 ZeroedRandom() { return Random(std::numeric_limits<Int32>::lowest(), std::numeric_limits<Int32>::max()); }
	template < > inline Int64 ZeroedRandom() { return Random(std::numeric_limits<Int64>::lowest(), std::numeric_limits<Int64>::max()); }
	template < > inline Float32 ZeroedRandom() { return Random(-1.0f, 1.0f); }
	template < > inline Float64 ZeroedRandom() { return Random(-1.0, 1.0); }
	
	void RandomSeed(UInt32 lSeed);
}
