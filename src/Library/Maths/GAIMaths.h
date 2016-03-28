
#pragma once

namespace GAI
{
	namespace Maths
	{
		template<class T, class Lambda> T CalculateGradient(T lX, Lambda lFunc, T lEpsilon = std::numeric_limits<T>::epsilon())
		{
			return (lFunc(lX + lEpsilon) - lFunc(lX - lEpsilon)) / (2 * lEpsilon);
		}
	}
 }
	