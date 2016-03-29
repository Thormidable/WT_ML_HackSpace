
#pragma once

namespace GAI
{
	namespace NeuronFunctions
	{
		template<class T> class Sigmoid
		{
		public:
			static inline void Function(T &lValue)
			{
				lValue = 1 / (1 + std::exp(-lValue));
			}
			static inline void Prime(T &lValue)
			{
				lValue = exp(-lValue) / std::pow(1 + exp(-lValue),2);
			}
		};

		template<class T> class Linear
		{
		public:
			static inline void Function(T &lValue)
			{
				(void)lValue;
			}
			static inline void Prime(T &lValue)
			{
				lValue = 1;
			}
		};

		template<class T> class TanH
		{
		public:
			static void Function(T& lValue)
			{				
				lValue = std::tanh(lValue);
			}
			static void Prime(T &lValue)
			{
				return std::atanh(lValue);
			}
		};
	}
}
