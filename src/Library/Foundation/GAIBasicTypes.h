
#pragma once

#define GAI_ALIGN_VALUE(Value, Alignment) ((Value) & ~((Alignment) - 1))

#define GAI_EMPTY(...)
#define GAI_DEFER(S) S GAI_EMPTY()

// Feel free to extend this macro to support additional expansions
#define GAI_EXPAND_INNER0(...) __VA_ARGS__
#define GAI_EXPAND_INNER1(...) GAI_EXPAND_INNER0(__VA_ARGS__)
#define GAI_EXPAND_INNER2(...) GAI_EXPAND_INNER1(__VA_ARGS__)
#define GAI_EXPAND_INNER3(...) GAI_EXPAND_INNER2(__VA_ARGS__)
#define GAI_EXPAND_INNER4(...) GAI_EXPAND_INNER3(__VA_ARGS__)
#define GAI_EXPAND(...) GAI_EXPAND_INNER4(__VA_ARGS__)

#define GAI_FOR_EACH_UINT_TYPE(M, ...) \
	GAI_DEFER(M)(GAI::UInt8, __VA_ARGS__) \
	GAI_DEFER(M)(GAI::UInt16, __VA_ARGS__) \
	GAI_DEFER(M)(GAI::UInt32, __VA_ARGS__) \
	GAI_DEFER(M)(GAI::UInt64, __VA_ARGS__)

#define GAI_FOR_EACH_INT_TYPE(M, ...) \
	GAI_DEFER(M)(GAI::Int8, __VA_ARGS__) \
	GAI_DEFER(M)(GAI::Int16, __VA_ARGS__) \
	GAI_DEFER(M)(GAI::Int32, __VA_ARGS__) \
	GAI_DEFER(M)(GAI::Int64, __VA_ARGS__)

#define GAI_FOR_EACH_INTEGER_TYPE(M, ...) \
	GAI_FOR_EACH_UINT_TYPE(M, __VA_ARGS__) \
	GAI_FOR_EACH_INT_TYPE(M, __VA_ARGS__)

#define GAI_FOR_EACH_FLOAT_TYPE(M, ...) \
	GAI_DEFER(M)(GAI::Float32, __VA_ARGS__) \
	GAI_DEFER(M)(GAI::Float64, __VA_ARGS__)

#define GAI_FOR_EACH_BASIC_DATA_TYPE(M, ...) \
	GAI_FOR_EACH_INTEGER_TYPE(M, __VA_ARGS__) \
	GAI_FOR_EACH_FLOAT_TYPE(M, __VA_ARGS__)

#define GAI_FOR_COUNT_0(...)
#define GAI_FOR_COUNT_1(M, ...) GAI_DEFER(M)(0, __VA_ARGS__)
#define GAI_FOR_COUNT_2(M, ...) GAI_FOR_COUNT_1(M, __VA_ARGS__) GAI_DEFER(M)(1, __VA_ARGS__)
#define GAI_FOR_COUNT_3(M, ...) GAI_FOR_COUNT_2(M, __VA_ARGS__) GAI_DEFER(M)(2, __VA_ARGS__)
#define GAI_FOR_COUNT_4(M, ...) GAI_FOR_COUNT_3(M, __VA_ARGS__) GAI_DEFER(M)(3, __VA_ARGS__)
#define GAI_FOR_COUNT_5(M, ...) GAI_FOR_COUNT_4(M, __VA_ARGS__) GAI_DEFER(M)(4, __VA_ARGS__)
#define GAI_FOR_COUNT_6(M, ...) GAI_FOR_COUNT_5(M, __VA_ARGS__) GAI_DEFER(M)(5, __VA_ARGS__)
#define GAI_FOR_COUNT_7(M, ...) GAI_FOR_COUNT_6(M, __VA_ARGS__) GAI_DEFER(M)(6, __VA_ARGS__)
#define GAI_FOR_COUNT_8(M, ...) GAI_FOR_COUNT_7(M, __VA_ARGS__) GAI_DEFER(M)(7, __VA_ARGS__)
#define GAI_FOR_COUNT_9(M, ...) GAI_FOR_COUNT_8(M, __VA_ARGS__) GAI_DEFER(M)(8, __VA_ARGS__)
#define GAI_FOR_COUNT_10(M, ...) GAI_FOR_COUNT_9(M, __VA_ARGS__) GAI_DEFER(M)(9, __VA_ARGS__)
#define GAI_FOR_COUNT_11(M, ...) GAI_FOR_COUNT_10(M, __VA_ARGS__) GAI_DEFER(M)(10, __VA_ARGS__)
#define GAI_FOR_COUNT_12(M, ...) GAI_FOR_COUNT_11(M, __VA_ARGS__) GAI_DEFER(M)(11, __VA_ARGS__)
#define GAI_FOR_COUNT_13(M, ...) GAI_FOR_COUNT_12(M, __VA_ARGS__) GAI_DEFER(M)(12, __VA_ARGS__)
#define GAI_FOR_COUNT_14(M, ...) GAI_FOR_COUNT_13(M, __VA_ARGS__) GAI_DEFER(M)(13, __VA_ARGS__)
#define GAI_FOR_COUNT_15(M, ...) GAI_FOR_COUNT_14(M, __VA_ARGS__) GAI_DEFER(M)(14, __VA_ARGS__)
#define GAI_FOR_COUNT_16(M, ...) GAI_FOR_COUNT_15(M, __VA_ARGS__) GAI_DEFER(M)(15, __VA_ARGS__)
#define GAI_FOR_COUNT_17(M, ...) GAI_FOR_COUNT_16(M, __VA_ARGS__) GAI_DEFER(M)(16, __VA_ARGS__)
#define GAI_FOR_COUNT_18(M, ...) GAI_FOR_COUNT_17(M, __VA_ARGS__) GAI_DEFER(M)(17, __VA_ARGS__)
#define GAI_FOR_COUNT_19(M, ...) GAI_FOR_COUNT_18(M, __VA_ARGS__) GAI_DEFER(M)(18, __VA_ARGS__)
#define GAI_FOR_COUNT_20(M, ...) GAI_FOR_COUNT_19(M, __VA_ARGS__) GAI_DEFER(M)(19, __VA_ARGS__)
#define GAI_FOR_COUNT_21(M, ...) GAI_FOR_COUNT_20(M, __VA_ARGS__) GAI_DEFER(M)(20, __VA_ARGS__)
#define GAI_FOR_COUNT_22(M, ...) GAI_FOR_COUNT_21(M, __VA_ARGS__) GAI_DEFER(M)(21, __VA_ARGS__)
#define GAI_FOR_COUNT_23(M, ...) GAI_FOR_COUNT_22(M, __VA_ARGS__) GAI_DEFER(M)(22, __VA_ARGS__)
#define GAI_FOR_COUNT_24(M, ...) GAI_FOR_COUNT_23(M, __VA_ARGS__) GAI_DEFER(M)(23, __VA_ARGS__)
#define GAI_FOR_COUNT_25(M, ...) GAI_FOR_COUNT_24(M, __VA_ARGS__) GAI_DEFER(M)(24, __VA_ARGS__)
#define GAI_FOR_COUNT_26(M, ...) GAI_FOR_COUNT_25(M, __VA_ARGS__) GAI_DEFER(M)(25, __VA_ARGS__)
#define GAI_FOR_COUNT_27(M, ...) GAI_FOR_COUNT_26(M, __VA_ARGS__) GAI_DEFER(M)(26, __VA_ARGS__)
#define GAI_FOR_COUNT_28(M, ...) GAI_FOR_COUNT_27(M, __VA_ARGS__) GAI_DEFER(M)(27, __VA_ARGS__)
#define GAI_FOR_COUNT_29(M, ...) GAI_FOR_COUNT_28(M, __VA_ARGS__) GAI_DEFER(M)(28, __VA_ARGS__)
#define GAI_FOR_COUNT_30(M, ...) GAI_FOR_COUNT_29(M, __VA_ARGS__) GAI_DEFER(M)(29, __VA_ARGS__)
#define GAI_FOR_COUNT_31(M, ...) GAI_FOR_COUNT_30(M, __VA_ARGS__) GAI_DEFER(M)(30, __VA_ARGS__)
#define GAI_FOR_COUNT_32(M, ...) GAI_FOR_COUNT_31(M, __VA_ARGS__) GAI_DEFER(M)(31, __VA_ARGS__)

#define GAI_FOR_COUNT(COUNT, M, ...) \
	GAI_FOR_COUNT_##COUNT(M, __VA_ARGS__)

#define GAI_BIT(bitIndex) (1 << (bitIndex))

#define GAI_NOT_FOUND ((std::numeric_limits<GAI::Size>::max)())

/*!
Ibex Public Namespace. Contains all Exposed Ibex API Functionality
*/
namespace GAI
{
	typedef int64_t				Int64;
	typedef int32_t				Int32;
	typedef int16_t				Int16;
	typedef int8_t				Int8;

	typedef uint64_t			UInt64;
	typedef uint32_t			UInt32;
	typedef uint16_t			UInt16;
	typedef uint8_t				UInt8;

	typedef double				Real;
	typedef long				Integer;

	typedef long double			Float80;
	typedef double				Float64;
	typedef float				Float32;

	typedef size_t				Size;
	typedef size_t				Offset;

	typedef wchar_t				WChar;

	typedef uint32_t			Flags;

	typedef bool				Bool;
	typedef char				Char;

	template<typename R>
	struct BoolSafeStaticCast
	{
		template<typename T>
		static inline R Cast(const T& value)
		{
			return static_cast<R>(value);
		}
	};

	template<>
	struct BoolSafeStaticCast<GAI::Bool>
	{
		template<typename T>
		static inline GAI::Bool Cast(const T& value)
		{
			return value != 0;
		}

		template<>
		static inline GAI::Bool Cast<Float32>(const Float32& value)
		{
			return value > std::numeric_limits<Float32>::epsilon();
		}

		template<>
		static inline GAI::Bool Cast<Float64>(const Float64& value)
		{
			return value > std::numeric_limits<Float64>::epsilon();
		}
	};

	template <class A, class B> struct CompareTypes
	{
		static const GAI::Bool result = false;
	};

	template <class A> struct CompareTypes<A, A>
	{
		static const GAI::Bool result = true;
	};

	typedef unsigned short tAPIImageType;
	typedef Float64 tDatabaseDataType;
	
}
