
#pragma once

//This tells the system not to use DLL interfaces for the API.
//#define GAI_STATIC_LINK 

// Linux
#if defined(__linux)
#	define GAI_HOST_LINUX 1
#	define GAI_PLATFORM_32BIT 1
#	pragma warning Distinction between 32bit and 64bit linux not yet implemented. Defaulting to 32bit.
#	define GAI_GCC 1
#endif // defined(__linux)

// Windows
#if defined(_MSC_VER)
#	define GAI_HOST_WINDOWS 1
#	if defined(__CUDACC__)
#		define GAI_PLATFORM_CUDA 1
#	elif defined(_WIN64)
#		define GAI_PLATFORM_64BIT 1
#	elif defined(_WIN32)
#		define GAI_PLATFORM_32BIT 1
#	else
#		error Platform not supported or not detected correctly
#	endif
#	define GAI_MSVC 1
#endif // defined(_MSC_VER)

// Catch errors
#if !defined(GAI_HOST_WINDOWS) && !defined(GAI_HOST_LINUX)
#	error Host operating system not supported or not detected correctly
#endif

// CUDA
#if defined(__NVCC__)
#	define GAI_CUDA 1
#endif

// API definition
#if GAI_HOST_WINDOWS && !defined(GAI_STATIC_LINK)
#	if defined(GAI_API_EXPORT)
#		define GAI_API __declspec(dllexport) 
#		define GAI_INSTANCED_TEMPLATE
#	else
#		define GAI_API __declspec(dllimport) 
#		define GAI_INSTANCED_TEMPLATE extern
#	endif
#else
#	define GAI_API
#	define GAI_INSTANCED_TEMPLATE
#endif // GAI_HOST_WINDOWS


#define GAI_POINTER_SIZE sizeof(void*)

#define GAI_PI 3.14156f
