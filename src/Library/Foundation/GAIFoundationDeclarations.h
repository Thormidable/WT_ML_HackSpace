#pragma once


#define GAI_WIDEN2(str) L ## str
#define GAI_WIDEN(str)  GAI_WIDEN2(str)

#define __WFILE__   GAI_WIDEN(__FILE__)
