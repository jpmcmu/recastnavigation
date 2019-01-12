/**
@defgroup detour Detour

Members in this module are wrappers around the standard math library
*/

#ifndef DETOURMATH_H
#define DETOURMATH_H

#include <math.h>
#include "../../../../include/Common/cFixedPoint.h"

#define DT_FLT_MAX 741455
#define DT_FLT_EPSILON 1e-5 
typedef fp64 dtFloat;

inline dtFloat dtMathFabsf(dtFloat x) { return fabs(x); }
inline dtFloat dtMathSqrtf(dtFloat x) { return sqrt(x); }
inline dtFloat dtMathFloorf(dtFloat x) { return floor(x); }
inline dtFloat dtMathCeilf(dtFloat x) { return ceil(x); }
inline dtFloat dtMathCosf(dtFloat x) { return cos(x); }
inline dtFloat dtMathSinf(dtFloat x) { return sin(x); }
inline dtFloat dtMathAtan2f(dtFloat y, dtFloat x) { return atan2(y, x); }

#endif
