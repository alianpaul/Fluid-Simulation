#ifndef UTILITY_H
#define UTILITY_H

#define MAX2(x, y)	(((x) < (y))? (y) : (x))
#define MIN2(x, y)	(((x) < (y))? (x) : (y))
#define POW2(x)		((x) * (x))
#define POW3(x)		((x) * (x) *(x))
#define CLAMP(x, min, max) MIN2( MAX2(x, min), max)

static const float	PI		= 3.1415926535897932384626433832795f;
static const float	TWO_PI	= 2.0f * PI;
static const int	THREAD_NUM = 8;

inline unsigned int NearestPowerOfTwoExponent(unsigned int iVal)
{
	if (iVal == 0)
		return 0;

	--iVal;

	int shift = 1;
	while ((iVal >> shift) > 0) 
	{
		++shift;
	}

	return shift;
}

inline unsigned int NearestPowerOfTwo(unsigned int iVal)
{
	return 1 << NearestPowerOfTwoExponent(iVal);
}

inline float RandomSpread( float fSpread )
{
	return fSpread * ( float(rand()) / float(RAND_MAX) -0.5f);
}

inline float sechf(const float & x) { return 1.0f / coshf(x); }

inline float finvsqrtf(const float & val)
{
	long    i = (long&)val;             // Exploit IEEE 754 inner workings.
	i = 0x5f3759df - (i >> 1);          // From Taylor's theorem and IEEE 754 format.
	float   y = (float&)i;              // Estimate of 1/sqrt(val) close enough for convergence using Newton's method.
	static const float  f = 1.5f;        // Derived from Newton's method.
	const float         x = val * 0.5f;  // Derived from Newton's method.
	y = y * (f - (x * y * y));        // Newton's method for 1/sqrt(val)
	y = y * (f - (x * y * y));        // Another iteration of Newton's method
	return y;
}

#endif