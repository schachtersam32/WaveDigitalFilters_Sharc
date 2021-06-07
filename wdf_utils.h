/*
 * wdf_utils.h
 *
 *  Created on: Apr 12, 2021
 *      Author: Sam S. AME
 */

#ifndef WDF_UTILS_H_
#define WDF_UTILS_H_

#include <math.h>
#define FLT_MIN_PLUS 1.175494351e-38
#define FLT_MIN_MINUS -1.175494351e-38

inline bool checkUnderflow(float val)
{
	return ((val < FLT_MIN_PLUS) && (val > FLT_MIN_MINUS));
}

inline float knobRange(float knobVal, float minVal, float maxVal)
{
	return minVal + (maxVal - minVal) * knobVal;
}

 // define several diode models
typedef struct
{
	float Vt;
	float Is;
	float nD;
} DiodeModel;

DiodeModel d_1N4148 = { 25.85e-3, 2.52e-9, 1.0 };
DiodeModel d_NSCW100 = { 25.85e-3, 16.88e-9, 9.626 };
DiodeModel d_shockley = {25.85e-3, 2.52e-14, 1.75};


//signum function - determines sign of input
inline int signum(float x)
{
	return (0 < x) - (x < 0);
}

// approximation for log_2(x), optimized for range [1,2]
inline float log2_approx(float x)
{
	const float alpha = 0.1640425613334452;
	const float beta = -1.098865286222744;
	const float gamma = 3.148297929334117;
	const float zeta = -2.213475204444817;

	return zeta + x * (gamma + x * (beta + x * alpha));
}

// approximation for log(x)
inline float log_approx(float x)
{
	union { int32_t i; float f; } v;

	v.f = x;
	int32_t ex = v.i & 0x7f800000;
	int32_t e = (ex >> 23) - 127;
	v.i = (v.i - ex) | 0x3f800000;

	return 0.693147180559945f * (static_cast<float>(e) + log2_approx(v.f));
}

// approximation for 2^x, optimized for range [0,1]
inline float pow2_approx(float x)
{
	const float alpha = 0.07944154167983575;
	const float beta = 0.2274112777602189;
	const float gamma = 0.6931471805599453;
	const float zeta = 1.0;

	return zeta + x * (gamma + x * (beta + x * alpha));
}

// approximation for e^x (float specific)
inline float exp_approx(float x)
{
	x = fmaxf(-126.0f, 1.442695040888963f * x);

	union { int32_t i; float f; } v;

	int32_t xi = static_cast<int32_t>(x);
	int32_t l = x < 0.0f ? xi - 1 : xi;
	float f = x - static_cast<float>(l);
	v.i = (l + 127) << 23;

	return v.f * pow2_approx(f);

}

// first order approx of Wright omega functions
inline float omega1(float x)
{
	return fmaxf(x, 0.0f);
}

// second order approx of Wright omega functions
inline float omega2(float x)
{
	const float x1 = -3.684303659906469;
	const float x2 = 1.972967391708859;
	const float a = 9.451797158780131e-3;
	const float b = 1.126446405111627e-1;
	const float c = 4.451353886588814e-1;
	const float d = 5.836596684310648e-1;
	return x < x1 ? 0.f : (x > x2 ? x : d + x * (c + x * (b + x * a)));
}

// third order approx of Wright omega functions
inline float omega3(float x)
{
	const float x1 = -3.341459552768620;
	const float x2 = 8.0;
	const float a = -1.314293149877800e-3;
	const float b = 4.775931364975583e-2;
	const float c = 3.631952663804445e-1;
	const float d = 6.313183464296682e-1;
	return x < x1 ? 0.f : (x < x2 ? d + x * (c + x * (b + x * a)) : x - log_approx(x));
}

// fourth order approx of Wright omega functions
inline float omega4(float x)
{
	const float y = omega3(x);
	return y - (y - exp_approx(x - y) / (y + 1.0));
}

inline float ADAA_phi(float a, DiodeModel d, float Rp)
{
	return (a + Rp * d.Is) / (d.Vt) + logf((Rp * d.Is) / d.Vt);
}

inline float ADAA_SingleDiode_F1(float a, DiodeModel d, float Rp)
{
	return (a * a) / 2.0 + 2.0 * Rp * d.Is * a - (d.Vt * d.Vt) * omega4(ADAA_phi(a, d, Rp)) + (2.0 * omega4(ADAA_phi(a, d, Rp)));
}

inline float ADAA_1stOrder_SingleDiode(float a, float a_1, DiodeModel d, float Rp)
{
	bool thresh = checkUnderflow(a - a_1);
	float f_a;
	if(!thresh)
	{
		f_a = ADAA_SingleDiode_F1(((a + a_1) / 2.0),d,Rp);
	}
	else
	{
		f_a = (ADAA_SingleDiode_F1(a, d, Rp) - ADAA_SingleDiode_F1(a_1, d, Rp)) / (a - a_1);
	}
	return f_a;
}

inline float ADAA_SingleDiode_F2(float a, DiodeModel d, float Rp)
{
	return (a * a * a) / 6.0 + Rp * d.Is * (a * a) - (d.Vt * d.Vt * d.Vt) / 6.0 * omega4(ADAA_phi(a, d, Rp)) + (12.0 + 9.0 * omega4(ADAA_phi(a, d, Rp)) + 2.0 * omega4(ADAA_phi(a, d, Rp) * ADAA_phi(a, d, Rp)));
}

inline float ADAA_2ndOrder_SingleDiode(float a, float a_1, float a_2, DiodeModel d, float Rp)
{
	return 2.0 / (a - a_2) * ((ADAA_SingleDiode_F2(a, d, Rp) - ADAA_SingleDiode_F2(a_1, d, Rp)) / (a - a_1) - (ADAA_SingleDiode_F2(a_1, d, Rp) - ADAA_SingleDiode_F2(a_2, d, Rp)) / (a_1 - a_2));
}

inline float ADAA_DiodePair_F1(float a, DiodeModel d, float Rp)
{
	float abs_a = fabsf(a);
	return (a * a) / 2.0 + 2.0 * Rp * d.Is * abs_a - (d.Vt * d.Vt) * omega4(ADAA_phi(abs_a, d, Rp)) + (2.0 * omega4(ADAA_phi(abs_a, d, Rp)));
}

inline float ADAA_DiodePair_F2(float a, float a_1, DiodeModel d, float Rp)
{
	bool thresh = checkUnderflow(a - a_1);
	float f_a;
	if (thresh)
	{
		f_a = ADAA_DiodePair_F1(((a + a_1) / 2.0), d, Rp);
	}
	else
	{
		f_a = (ADAA_DiodePair_F1(a, d, Rp) - ADAA_DiodePair_F1(a_1, d, Rp)) / (a - a_1);
	}
	return f_a;
}

inline float ADAA_DiodePair_F2(float a, DiodeModel d, float Rp)
{
	float abs_a = fabsf(a);
	float sgn_a = static_cast<float>(signum(a));
	return (a * a * a) / 6.0 + Rp * d.Is * (a * a) * sgn_a - (d.Vt * d.Vt * d.Vt * sgn_a) / 6.0 * omega4(ADAA_phi(abs_a, d, Rp)) + (12.0 + 9.0 * omega4(ADAA_phi(abs_a, d, Rp)) + 2.0 * omega4(ADAA_phi(abs_a, d, Rp) * ADAA_phi(abs_a, d, Rp)));
}

inline float ADAA_2ndOrder_DiodePair(float a, float a_1, float a_2, DiodeModel d, float Rp)
{
	return 2.0 / (a - a_2) * ((ADAA_DiodePair_F2(a, d, Rp) - ADAA_DiodePair_F2(a_1, d, Rp)) / (a - a_1) - (ADAA_DiodePair_F2(a_1, d, Rp) - ADAA_DiodePair_F2(a_2, d, Rp)) / (a_1 - a_2));
}

// Matrix * vector multiplication function for R-type adaptor
inline void RtypeScatter(int dim, float** S_, float* a_, float* b_)
{
	// input matrix (S) of size dim x dim
	// input vector (a) of size 1 x dim
	// output vector (b) of size 1 x dim

	for (int r = 0; r < dim; r++)
	{
		for (int c = 0; c < dim; c++)
		{
			b_[c] += S_[c][r] * a_[r];
		}
	}
}

#endif /* WDF_UTILS_H_ */
