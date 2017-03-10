#ifndef FFT_h
#define FFT_h
#include "BigInteger.h"

class FFTFunctor{
public:
	typedef double Real;
	typedef vector<Real> Vectype;
	typedef Vectype&  RefVec;

	const double PI = 3.1415926535897932384626;

	FFTFunctor(size_t fft_length);
	~FFTFunctor();
	 
	void fft(const RefVec Coef, const RefVec FFT, short sign = 1);
	void InverseFFT(const RefVec FFT, RefVec&);
	unsigned int reverseBits(unsigned int num, int count);

	size_t  FFTLength;//FFTµÄ³¤¶È
};	

#endif