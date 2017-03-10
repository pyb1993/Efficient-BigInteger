#include "FFTFunctor.h"
#include "cassert"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "iostream"
//FFTFunctor 默认所有的长度是2的幂次,由caller来保证长度和补0

FFTFunctor::FFTFunctor(size_t fft_length) :
				 FFTLength(fft_length)
{
	assert((FFTLength & (FFTLength - 1)) == 0);// 小 trick,确保长度是2的幂次	
}

FFTFunctor::~FFTFunctor(){}

void FFTFunctor::fft(const RefVec Coef, RefVec FFT, short sign){

		//分别处理实数和虚数
		if (sign == -1){
			for (int i = 0; i < Coef.size(); ++i)
				FFT[i] = Coef[i];
		}
		else{
			for (int i = 0; i < FFTLength; ++i)
			{
				FFT[2 * i] = (i < Coef.size()) ? Coef[i] : 0;
				FFT[2 * i + 1] = 0;
			}
	}
	//统计位数数量
	int t = FFTLength -1;
	int count = 0;
	while (t){ 
		t = t &(t - 1);
		++count;
	}
	for (int i = 1; i < FFTLength - 1; ++i){
		auto j = reverseBits(i, count);
		if (j > i)
		{
			swap(FFT[2 * i], FFT[2 * j]);
			swap(FFT[2 * i+1], FFT[2 * j+1]);
		  }
	}

	//interval代表当前FFT的子序列的长度
	for (size_t interval = 1; interval < FFTLength; interval*= 2) // 执行 log2(nn) 次外循环
		       {
				cout << "outer loop\n";
		         double theta = (2*sign) * PI / (2*interval);
		         double alpha = sin(0.5 * theta);
		         alpha = -2 * alpha * alpha;
		         double beta = sin(theta);
		         double OmegaR = 1;
		         double OmegaI = 0;
		         for (int m = 0;  m < interval; ++m)//m代表每一个子变换的起始位置
			         {
					   cout << " middle loop  "<<"("<<OmegaR<<","<<OmegaI<<")\n";
			           for (int i = m; i < FFTLength; i += 2*interval)// FFT[i]代表每一个子序列在当前位置的FFT系数
				           {
				             size_t j = i + interval; // 下面是 Danielson-Lanczos 公式
							 cout << "   innner loop : i = " << i << "  j = " << j << endl;
				             Real tmpR = OmegaR * FFT[2*j] - OmegaI * FFT[2*j + 1];
							 Real tmpI = OmegaR * FFT[2*j + 1] + OmegaI * FFT[2*j];
				             FFT[2*j] = FFT[2*i] - tmpR;
				             FFT[2*j + 1] = FFT[2*i + 1] - tmpI;
				             FFT[2*i] += tmpR;
				             FFT[2*i + 1] += tmpI;
				           }
					   //注意上面为什么是这样顺序的两层循环?原因是方便Wk的计算!
					   auto tmp = OmegaR;
			           OmegaR = OmegaR * alpha - OmegaI *beta  + OmegaR; // 三角递归
			           OmegaI = OmegaI * alpha +  tmp* beta + OmegaI;
			         }
		       }
	}



void FFTFunctor::InverseFFT(const RefVec FFT, RefVec Coef)
{
	long i;
	Real invNFFT = 1. / (Real)FFTLength, tmp;
	Vectype IFFT(2*FFTLength, 0);
	fft(FFT,IFFT,-1);
	auto len = Coef.size();
	for (i = 0; i<2 * len; i+=2) 
	{
		/* 四舍五入,还原序列*/
		tmp = invNFFT*IFFT[i];
		Coef[i/2] = (floor(0.5 + tmp));
	}

}



unsigned int FFTFunctor::reverseBits(unsigned int num, int count)
{
	unsigned int reverse_num = 0;

	//num >>= 1;
	while (count--)
	{
		reverse_num <<= 1;
		reverse_num |= num & 1;
		num >>= 1;
	}
	return reverse_num;
}
