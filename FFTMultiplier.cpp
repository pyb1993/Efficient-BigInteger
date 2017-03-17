#include "FFTMultiplier.h"
#include "cassert"

#define ASSIGN(_array,x,start,num) do{for(int i =start;i< num;++i) _array[i] = x;}while(0)
#define COPY(_array,_assigner,num) ASSIGN(_array,_assigner[i],0,num)

FFTMultiPlier::FFTMultiPlier(const BigInteger& lhs, const BigInteger& rhs){

	assert(lhs.WIDTH == rhs.WIDTH);
	//处理lhs和rhs获得BASE=10000的结果
	if (lhs.WIDTH != 4){
		//ASize = lhs.length() << 1;
		//BSize = rhs.length() << 1;
		ACoef = Vectype();
		BCoef = Vectype();
		Convert_Base(ACoef, lhs);//给ACoef分配内存并且转换
		Convert_Base(BCoef, rhs);
	}
	ASize = ACoef.size();
	BSize = BCoef.size();
	auto len = ASize + BSize;
	FFTLength = 2;
	while (FFTLength < len)
		FFTLength <<= 1;

	
	ArrayFFT0 = Vectype(2 * FFTLength);//存储lhs的FFT结果
	ArrayFFT1 = Vectype(2 * FFTLength);//存储rhs的FFT结果
	Result = Vectype(len,0);
	}

void FFTMultiPlier::Convolution(const RefVec A, const RefVec B, RefVec C)
{
	long i;
	for (i = 0; i<FFTLength; i++) {
		C[2 * i] = A[2 * i] * B[2 * i] - A[2 * i + 1] * B[2 * i + 1];//tmpR
		C[2 * i + 1] = A[2 * i] * B[2 * i + 1] + A[2 * i + 1] * B[2 * i];//tmpI
	}
}
BigInteger FFTMultiPlier::MulWithFFT(){
	FFTFunctor F(FFTLength);
	F.fft(ACoef, ArrayFFT0, 1);
	F.fft(BCoef, ArrayFFT1, 1);
	vector<Real> CovResult(FFTLength*2,0);
	Convolution(ArrayFFT0, ArrayFFT1,CovResult);//转换结果
	F.InverseFFT(CovResult, Result);
	UpdateResult( );
	return BigInteger(Result,ResultSize);
}
void FFTMultiPlier::Convert_Base( RefVec Ret, const BigInteger& lhs)
{	
	auto len = lhs.length();
	assert(len > 0);
	Ret.reserve(lhs.length() << 2);
	if (lhs.WIDTH != FFTWidth){
		auto s = lhs.s;
		for (size_t i = 0; i < len - 1; ++i){
			auto val = (*s)[i];
			for (size_t j = 0; j < lhs.WIDTH / FFTWidth; ++j)
			{
				Ret.push_back(val & FFTMax);
				val >>= FFTWidth;
			}
		}

		auto remain = (*s)[len - 1];
		while (remain){
			Ret.push_back(remain & FFTMax);
			remain >>= FFTWidth;
		}
	}
}
//处理进位
void FFTMultiPlier::UpdateResult()
{
	assert(Result.size() > 0);
	//合并
	long  carry = 0;
	for (size_t j = 0; j < Result.size(); ++j){
		long tmp = Result[j] + carry;
		carry = tmp >> FFTWidth;
		Result[j] = tmp & FFTMax;
	}
	int i = Result.size() - 1;
	while (i >= 0 && Result[i] == 0) { --i; }
	ResultSize = i>=0 ? i + 1:1 ;//去掉0后的size
}
FFTMultiPlier::~FFTMultiPlier(){

}