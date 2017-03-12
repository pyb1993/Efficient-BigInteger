#include "BigInteger.h"
#include "FFTFunctor.h"
#include "vector"
using namespace std;

class FFTMultiPlier
{
public:
	static const int FFTWidth = 8;
	static const int FFTBase  = 1 << FFTWidth;//
	static const int FFTMax   = FFTBase - 1;// 0xff
	typedef double Real;
	typedef vector<Real> Vectype;
	typedef Vectype&  RefVec;


	FFTMultiPlier(const BigInteger& lhs, const BigInteger& rhs);
	~FFTMultiPlier();
	BigInteger MulWithFFT();


private:
	void Convolution(const RefVec A, const RefVec B, RefVec C);
	void Convert_Base(RefVec Ret, const BigInteger&);
	void UpdateResult();
	//下面个元素是程序热点访问的对象,所以考虑使用堆分配空间
	Vectype ACoef;
	Vectype BCoef;
	Vectype ArrayFFT0;
	Vectype ArrayFFT1;

	Vectype CovResult;
	Vectype Result;

	long ASize;
	long BSize;
	long FFTLength;
	long ResultSize;

};