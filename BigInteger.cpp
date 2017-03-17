#include "BigInteger.h"
#include "StringBigdata.h"
#include "sstream"
#include "cassert"
#include "algorithm"
#include "bitset"
#include "FFTMultiplier.h"
#include <utility>

using namespace std;

map<size_t, BigInteger> BigInteger::TenPower = map<size_t, BigInteger>();
map<size_t,map<size_t, BigInteger>> BigInteger::PreRepricoal = map<size_t,map<size_t, BigInteger>> ();

BigInteger::BigInteger()
{
	s = new datatype();
}

BigInteger::BigInteger(const BigInteger& rhs) :\
s(new datatype(*rhs.s)), negative(rhs.negative){   }

//convert data fom FFT coeficients to BigInteger
BigInteger::BigInteger(const vector<double>& Coef, size_t len,unsigned int other_width)
:s( new datatype)
{
	s->reserve(len >> 2);
	maxtype carry = 0;
	size_t i = 0;
	for (; i+3 < len; i += 4){
		maxtype x = carry + (maxtype)Coef[i];// 1
		x += ( (maxtype(Coef[i + 1])) << (other_width));// 8
		x += ((maxtype(Coef[i + 2])) << (2*other_width));//16
		x += ((maxtype(Coef[i + 3])) << (3*other_width));//24

		if (x > BASE)
		{
			carry = x>>WIDTH;
			x &= Max;
		}
		s->push_back(x & Max);
	}
		for (int k = 0; i < len;++i,++k)
			carry += (maxtype)Coef[i]<<(other_width*k) ;

	while (carry) {
		s->push_back(carry&Max);
		carry >>= WIDTH;
	}

	if (length() == 0)//corner case 
		s->push_back(0);

}

//move construction
BigInteger::BigInteger(BigInteger&& rhs) NOEXCEPT  :
negative(rhs.negative), s(rhs.s)
{
	rhs.s = nullptr;
}

BigInteger::BigInteger( unsigned long long& num)
{
	s = new datatype();
	*this = num; 
}

BigInteger::BigInteger(const  long long& num)
{
	s = new datatype();
	*this = num;
}

BigInteger::BigInteger(const string& str,size_t sublen):
substr_len(sublen)
{
	s = new datatype();
	*this = str;
}

BigInteger& BigInteger::operator= (const BigInteger& rhs){
	if (this != &rhs){
		BigInteger  temp (rhs);
	    //swap 资源,自动释放局部变量temp
		swap(*this, temp);
	}
	return *this;
}


BigInteger& BigInteger::operator= ( long long num){
	if (num < 0){
		negative = true;
		num = -num;
	}

	do{
		s->push_back(num&Max);
		num >>= WIDTH;
	} while (num > 0);
	return *this;
}
#if 1
BigInteger& BigInteger:: operator= (unsigned long long& _num){
	auto num = _num;
	do{
		s->push_back(num & Max);
		num >>= WIDTH;
	} while (num > 0);
	return *this;
}
#endif
BigInteger& BigInteger:: operator= (const string& str){
	
	/****Deal with the Binary string : bin,hex******/
	int str_len = str.length();
	assert(str_len != 0);

	size_t bound = 1025;
	Radix radix;
	int beg;

	/*******For short decimal string OR Binary string, directly conversion*******/
	get_type(str, radix, negative, beg);
	if (radix != Dec || str.length() < bound){
		
		if (radix == Dec)
			parse_binary(convert_to_binary(str.begin()+beg,str.length()-beg), Bin);
		else
			parse_binary(beg ? str.substr(beg,-1) : str,radix);

		remove_redundant_zeros();
		return *this;
	}
	
	/***********Deal with the Decimal string***************************/
	size_t i = 0;
	size_t _len = substr_len;
	auto strlen = str.length();

	if (str[0] == '-')// skip sign
		++i;

	BigInteger& result = *this;
	result = BigInteger(str.substr(i, _len),substr_len>>1);
	i += _len;

	// each iteration, we convert _len digits to BigInteger by ret = ret*(10.....0) + remain
	// we store the BigInteger 10000...0 when we compute it for reuse it efficiently
	while (i + _len < strlen)
	{
		if (TenPower.find(_len) == TenPower.end()){
			TenPower[_len] = BigInteger("1" + string(_len, '0'), _len >> 1);//recursively convert, but the _len become half
		}

		result  *= TenPower[_len]; 
		result  += BigInteger(str.substr(i, _len),_len>>1);
		i += _len;
	}

	
	if (i  < strlen)
	{
		auto remain = str.substr(i, -1);
		result *= BigInteger("1" + string(remain.length(), '0'),_len>>1); //* 1000
		result += remain;
	}
		remove_redundant_zeros();
		return *this;
}

BigInteger& BigInteger:: operator= (BigInteger&& rhs)NOEXCEPT{
	if (this != &rhs){
		delete s;
		s = rhs.s;
		negative = rhs.negative;
		rhs.s = nullptr;
	}
	return *this;
}



BigInteger:: ~BigInteger(){
	delete s;
}

void add(const BigInteger& lhs, long num, BigInteger& ret){

	auto Max = BigInteger::Max;
	auto WIDTH = BigInteger::WIDTH;

	unsigned long long  carry = num;
	unsigned long long  tmp = 0;
	auto len = lhs.length();
	auto &ls = *lhs.s;
	ret.resize(len);

	size_t i = 0;
	while (i < len)
	{
		tmp = carry + ls[i];
		ret[i++] = tmp & Max;
		carry = tmp >> WIDTH;
	}
	if (carry)
		ret.push_back(carry);
}

/****consider sign*****/
void _add(const BigInteger& lhs, long num, BigInteger& ret)
{
	if (lhs.negative == false)
	{
		if (num > 0)
			add(lhs, num, ret);
		else
			subtract(lhs, num, ret);
	}

	else
	{
		if (num > 0)
			add(lhs, -num, ret);
		else
			subtract(-num, -lhs, ret);
	}
}


BigInteger BigInteger:: operator + (const long long& num) const
{
	BigInteger ret;
	add(*this, num, ret);
	return ret;
}

BigInteger BigInteger:: operator + (const BigInteger& rhs) const{
	
	BigInteger ret;
	_add(*this, rhs, ret);

	ret.remove_redundant_zeros();
	return ret;
}

void BigInteger::operator += (const BigInteger& rhs)
{
	_add(*this, rhs, *this);
	this->remove_redundant_zeros();
}

/*******more efficient than operator + , avoid copy******************/
void BigInteger::operator += (long num)
{
	_add(*this, num, *this);
}


/*********consider sign***********************************/
void _subtract(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret){
	
	if (lhs.negative == false)
	{
		if (rhs.negative == false)
			subtract(lhs, rhs, ret);
		else
			add(lhs,-rhs,ret);// + - - = + + +
	}
	else
	{
		if (rhs.negative == false)
		{
			add(-lhs, rhs, ret);// - - + = - (+ + +)
			ret.negative = true;
		}
		else
			subtract(-rhs, -lhs, ret);// - - - = + - +
	}
}

/********help function for operator -,use complementary *******************************/
void subtract(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret, bool neg)
{
	typedef BigInteger::maxtype maxtype;
	auto WIDTH =  BigInteger::WIDTH;

	auto llen = lhs.length();
	auto rlen = rhs.length();

	if (lhs < rhs){
		subtract(rhs, lhs, ret, true);
		ret.set_resversable();
		return;
	}

	auto Max = BigInteger::Max;
	ret.resize(llen);

	maxtype subtracted;
	maxtype carry = 1;
	auto &ls = *lhs.s;
	auto &rs = *rhs.s;

	size_t i = 0;
	while (i < rlen )
	{
		subtracted = carry + ls[i] + (~rs[i]);
		carry = subtracted >> WIDTH;
		ret[i++] = subtracted & Max;
	}
	while (i < llen - 1)
	{
		subtracted = ls[i] + (Max) + carry;
		carry = subtracted >> WIDTH;
		ret[i++] = subtracted & Max;
	}
	if (i < llen){
		subtracted = lhs[i] + Max + carry;
		carry = subtracted >> WIDTH;
		ret[i++] = (subtracted & Max) - (carry == 0);
	}

}


BigInteger BigInteger::operator - (const BigInteger& rhs) const
{
	BigInteger ret;
	_subtract(*this, rhs,ret);
	ret.remove_redundant_zeros();
	return ret;
}

/*******more efficient than operator - , avoid copy******************/
void BigInteger::operator -= (const BigInteger& rhs)
{

	_subtract(*this, rhs, *this);
	this->remove_redundant_zeros();
}

BigInteger BigInteger:: operator *(const long & num) const{
	BigInteger ret;
	mul(*this, num, ret);
	ret.negative = (negative != (num < 0));
	return ret;
}

/************help function, consider sign****************************/
void _add(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret)
{
	if (lhs.negative == false)
	{
		if (rhs.negative == false)
			add(lhs, rhs, ret);
		else
			subtract(lhs, -rhs, ret);
	}

	else
	{
		if (rhs.negative == false)
			subtract(rhs, lhs, ret);
		else
		{
			add(lhs, rhs, ret);
			ret.negative = true;
		}
	}
}
/*********help function for operator +********************
which is a bit heavy for the consideration of performance and effiency */
void add(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret){

	typedef BigInteger::maxtype maxtype;
	auto Max = BigInteger::Max;
	auto WIDTH = BigInteger::WIDTH;

	size_t llen = lhs.length(), rlen = rhs.length();
	if (llen < rlen)
	{
		add(rhs, lhs, ret);
		return;
	}

	maxtype x = 0;
	maxtype carry = 0;

	ret.resize(llen);
	auto &ls = *lhs.s;
	auto &rs = *rhs.s;

	int i = 0;
	while (i < rlen){
		x = carry + ls[i] + rs[i];
		carry = x >> WIDTH;
		ret[i++] = x & Max;
	}
	
	while (i < llen && carry){
		x = carry + ls[i];
		carry = x >> WIDTH;
		ret[i++] = x & Max;
	}
	while (i < llen){
		ret[i] = ls[i];
		i++;
	}

	if (carry)
		ret.push_back(carry);
}

/***********help function for small multiplication*********************/
void mul(const BigInteger& lhs, long num, BigInteger& ret)
 {
	auto WIDTH = BigInteger::WIDTH;
	auto Max   = BigInteger::Max;
	typedef BigInteger::eletype eletype;

	auto len = lhs.length();
	auto &ls = *(lhs.s);

	ret.resize(len);
	unsigned long long carry = 0;
	unsigned long long tmp = 0;
	size_t i = 0;
	for(auto& val:ls){
		tmp = carry + (unsigned long long)num * val;//不会溢出
		carry = tmp >> WIDTH;
		ret[i++] = (tmp & Max);
	}
	while (carry){
		ret.push_back(carry);
		carry >>= WIDTH;
	}
 }


/******help function for multiplication****************************/
/******the main framework is Karatsuba which is used to split BigInteger into small pieces 
	   the non-recursive FFT Base multiplication is used to compute BigInteger smaller than 65536 digits **********************************/

void mul(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret) 
{
	//classical KR aloorithm
	typedef BigInteger::maxtype maxtype;
	size_t	WIDTH = BigInteger::WIDTH;

	auto llen = lhs.length(), rlen = rhs.length();
	assert(llen > 0 && rlen > 0);
	if (llen < 65536 && rlen < 65536)//avoid error from FFT
	{
		auto F = FFTMultiPlier(lhs, rhs);
		ret = F.MulWithFFT();
		return;
	}
	//确保a大于b
	if (llen < rlen)
	{
		mul(rhs, lhs, ret);
		return;
	}

	if (lhs == 0 || rhs == 0)
		ret  = BigInteger("0");

	if (rlen * 2 > llen){
		size_t sub_len = llen >> 1;   //低位的长度,确保低位长度一致

		size_t mid1 = llen - sub_len, mid2 = rlen - sub_len;//注意这里的trick
		auto a1 = lhs.sub_integer_rough(0, mid1), a0 = lhs.sub_integer_rough(mid1);
		auto b1 = rhs.sub_integer_rough(0, mid2), b0 = rhs.sub_integer_rough(mid2);
		size_t weight = sub_len * WIDTH;//权重 = 32*sublen

		auto middle_res = (a1 + a0)*(b1 + b0);//(a1+a0)(b1+b0)

		a1 *= b1;// z2 = a1*b1
		a0 *= b0;// z0 = a0*b0


		middle_res -= a1 - a0;

		middle_res <<= weight;
		a1 <<= (weight << 1); //  z2 * b^xx
		a1 += middle_res;
		a1 += a0;// + z0
		ret = a1;
	}
	else{
		//b的位数太小,不足以分割
		size_t mid = llen >> 1;
		auto a1 = lhs.sub_integer_rough(0, mid), a0 = lhs.sub_integer_rough(mid);
		size_t weight = a0.length() * WIDTH;//权重 = 32*低位长度

		a1 *= rhs;
		a0 *= rhs;
		a1 <<= weight;
		a1 += a0;
		ret = a1;
	}
	ret.remove_redundant_zeros();
}

BigInteger BigInteger::operator * (const BigInteger& rhs) const
{
	BigInteger ret;
	ret.reserve(length() + rhs.length());
	mul(*this, rhs, ret);
	ret.negative = (negative != rhs.negative);
	return ret;
}

/*****************power function  **********************************/
BigInteger BigInteger::pow(const long long& rhs) const
 {
	if (*this == 0) return 0;
	auto highest_bit = (static_cast<maxtype>(1) << 63);
	while ((highest_bit & rhs) == 0)
		highest_bit >>= 1;

	BigInteger ret = *this;
	while (highest_bit > 1){
	
		highest_bit >>= 1;
		ret *= ret;
		if (rhs & highest_bit){
			ret *= *this;
		}
	}
	return ret;
 }

void BigInteger:: operator *= (long num)
{
	mul(*this, num, *this);
	negative = (negative != (num < 0));
}

void BigInteger:: operator *= (const BigInteger& rhs)
{
	mul(*this, rhs, *this);
	negative = (negative != rhs.negative);
}
/*****************help function for shift left, no mater whether ret and lhs is same  **********************************/
void shift_left(const BigInteger& lhs, BigInteger& ret, size_t shift_num){
	
	if (lhs == 0)
	{
		ret = BigInteger(0);
		return;
	}

	auto WIDTH = BigInteger::WIDTH;
	auto Max = BigInteger::Max;
	typedef BigInteger::maxtype maxtype;

	unsigned int shift_length = shift_num >> 5;//move long distance
	unsigned int shift_remain = shift_num - shift_length*WIDTH;//move remain distance
	ret.negative = lhs.negative;
	auto len = lhs.length();
	ret.resize(shift_length + len);
	auto &ls = *(lhs.s);

	if (shift_length)
	{
		for (size_t j = len; j--;)
			ret[j + shift_length] = ls[j];
	}

	maxtype carry = 0;
	size_t i = shift_length;
	for (; i < len + shift_length; ++i)
	{
		maxtype  val = (maxtype(ls[i - shift_length]) << shift_remain) | carry;//long long转换,为了多移位
		ret[i] = val & Max;
		carry = val >> WIDTH;
	}
	if (carry > 0)
		ret.push_back(carry);

}
BigInteger BigInteger:: operator <<(size_t shift_num) const
{
	BigInteger ret;
	shift_left(*this, ret, shift_num);
	return ret;
}

void BigInteger:: operator <<= (size_t shift_num) 
{
	shift_left(*this, *this, shift_num);
}


/*******help function for shift right >> *******************************/
//corner case: we need to consider ret is lhs itself!
void shift_right(const BigInteger& lhs, BigInteger& ret, size_t shift_num){

	if (lhs == 0){
		ret =  BigInteger(0);
		return ;
	}

	auto WIDTH = BigInteger::WIDTH;
	auto Max = BigInteger::Max;
	typedef BigInteger::maxtype maxtype;

	size_t shift_length = shift_num >> 5;
	size_t shift_remain = shift_num - shift_length*WIDTH;

	auto &ls = *lhs.s;
	ret.negative = lhs.negative;
	auto len = lhs.length();

	if (shift_length >= len){
		ret = BigInteger(0);
		return;
	}

	if (ret.length() == 0)
		ret.resize(len - shift_length);
	
	if (shift_length)
	{
		for (size_t j = shift_length; j < len; ++j)//整体移动
			ret[j - shift_length] = ls[j];
	}


	unsigned int mask = (1 << shift_remain) - 1;
	unsigned int invshift = WIDTH - shift_remain;
	maxtype carry = 0;

	if (ret.length() != len - shift_length)
	{
		auto s = ret.s;
		s->erase(s->begin() + len - shift_length, s->end());
		for (size_t i = ret.length(); i--;)
		{
			maxtype val = (ret[i] >> shift_remain) | carry;
			carry = (lhs[i] & mask) << invshift;
			ret[i] = val;
		}
	}
	else{
		for (size_t i = len - shift_length; i--;)
		{
			maxtype val = (lhs[i + shift_length] >> shift_remain) | carry;
			carry = (lhs[i + shift_length] & mask) << invshift;
			ret[i] = val;
		}
	}
	ret.remove_redundant_zeros();
}

BigInteger BigInteger:: operator >>(size_t shift_num) const
{
	BigInteger ret;
	shift_right(*this, ret, shift_num);
	return ret;
}

void BigInteger:: operator >>=(size_t shift_num) 
{
	shift_right(*this, *this, shift_num);
}


BigInteger BigInteger::operator / (const BigInteger& rhs) const
{
		BigInteger q, ret;
		div(*this,rhs,q,ret);
		q.negative = (negative != rhs.negative);
		return q;
}

/********help function to compute division by modified newton-raphason algorithm***********************/
void div(const BigInteger& lhs, const BigInteger& rhs, BigInteger& quotient, BigInteger& ret,bool byconst)
{
	if ( rhs == 0 )
		throw overflow_error("Divide by zero exception");
	if (lhs < rhs)
	{
		quotient = BigInteger(0);
		ret = lhs;
		return;
	}

	const BigInteger& dividend = lhs;
	const BigInteger& divisor  = rhs;

	size_t k = dividend.numBits() + divisor.numBits();
	BigInteger x;
	
	if (byconst == false || (x = find_appropriate_rep(k,divisor.length())) == 0)// not find
	{
		BigInteger pow(1);
		auto pow2 = pow << (k + 1);
		x = pow << (k - divisor.numBits());
		BigInteger lastx[] = { 0, 0 };
		bool c = 0;
		int count = 0;

		while (1)
		{
			if (k < 1000000 || count++ > 6)
			{
				x *= pow2 - x * divisor;
				x >>= k;
			}

			else
				x = x.truncate_mul(pow2 - x.truncate_mul(divisor)) >> k;

			if (x == lastx[c] || x == lastx[!c]) break;
			lastx[c] = x;
			c = !c;
		}
		if (byconst)
			BigInteger::PreRepricoal[divisor.length()][k] = x;
	}

	quotient = (dividend * x) >> k;
	ret = dividend - ( quotient * divisor );
	while ( ret >= divisor )
	{
		quotient += 1;
		ret -= divisor;
	}
}

BigInteger BigInteger::operator / (const long num) const
{
	auto len = length();
	assert(len > 0 && num != 0);
	BigInteger ret;
	ret.reserve(len);
	const auto& lhs = *s;
	maxtype remain = 0;
	for (size_t i = len; i--;){
		
		maxtype val = ((remain << WIDTH) + lhs[i]);
		maxtype x = val / num;
		remain = val % num;
		ret.push_back(x);
	}

	reverse(ret.s->begin(), ret.s->end());
	ret.remove_redundant_zeros();
	return ret;
}



BigInteger BigInteger::operator % (const BigInteger& rhs) const
{

	BigInteger q, ret;
	div(*this,rhs,q,ret);
	
	if (negative && !rhs.negative)// -7 % -5 = -2
		return ret - rhs;
	else if (negative && !rhs.negative)// -7 % 5 = 3 
		return rhs - ret;
	else if (negative && rhs.negative)// -7 % -5 = -2
		return -ret;
	else
		return ret;// 7 % 5 = 2

}
/********the number of bits in BigInteger*******************************/
size_t BigInteger::numBits() const {
	auto len = length();
	auto val = (*s)[len - 1];
	size_t num = 0;
	while (val)
	{
		val >>= 1;
		num += 1;
	}

	return (len - 1)*WIDTH + num;
}

size_t BigInteger::length()const { return s->size(); }
void  BigInteger::push_back(const eletype& x) { s->push_back(x); }
void BigInteger::resize(size_t size){ s->resize(size); }
void BigInteger::reserve(size_t size){ s->reserve(size); }

// get part of the BigInteger roughly, per 32 bits as one element
BigInteger BigInteger::sub_integer_rough(size_t start, unsigned int len)const{
	//note the start is logical concept, we need to recompute the position of start
	auto _length = length();
	assert(start>=0 && start < _length && _length != 0);

	BigInteger ret;

	size_t new_end = _length - start;
	size_t new_start;
	if (len == npos || len >= new_end)
		new_start = 0;
	else
		new_start = new_end - len;

	auto iter    = s->begin() + new_start;
	auto e_iter = s->begin() + new_end;
	ret.s = new datatype (iter, e_iter);
	return ret;
}


// get decimal representation by double dabble algorithm, for small BigInteger
BigInteger::operator std::string() const{
	return str();
}
// get decimal representation by double dabble algorithm, for small BigInteger
string BigInteger::str() const {
	return str_from_vec();
}

/**********convert string to BigInteger, str() is a history function which need to be modified*********************************/
string BigInteger::toString(Radix r) const{
	string str;
	auto &ls = *s;
	auto len = length();
	if (r == Bin)
	{ 
		auto last_len = numBits() % WIDTH;
		str += convert_to_binary(ls[len - 1], last_len);
		for (size_t i  = len -1 ;i--;) 
			str += convert_to_binary(ls[i], 32);
	}
	else if (r == Dec){
		if (length() <= 200)
			return this->str();
		return convert_to_dec(*this);
	}

	return str;
}


bool BigInteger:: operator < (const BigInteger& other)const{
	
	if (negative != other.negative)
		return negative < other.negative;

	if (length() != other.length()){
		return length() < other.length();
	}
	
	for (int i = length() - 1; i >= 0; i--){
		if (  (*s)[i] != other[i])
			return (*s)[i] < other[i];
	}
	return false;
}

bool BigInteger:: operator >= (const BigInteger& other)const{
	return !(*this < other);
}

bool BigInteger:: operator > (const BigInteger& other)const{
	return other < *this;
}

bool BigInteger::operator==(const BigInteger& other)const{
	return negative == other.negative && *s == *(other.s);
	}
bool BigInteger::operator !=(const BigInteger& other)const{
	return !(*this == other);
}

bool BigInteger::operator==(long num)const{
	return (length() == 1) && ((negative? -1: 1) * (*s)[0] == num);
}
bool BigInteger::operator!=(long num)const{
	return !(*this == num);
}

/*******swap function, no copy*************/
inline void swap(BigInteger& lhs, BigInteger& rhs)
{
	using std::swap;
	swap(lhs.s, rhs.s);
	swap(lhs.negative, rhs.negative);
}

ostream& operator << (ostream& out, const BigInteger& big)
{
	assert(big.length() > 0);
	out << big.str();

	return out;
}

//optimized binary to integer, not sage
#define unloop(u) case (u) : x += (1<<(u-1))*(*str++ - '0')
#define expand_four(u) unloop(u);unloop(u-1);unloop(u-2);unloop(u-3)
#define expand_16(u) expand_four(u);expand_four(u-4);expand_four(u - 8); expand_four(u - 12)

inline unsigned long fast_btoi(const char* str, unsigned int len)
{
	unsigned long x = 0;
	switch (len)
	{
		expand_16(32);
		expand_16(16);
	}
	return x;
}
#define hex_unloop(u)case (u) : x += (1<<(4*(u-1)))*\
	(*str - (( (*str) >= '0' && (*str) <= '9') ? '0' : ('A' - 10))); \
	str++
#define hex_expand_four(u) hex_unloop(u);hex_unloop(u-1);hex_unloop(u-2);hex_unloop(u-3)

/********optimized hex to integer ,not safe******************/
inline unsigned long fast_htoi(const char* str, unsigned int len)
{
	unsigned long x = 0;
	switch (len)
	{

		hex_expand_four(8);
		hex_expand_four(4);
	}
	return x;
}

//optimized atoi,not safe
inline unsigned long long fast_atoi(const char* str, unsigned int len)
{
		unsigned long long value_ = 0;
		switch (len) { // handle up to 10 digits, assume we're 32-bit
		case  10:   value_ += (*str++ - '0') * 1000000000;
		case  9:    value_ += (*str++ - '0') * 100000000;
		case  8:    value_ += (*str++ - '0') * 10000000;
		case  7:    value_ += (*str++ - '0') * 1000000;
		case  6:    value_ += (*str++ - '0') * 100000;
		case  5:    value_ += (*str++ - '0') * 10000;
		case  4:    value_ += (*str++ - '0') * 1000;
		case  3:    value_ += (*str++ - '0') * 100;
		case  2:    value_ += (*str++ - '0') * 10;
		case  1:    value_ += (*str++ - '0');
		}
		return value_;
}
BigInteger operator - (const BigInteger& big)
{
	BigInteger reverse = big;
	reverse.negative ^= 1;
	return reverse;
}

BigInteger::eletype& BigInteger::operator[](size_t index){
	return (*s)[index];
}
BigInteger::eletype BigInteger::operator[](size_t index)const{
	return (*s)[index];
}


/*construct data from binary representation*/
void BigInteger::parse_binary(const string& str,Radix radix)
{
	typedef unsigned long(*pfunType)(const char*, unsigned int);

	static const int steps[] = {32, 8}; 
	static const int perbits[] = { 1,4 };
	static const pfunType convert[] = {fast_btoi,fast_htoi};
	auto step = steps[radix], perbit = perbits[radix];
	auto fast_convert = convert[radix];
	
	s->clear();
	s->reserve(str.size() / (WIDTH/perbit));
	
	int len = step;
	int start = str.length();
	auto ch = str.c_str();

	while (start  > len)
	{
		start -= len;
		s->push_back(fast_convert(ch + start, len));//store 32 bits
	}
	s->push_back(fast_convert(ch, start));// stor remain bits
}

BigInteger BigInteger::truncate_mul(const BigInteger& rhs)
{

	auto llen = length();
	auto rlen = rhs.length();
	BigInteger ret;

	if (max(llen, rlen) < 10)
		return *this*rhs;

	int minbits = 1;
	int lnumbits = numBits(), rnumbits = rhs.numBits();
	size_t lalpha = max(minbits,lnumbits - 100);
	size_t ralpha = max(minbits, rnumbits - 100);

	ret = ((*this)>>lalpha)*(rhs>>ralpha);
	ret <<= lalpha+ralpha;
	return ret;
}

string  BigInteger::str_from_vec()  const
{
	int n = length();
	int nbits = 32 * n;         /* length of arr in bits */
	int nscratch = nbits / 3;   /* length of scratch in bytes */
	string scratch(1 + nscratch, 0);
	int i, j, k;
	int smin = 1;    /* speed optimization */

	for (i = n - 1; i >= 0; --i) {
		for (j = 31; j >= 0; --j) {
			int shifted_in = ((*s)[i] & (1 << j)) ? 1 : 0;//move

			for (k = smin; k >= 0; --k)
				scratch[k] += (scratch[k] >= 5) ? 3 : 0;//add 3 when larger than 5

			/**shift right one bit**/
			if (scratch[smin] >= 8)
				smin += 1;
			for (k = smin; k > 0; --k) {
				scratch[k] <<= 1;
				scratch[k] &= 0xF;//maintain 4 low digits
				scratch[k] |= (scratch[k - 1] >= 8);//consider the carry from low digits
			}

			/* 处理 */
			scratch[0] <<= 1;
			scratch[0] &= 0xF;
			scratch[0] |= shifted_in;
		}
	}

	/* Convert the scratch space from BCD digits to ASCII. */
	for (k = nscratch - 1; k >= 0; --k)
		scratch[k] += '0';
	for (k = nscratch - 1; k >= 0 && scratch[k] == '0'; --k){}
	k = nscratch - k;
	if (negative)
	{
		*(scratch.rbegin() + k - 1) = '-';
		k -= 1;
	}
	if (k == scratch.length())
		k -= 1;
	return string(scratch.rbegin() + k, scratch.rend());
}
		 
void BigInteger::	set_resversable(){
	negative = !negative;
}

/*************remove redundant zeros*************************/
void inline BigInteger::remove_redundant_zeros()
 {

	if ((*s)[length()-1] != 0)
		return;

	auto iter = find_if(s->rbegin(), s->rend(), [&](int x){return x != 0; });
	auto len = distance(iter, s->rend());

	if (len == 0)
	{
		s->erase(s->begin() + 1, s->end());	//corner case : contain at least one 0
		negative = false;// corner case: zero should be non-negative
	}
	else
		s->erase(s->begin() + len, s->end());
 }

/********convert long long integer to binary string************************************/
string convert_to_binary(unsigned long long val,size_t len)
{
	string ret(len,'0');
	int i =  len;
	while (i){
		ret[--i] = val & 1 ? '1':'0';
		val >>= 1;
	}
	return ret;
}
/******convert small decimal string to binary*********************************************************/
string convert_to_binary(string::const_iterator _s,size_t len)
{

	const static string str[] = { "0", "1" };
	string s(len, '0');
	string binary;
	binary.reserve(len * 3);
	int i = 0;
	for_each(_s, _s + len, [&](const char& c){s[i++] = c - '0'; });

	while (s != "0")//simulate divide by 2
	{
		int t = 0, old_t = 0;
		for (auto& ch : s)
		{
			t = ((old_t * 10 + ch) & 1);
			ch = (ch + old_t * 10) >> 1;
			old_t = t;
		}
		binary += str[t];
		if (s[0] == 0)
			s = remove_pre_zero(s);
	}
	reverse(binary.begin(),binary.end());
	return binary;	
}
/***************assistant function to convert huge BigInteger to decimal***************************/
string  convert_to_dec(const BigInteger& rhs)
{
	BigInteger ret(rhs);
	string str(rhs.numBits()/3,'\0');
	size_t len = min(unsigned(20000), 11 * ret.length());
	BigInteger divisor = BigInteger("1"+string(len,'0'));
	size_t top = str.size() - 1;
	while ( ret != 0 )
	{
		BigInteger remain,_ret;
		div(ret, divisor, _ret, remain,true);
		BigInteger q(remain);
		do
		{
			 auto _q = q / 1000000000;
			 string r = (q - _q * 1000000000).str();
			 q = _q;
			 for (auto iter = r.rbegin(); iter != r.rend(); ++iter)
				 str[top--] = *iter;
		} while (q != 0);
		ret = _ret;
	}
	str = remove_pre_zero(str);
	return str;
 }
 /*************get radix, negative of str**************************/
 inline void get_type(const string& str,Radix& radix,bool& negative,int& beg)
 {
	 beg = 0;
	 auto _len = str.length();

	 if (str[0] == '-'){
		 beg += 1;
		 negative = true;
	 }

	 if (_len >= 3 && str.substr(beg, 2) == "0b"){
		 radix = Bin;
		 beg += 2;
	 }
	 else if (_len >= 3 && str.substr(beg, 2) == "0x") {
		 radix = Hex;
		 beg += 2;
	 }
		 else 	 
			 radix = Dec;
	 }
 
/*************************************************/
 // find approriate repricoal to compute fixed division
 /*************************************************/
 BigInteger find_appropriate_rep(size_t& k,size_t len)
 {
	 for (auto& pair : BigInteger::PreRepricoal[len])
	 {
		 size_t bits = pair.first;
		 if (bits >= k)
		 {
			 k = bits;
			 return pair.second;
		 }
	 }
	 return 0;
 }
