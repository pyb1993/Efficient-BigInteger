
#include "BigInteger.h"
#include "StringBigdata.h"
#include "sstream"
#include "cassert"
#include "algorithm"
#include "bitset"
#include "FFTMultiplier.h"
#include <utility>

using namespace std;

map<size_t, map<size_t, BigInteger>> BigInteger::PreRepricoal = map<size_t, map<size_t, BigInteger>>();

BigInteger::BigInteger()
{
	s = new datatype();
}

BigInteger::BigInteger(const BigInteger& rhs) :\
s(new datatype(*rhs.s)), negative(rhs.negative){   }

//convert data fom FFT coeficients to BigInteger
BigInteger::BigInteger(const vector<double>& Coef, size_t len, unsigned int other_width)
:s(new datatype)
{
	s->reserve(len >> 2);
	maxtype carry = 0;
	size_t i = 0;
	for (; i + 3 < len; i += 4){
		maxtype x = carry + (maxtype)Coef[i];// 1
		x += ((maxtype(Coef[i + 1])) << (other_width));// 8
		x += ((maxtype(Coef[i + 2])) << (2 * other_width));//16
		x += ((maxtype(Coef[i + 3])) << (3 * other_width));//24

		if (x > BASE)
		{
			carry = x >> WIDTH;
			x &= Max;
		}
		s->push_back(x & Max);
	}
	for (int k = 0; i < len; ++i, ++k)
		carry += (maxtype)Coef[i] << (other_width*k);

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

BigInteger::BigInteger(unsigned long long& num)
{
	s = new datatype();
	*this = num;
}

BigInteger::BigInteger(const  long long& num)
{
	s = new datatype();
	*this = num;
}

BigInteger::BigInteger(const string& str, size_t sublen) :
substr_len(sublen)
{
	s = new datatype();
	*this = str;
}

BigInteger& BigInteger::operator= (const BigInteger& rhs){
	if (this != &rhs){
		BigInteger  temp(rhs);
		//swap 资源,自动释放局部变量temp
		swap(*this, temp);
	}
	return *this;
}


BigInteger& BigInteger::operator= (long long num){
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
BigInteger& BigInteger:: operator= (unsigned long long& _num){
	auto num = _num;
	do{
		s->push_back(num & Max);
		num >>= WIDTH;
	} while (num > 0);
	
	return *this;
}



BigInteger fromDec(const string& str,size_t sublen = 16384){
	
	BigInteger result;
	size_t str_len = str.length();
	static map<size_t, BigInteger> TenPower = map<size_t, BigInteger>();

	if (str.length() <= 2048)// directly parse for small decimal string
		return fromBinary(convert_to_binary(str), Bin);
	
	else
	{
	/***********Deal with the Long Decimal string***************************/
		size_t i = 0;
		if (str[0] == '-') ++i;// skip sign

		result = fromDec(str.substr(i, sublen), sublen >> 1);
		i += sublen;

	// each iteration, we convert sublen digits to BigInteger by ret = ret*(10.....0) + remain
	// we store the BigInteger 10000...0 when we compute it for reuse
		while (i + sublen < str_len)
		{
			if (TenPower.find(sublen) == TenPower.end()){
				TenPower[sublen] = fromDec("1" + string(sublen, '0'), sublen >> 1);//recursively convert,with half sublen
			}

			result *= TenPower[sublen];
			result += fromDec(str.substr(i, sublen), sublen >> 1);
			i += sublen;
		}

		if (i  < str_len)
		{
			auto remain = str.substr(i, -1);
			result *= fromDec("1" + string(remain.length(), '0'), sublen >> 1); //* 1000
			result += remain;
		}
		result.remove_redundant_zeros();
		return result;
	
	}
}

BigInteger& BigInteger:: operator= (const string& str){

	/****Deal with the Binary string :bin,hex******/
	int str_len = str.length();
	assert(str_len != 0);

	size_t bound = 1025;
	Radix radix;
	int start;

	/*******For short decimal string OR Binary string, directly conversion*******/
	get_type(str, radix, negative, start);

	if (radix == Dec){
		auto ng = negative;
		*this = fromDec(str.substr(start, str.length() - start));
		negative = ng;//recover the negative
	}

	else {
		auto ng = negative;
		*this = fromBinary(start ? str.substr(start, -1) : str, radix);
		negative = ng;//recover the negative
	}

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
			add(lhs, -rhs, ret);// + - - = + + +
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

/********help function for operator -*******************************/
void subtract(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret, bool neg)
{
	typedef BigInteger::maxtype maxtype;
	auto base = BigInteger::BASE;

	auto llen = lhs.length();
	auto rlen = rhs.length();

	if (lhs < rhs)
	{
		ret = rhs;
		subtract(rhs, lhs, ret);
		ret.set_resversable();
		return;
	}
		
		auto& ls = *lhs.s;
		auto& rs = *rhs.s;

		int i = 0;
		unsigned long borrow = 0;
		//use a trick: the overflow will automatically do work of borrow
		for (; i < rlen; i++)
		{
			ret[i] = ls[i] - rs[i] - borrow;// imporvement to long long
			borrow = (ls[i] < rs[i]);
		}
		if (borrow)
		{
			for (; ls[i] == 0; i++) --ret[i];
			--ret[i];
		}
		ret.remove_redundant_zeros();
	}


BigInteger BigInteger::operator - (const BigInteger& rhs) const
{
	BigInteger ret(*this);
	_subtract(*this, rhs, ret);
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
	auto Max = BigInteger::Max;
	typedef BigInteger::eletype eletype;

	auto len = lhs.length();
	auto &ls = *(lhs.s);

	ret.resize(len);
	unsigned long long carry = 0;
	unsigned long long tmp = 0;
	size_t i = 0;
	for (auto& val : ls){
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
		ret = BigInteger("0");

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

BigInteger BigInteger:: operator <<(size_t shift_num) const
{
	BigInteger ret(*this);
	ret <<= shift_num;
	return ret;
}

void BigInteger:: operator <<= (size_t shift_num)
{
	auto& lhs = *s;
	size_t shift_length = shift_num >> 5;//move long distance
	size_t shift_remain = shift_num & (WIDTH-1);//move remain distance
	size_t inv_shift = WIDTH - shift_remain;
	auto len = length();
	resize(shift_length + len);
	
	for (size_t j = len; shift_length && j--;)
		lhs[j + shift_length] = lhs[j];
	for (size_t j = 0; j < shift_length; ++j)
		lhs[j] = 0;

	len = length();
	size_t carry = 0,last = 0;

	for (size_t i = 0; i < len; ++i)
	{
		last = (lhs[i] >> inv_shift);
		lhs[i] = (lhs[i] << shift_remain) | carry;//long long转换
		carry = last;
	}

	if (carry > 0)
		lhs.push_back(carry);

}


BigInteger BigInteger:: operator >>(size_t shift_num) const
{
	BigInteger ret(*this);
	ret >>= shift_num;
	return ret;
}

void BigInteger:: operator >>=(size_t shift_num)
{

	auto &lhs = *s;
	size_t shift_length = shift_num >> 5;
	size_t shift_remain = shift_num & (WIDTH-1);
	size_t invshift = WIDTH - shift_remain;
	size_t mask = (1 << shift_remain) - 1;
	auto len = length();


	for (size_t j = shift_length; shift_length && j < len; ++j)//整体移动
		lhs[j - shift_length] = lhs[j];

	for (size_t j = len - shift_length ; j < len; ++j)
		lhs[j] = 0;
	len -= shift_length;


	for (int i = 0; i + 1 < len; ++i){
		lhs[i] >>= shift_remain;
		lhs[i] |= (lhs[i+1] & mask) << invshift;
	}

	lhs[len - 1] >>= shift_remain;
	remove_redundant_zeros();
}


BigInteger BigInteger::operator / (const BigInteger& rhs) const
{
	BigInteger q, ret;
	knuthDiv(*this, rhs, q, ret);
	q.negative = (negative != rhs.negative);
	return q;
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
	knuthDiv(*this, rhs, q, ret);

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

// get part of the BigInteger, per 32 bits as one element
BigInteger BigInteger::sub_integer_rough(size_t start, unsigned int len)const{
	//note the start is logical concept, we need to recompute the position of start
	auto _length = length();
	assert(start >= 0 && start < _length && _length != 0);

	BigInteger ret;

	size_t new_end = _length - start;
	size_t new_start;
	if (len == npos || len >= new_end)
		new_start = 0;
	else
		new_start = new_end - len;

	auto iter = s->begin() + new_start;
	auto e_iter = s->begin() + new_end;
	ret.s = new datatype(iter, e_iter);
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
string BigInteger::toString(Radix r) const
{
	string str;
	auto &ls = *s;
	auto len = length();
	if (r == Bin)
	{
		auto last_len = numBits() % WIDTH;
		str += convert_to_binary(ls[len - 1], last_len);
		for (size_t i = len - 1; i--;)
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
		if ((*s)[i] != other[i])
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
	return (length() == 1) && ((negative ? -1 : 1) * (*s)[0] == num);
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
	(*str - (((*str) >= '0' && (*str) <= '9') ? '0' : ('A' - 10))); \
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
 BigInteger fromBinary(const string& str, Radix radix)
{
	typedef unsigned long(*pfunType)(const char*, unsigned int);
	BigInteger ret;
	auto s = ret.s;

	static const int steps[] = { 32, 8 };
	static const int perbits[] = { 1, 4 };
	static const pfunType convert[] = { fast_btoi, fast_htoi };
	auto step = steps[radix], perbit = perbits[radix];
	auto fast_convert = convert[radix];

	s->clear();
	s->reserve(str.size() / (BigInteger::WIDTH / perbit));

	int len = step;
	int start = str.length();
	auto ch = str.c_str();

	while (start  > len)
	{
		start -= len;
		s->push_back(fast_convert(ch + start, len));//store 32 bits
	}
	s->push_back(fast_convert(ch, start));// stor remain bits
	ret.remove_redundant_zeros();
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

void BigInteger::set_resversable(){
	negative = !negative;
}

/*************remove redundant zeros*************************/
void inline BigInteger::remove_redundant_zeros()
{

	if ((*s)[length() - 1] != 0)
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
string convert_to_binary(unsigned long long val, size_t len)
{
	string ret(len, '0');
	int i = len;
	while (i){
		ret[--i] = val & 1 ? '1' : '0';
		val >>= 1;
	}
	return ret;
}

/******convert small decimal string to binary*********************************************************/
string convert_to_binary(string s)
{
	const static string str[] = { "0", "1" };
	string binary;
	binary.reserve(s.length() * 3);
	string tmps;
	for_each(s.begin(),s.end(), [&](const char& c){tmps.push_back( c - '0'); });

	while (tmps != "0")//simulate divide by 2
	{
		int t = 0, old_t = 0;
		for (auto& ch : tmps)
		{
			t = ((old_t * 10 + ch) & 1);
			ch = (ch + old_t * 10) >> 1;
			old_t = t;
		}
		binary += str[t];
		if (tmps[0] == 0)
			tmps = remove_pre_zero(tmps);
	}

	reverse(binary.begin(), binary.end());
	return binary;
}

/***************assistant function to convert huge BigInteger to decimal***************************/
string  convert_to_dec(const BigInteger& rhs)
{
	BigInteger ret(rhs);
	string str(rhs.numBits() / 3, '\0');
	size_t len = min(unsigned(20000), 11 * ret.length());
	BigInteger divisor = BigInteger("1" + string(len, '0'));
	size_t top = str.size() - 1;
	while (ret != 0)
	{
		BigInteger remain, _ret;
		knuthDiv(ret, divisor, _ret, remain);
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
inline void get_type(const string& str, Radix& radix, bool& negative, int& beg)
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


void knuthDiv(const BigInteger& dividend, const BigInteger& divisor, BigInteger& quotient, BigInteger& remainder)
{
	using maxtype = BigInteger::maxtype;
	auto base = BigInteger::BASE;
	auto width = BigInteger::WIDTH;
	auto Max = BigInteger::Max;
	if (divisor == 0){
		throw exception("dividened by zero");
	}
	if (dividend < divisor){
		quotient = 0;
		remainder = dividend;
		return;
	}
	else if (divisor.length() == 1){
		quotient = dividend / divisor[0];
		remainder = dividend - quotient*divisor[0];
		return;
	}

	// Hacker's Delight Implementation of knuth algorithm D
	//normalize to make the most significant element of divisor greater than base/2
	auto last = divisor.s->back();
	int cnt = 1;
	while (last = last >> 1) { ++cnt; }
	cnt = width - cnt;

	auto divd_norm = dividend << cnt;
	auto divs_norm = divisor << cnt;

	int m = divd_norm.length();
	int n = divs_norm.length();
	divd_norm.push_back(0);
	maxtype q, r;

	quotient.resize(m - n + 1);
	for (int j = m - n; j >= 0; --j)// Main Loop
	{
		// Compute estimate qhat of q[j].
		q = (divd_norm[j + n] * base + divd_norm[j + n - 1]) / divs_norm[n - 1];
		r = (divd_norm[j + n] * base + divd_norm[j + n - 1]) - q*divs_norm[n - 1];
		while (q >= base || q*divs_norm[n - 2] > base*r + divd_norm[j + n - 2])
		{
			q -= 1;
			r += divs_norm[n - 1];
			if (r >= base) break;
		}

		// Multiply and subtract.
		maxtype p = 0;
		long long t = 0, k = 0;

		for (int i = 0; i < n; i++) {
			p = q * divs_norm[i];
			t = divd_norm[i + j] - k - (p & Max);
			divd_norm[i + j] = t;
			k = (p >> width) - (t >> width);
		}

		t = divd_norm[j + n] - k;
		divd_norm[j + n] = t;
		// If we subtracted too much, add back.
		if (t < 0) {
			q -= 1;
			k = 0;
			for (int i = 0; i < n; i++) {
				t = divd_norm[i + j] + divs_norm[i] + k;
				divd_norm[i + j] = t;
				k = t >> width;
			}
			divd_norm[j + n] += k;
		}
		quotient[j] = q; // Store quotient digit.
	} // End j.

	// unnormalize the remainder and pass it back

	remainder.resize(n);
	for (int i = 0; i < n; i++)
		remainder[i] = (divd_norm[i] >> cnt) | (divd_norm[i + 1] << (width - cnt));

	quotient.remove_redundant_zeros();
	remainder.remove_redundant_zeros();
}
