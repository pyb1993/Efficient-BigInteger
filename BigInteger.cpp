#include "BigInteger.h"
#include "StringBigdata.h"
#include "sstream"
#include "cassert"
#include "algorithm"
#include "bitset"
#include "FFTMultiplier.h"
#include <utility>

using namespace std;

BigInteger::BigInteger()
{
	s = new datatype();
}

BigInteger::BigInteger(const BigInteger& rhs) :\
s(new datatype(*rhs.s)), negative(rhs.negative){   }

//接收FFT的结果
BigInteger::BigInteger(const vector<double>& Coef, size_t len,short other_width)//用来接收一个数组,专门处理FFT的结果
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

	if (length() == 0)//cornercase 
		s->push_back(0);

}

//移动构造函数
BigInteger::BigInteger(BigInteger&& rhs) NOEXCEPT  :
negative(rhs.negative), s(rhs.s)
{
	rhs.s = nullptr;//这样对其调用析构函数是安全的
}

#if 1
BigInteger::BigInteger( unsigned long long& num)
{
	s = new datatype();
	*this = num; 
}
#endif
BigInteger::BigInteger(const  long long& num)
{
	s = new datatype();
	*this = num;
}

BigInteger::BigInteger(const string& str)
{
	s = new datatype();
	*this = str;
}

BigInteger& BigInteger:: operator= (const BigInteger& rhs){
	if (this != &rhs){
		BigInteger  temp (rhs);
	    //swap 资源,自动释放局部变量temp
		swap(*this, temp);
	}
	return *this;
}


BigInteger& BigInteger:: operator= ( long long num){
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

	size_t bound = 5001;
	Radix radix;
	int beg;

	/*******For short decimal string OR Binary string, directly conversion*******/
	get_type(str, radix, negative, beg);
	if (radix != Dec || str.length() < bound){
		
		if (radix == Dec)
			parse_binary(convert_to_bin( beg ? str.substr(beg,-1):str), Bin);
		else
			parse_binary(beg ? str.substr(beg,-1) : str,radix);
		remove_redundant_zeros();
		return *this;
	}
	
	/***********Deal with the Decimal string***************************/
	static const long ten_power[10]= { 1, 10, 100, 1000, 10000, 100000, \
		1000000, 10000000, 100000000, 1000000000 };
	size_t i = 0;
	size_t _len = 5000;
	auto strlen = str.length();

	if (str[0] == '-')
		++i;

	BigInteger& result = *this;
	result = str.substr(i, _len);
	i += _len;

	while (i + _len < strlen)
	{
		result = result * BigInteger("1"+string(_len,'0')); //* _len
		//result = result *ten_power[_len]; //* 10
		result += BigInteger(str.substr(i, _len));
		i += _len;
	}

	//处理余下的计算
	if (i  < strlen)
	{
		auto remain = str.substr(i, -1);
		result = result * BigInteger("1" + string(remain.length(), '0')); //* 1000
		//result = result* ten_power[remain.length()];
		result += remain;
	}
		remove_redundant_zeros();
		return *this;
}

//移动赋值运算符
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


BigInteger BigInteger:: operator + (const long long& num) const
{
	assert(num >= 0);//暂时只处理正数
	
	BigInteger ret;
	unsigned long long  carry = num;//long long carry转换为无符号时有问题!
	unsigned long long  tmp = 0;
	auto len = length();
	auto &ls = *s;
	ret.resize(len);

	size_t i = 0;
	while (i < len)
	{
		tmp = carry + ls[i];
		ret[i++] =  tmp & Max;
		carry = tmp >> WIDTH;
	}
	if (carry)
		ret.push_back(carry);

	return ret;
}

BigInteger BigInteger:: operator + (const BigInteger& rhs) const{

	int llen = length(), rlen = rhs.length();

	if (llen < rlen)
		return rhs + *this;
	
	maxtype x = 0;
	maxtype carry = 0;
	BigInteger ret;
	
	ret.resize(llen);
	auto &ls = *s;
	auto &rs = *rhs.s;

	int i = 0;
	while (i < rlen){
		x = carry + ls[i] + rs[i];
		carry = x >> WIDTH;
		ret[i++] = x & Max;
	}

	while ( i < llen){
	
		x = carry + (ls)[i];
		carry = x >> WIDTH;
		ret[i++] = x & Max;
	}

	if (carry)
		ret.push_back(carry);
	
	if (ret[ret.length() - 1] == 0)
		ret.remove_redundant_zeros();
	return ret;
}

void BigInteger::operator += (const BigInteger& rhs)
{

	int llen = length(), rlen = rhs.length();
	auto maxlen = max(llen, rlen);
	auto minlen = min(llen, rlen);
	s->resize(maxlen);

	maxtype x = 0;
	maxtype carry = 0;

	auto &ls = *s;
	auto &rs = *rhs.s;

	int i = 0;
	while (i < minlen){
		x = carry + ls[i] + rs[i];
		carry = x >> WIDTH;
		ls[i++] = x & Max;
	}

	if (rlen < llen)
	{
		while (i < llen){
			x = carry + (ls)[i];
			carry = x >> WIDTH;
			ls[i++] = x & Max;
		}
	}
	else// rlen >= llen
	{
		while (i < rlen){
			x = carry + (rs)[i];
			carry = x >> WIDTH;
			ls[i++] = x & Max;
		}
	}

	if (carry)
		ls.push_back(carry);

	if (ls[length() - 1] == 0)
		remove_redundant_zeros();
}



 BigInteger subtract(const BigInteger& lhs, const BigInteger& rhs, bool neg)
{
	typedef BigInteger::maxtype maxtype;
	//assert(lhs >= 0 && rhs >= 0);//暂时只处理正数
	
	auto llen = lhs.length();
	auto rlen = rhs.length();

	if (lhs < rhs)
		return subtract(rhs, lhs, true);

	auto Max = BigInteger::Max;
	BigInteger ret;
	ret.resize(llen);

	long long subtracted;
	long long carry = 0;//因为元素是无符号,所以这里不能出现负号
	auto &ls = *lhs.s;
	auto &rs = *rhs.s;

	size_t i = 0;
	while(i < rlen){	
		subtracted = (long long)ls[i] - (long long)rs[i] - carry ;//无符号会导致溢出,必须转换64
		carry = subtracted < 0;
		ret[i++] = subtracted & Max;
	}

	while (i < llen){
		subtracted =  (long long)lhs[i] - carry;
		carry = subtracted < 0;
		ret[i++] = subtracted & Max;
	}
	if (neg)
		ret.set_resversable();

	return ret;
}

BigInteger BigInteger::operator - (const BigInteger& rhs) const
{
    //暂时只处理正数
	auto ret = subtract(*this, rhs);
	if (ret[ret.length() - 1] == 0)
		ret.remove_redundant_zeros();
	return ret;
}


BigInteger BigInteger:: operator *(const long & num) const{
	assert(num >= 0);

	auto len = length();
	BigInteger ret;
	ret.reserve(len);
	unsigned long long carry = 0;
	unsigned long long tmp = 0;
	for (size_t i = 0; i < len; ++i){
		tmp = carry + (unsigned long long)num * (*s)[i];//不会溢出
		carry = tmp >> WIDTH;
		ret.push_back(tmp & Max);
	}
	while (carry){
		ret.push_back(carry);
		carry >>= WIDTH;
	}
	return ret;
}

#if 1
BigInteger BigInteger::operator *(const BigInteger& rhs) const
{

	//乘法的KR算法
	auto llen = length(), rlen = rhs.length();
	assert(llen > 0 && rlen > 0);
	if (llen < 65536 && rlen < 65536)//avoid error from FFT
	{
		auto F = FFTMultiPlier(*this, rhs);
		return F.MulWithFFT();
	}

	//确保a大于b
	if (llen < rlen){
		return rhs* (*this);
      }
		
	if (llen== 1 && rlen == 1){
		auto x = maxtype((*s)[0]) * maxtype(rhs[0]);
		return x;
		}
	if (*this == BigInteger("0") || rhs == BigInteger("0"))
		return BigInteger("0");

		BigInteger ret;
		
		if (rlen * 2 > llen){
			size_t sub_len = llen >> 1;   //低位的长度,确保低位长度一致
			
			size_t mid1 = llen - sub_len, mid2 = rlen - sub_len;//注意这里的trick
			auto a1 =  sub_integer_rough(0, mid1),      a0 = sub_integer_rough(mid1);
			auto b1 = rhs.sub_integer_rough(0, mid2), b0 = rhs.sub_integer_rough(mid2);
			size_t weight =  sub_len * WIDTH;//权重 = 32*sublen

			auto z2 = a1 * b1;// z2 = a1*b1
			auto z0 = a0 * b0;// z0 = a0*b0

			auto middle_res = (a1+ a0)*(b1+b0);//(a1+a0)(b1+b0)

			auto z1 = middle_res - z2 - z0;

			z1	=  z1 << (weight);
			ret = z2 << (weight << 1);
		    ret	+= z1;
			ret  += z0;
		}
		else{//b的位数太小,不足以分割
			size_t mid = llen >> 1;
			auto a1 = sub_integer_rough(0, mid), a0 = sub_integer_rough(mid);
			size_t weight = a0.length() * WIDTH;//权重 = 32*低位长度

			auto z1 = a1 * rhs;
			auto z0 = a0 * rhs;
			ret =  z1 << weight;
			ret += z0;
		}
		ret.remove_redundant_zeros();
		return ret;
}
#endif

void BigInteger:: operator *= (const BigInteger& rhs){
	*this = *this * rhs;
}




BigInteger BigInteger:: operator <<(size_t shift_num) const
{

	if (*this == 0)
		return BigInteger(0);

	unsigned int shift_length= shift_num >>5 ;//整体移动的距离
	unsigned short shift_remain= shift_num - shift_length*WIDTH;//剩下的移动的距离
	BigInteger ret(*this);
	auto len = length();
	ret.resize(shift_length + len);
	if (shift_length)
	{
		for (int j = len - 1; j >= 0; --j)//整体移动
			ret[j + shift_length] = ret[j];

		for (size_t j = 0; j < shift_length; ++j)
			ret[j] = 0;
	}
	int shift_amount = shift_remain;
	unsigned long long carry = 0;

	auto i = shift_length;
	for (; i < ret.length(); ++i)
	{
		unsigned long  long val = ( unsigned long long(ret[i]) << shift_amount) + carry;//long long转换,为了多移位
		ret[i] =  val & Max;
		carry = val >> WIDTH;
   }
	while (carry > 0){
		ret.push_back(carry & Max);
		carry >>= WIDTH;
	}
	return ret;
}

BigInteger BigInteger:: operator >>(size_t shift_num) const{

	return BigInteger();
}


#if 0
BigInteger BigInteger::Inverse() const{

	auto len = length();
	if (len == 1) return BigInteger(1 / (*s)[0]);

	double InvBase = 1. / BASE;
	double x = (*s)[len - 1] + InvBase* (*s)[len - 2];
	x = BASE / x;
	
	auto tmp = BigInteger(*this)*BigInteger(x);


}
#endif
#if 0
BigInteger BigInteger::operator / (const BigInteger& rhs) const{

	
	B->Coef[1] = floor(x);
	B->Coef[0] = floor((x - B->Coef[1])*BASE);
	B->Size = 2;
}
#endif

#if 0//太慢了
BigInteger BigInteger::operator %(const BigInteger& rhs) const{

	//默认a,b不存在前置的0
	assert(*this > 0 && rhs > 0);// 均为正数
	auto dllen = digits_len(), drlen = rhs.digits_len();//位数
	if (*this < rhs)//a一定小于b
		return *this;
	
	bool k = 0;
	auto temp = sub_integer(0, drlen);
	BigInteger record[2] = { temp, temp };
	for (auto i = drlen; i < dllen + 1; ++i)
	{
		while (1)//不断的试探
		{
			record[k] = record[!k] - rhs;// tmp - b
			if (record[k].negative == true)
				break;
			k = !k;
		}
	if (i != dllen)
			record[!k].append_digits(this->sub_integer(i, 1));//这里的开销是O(n)....
	}
	return record[!k];
}
#endif

BigInteger BigInteger::operator %(const BigInteger& rhs) const{
	auto s1 = *this;
	auto s2 = rhs;
	auto ret = BigInteger(Mod(s1,s2));
	ret.remove_redundant_zeros();
	return ret;
}

size_t BigInteger::digits_len()const{
	auto len1 =  (length() - 1) * WIDTH;
	auto remain = to_string((*s)[length()-1]).length() % WIDTH;
	return len1 + remain;
}
size_t BigInteger::length()const { return s->size(); }
void  BigInteger::push_back(const eletype& x) { s->push_back(x); }
void BigInteger::resize(size_t size){ s->resize(size); }
void BigInteger::reserve(size_t size){ s->reserve(size); }

BigInteger BigInteger::sub_integer_rough(size_t start, unsigned int len)const{
	//注意这里的 start是逻辑上的,我们需要倒过来计算
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

#if 1
BigInteger BigInteger::sub_integer(size_t start, size_t len)const{
	//最快的实现方式
	auto str = this->str();
	return BigInteger(str.substr(start, len));
}
#else
BigInteger BigInteger::sub_integer(size_t start, size_t len)const{
	//最快的实现方式
	auto str = this->str();
	return BigInteger(str.substr(start, len));
}
#endif

#if 0
void BigInteger::append_digits(const string& _str){
	string str = *this;
	str += _str;
	this->parse_binary(str);//避免重复拷贝_str
}

	void BigInteger::append_zeros(size_t num){
	append_digits(string(num, '0'));
	}
#endif

//用来强制转换  关系是 << 用来检测 str()是不是
BigInteger::operator std::string() const{
	return str();
}

string BigInteger::str() const {
	return str_from_vec();
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

//经过优化的字符转换,不考虑安全性
#define unloop(u) case (u) : x += (1<<(u-1))*(*str++ - '0')
#define expand_four(u) unloop(u);unloop(u-1);unloop(u-2);unloop(u-3)
#define expand_16(u) expand_four(u);expand_four(u-4);expand_four(u - 8); expand_four(u - 12)

inline unsigned long fast_btoi(const char* str, unsigned char len)
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


inline unsigned long fast_htoi(const char* str, unsigned char len)
{
	unsigned long x = 0;
	switch (len)
	{

		hex_expand_four(8);
		hex_expand_four(4);
	}
	return x;
}



//经过优化的atoi,不考虑安全性
inline unsigned long long fast_atoi(const char* str, int len)
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
BigInteger operator - (const BigInteger& big){
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
	typedef unsigned long(*pfunType)(const char*, unsigned char);

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

	while (start  > len){
		start -= len;
		s->push_back(fast_convert(ch + start, len));//一次存32bit
	}
	s->push_back(fast_convert(ch, start));
	

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
			int shifted_in = ((*s)[i] & (1 << j)) ? 1 : 0;//移动的距离

			for (k = smin; k >= 0; --k)
				scratch[k] += (scratch[k] >= 5) ? 3 : 0;//超过5就要加3

			/*原始的binary向右移动 1 bit*/
			if (scratch[smin] >= 8)
				smin += 1;
			for (k = smin; k > 0; --k) {
				scratch[k] <<= 1;
				scratch[k] &= 0xF;//保留低四位
				scratch[k] |= (scratch[k - 1] >= 8);//考虑低位的进位
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
	
	  


  unsigned short BigInteger::highest_bitsnum() const{
	  auto x = (*s)[length()];
	  unsigned short numbits = 0;
	  while (x){
		  x >>= 1;
		  numbits += 1;
	  }
	  return numbits;
  }

void BigInteger::	set_resversable(){
	negative = !negative;
}

void inline BigInteger::remove_redundant_zeros()//去掉多余的0
 {
	//这里有个细节;至少要保留一个0
	auto iter = find_if(s->rbegin(), s->rend(), [&](int x){return x != 0; });
	auto len = distance(iter, s->rend());
	if (len == 0)
		s->erase(s->begin() + 1, s->end());
	else
		s->erase(s->begin() + len, s->end());
 }

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


inline int get_min_binary_power(int x){
	assert(x >= 1);
	int i = 1;
	do
	{
	  i *= 2;
	} while (i <= x);

	return i >> 1;
}
 inline void bin_shift_multiply(BigInteger& result, const string& binary)
{
	 auto blen = binary.length();
	 auto base = result;
	 result = BigInteger(0);
	 for (size_t k = 0; k < blen; ++k){
		 if (binary[k] == '1')
			 result = result + (base << (blen - k - 1));
	 }
}

 /*************get radix, is negative of str**************************/
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
 
