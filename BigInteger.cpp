#include "BigInteger.h"
#include "StringBigdata.h"
#include "sstream"
#include "cassert"
#include "algorithm"
#include <utility>

using namespace std;

BigInteger::BigInteger(const int Base, const int Width):\
BASE(Base), WIDTH(Width)
{
	s = new datatype();
}

BigInteger::BigInteger(const BigInteger& rhs, const int Base, const int Width) :\
s(new datatype(*rhs.s)), negative(rhs.negative), BASE(Base),WIDTH(Width){   }

//接收FFT的结果
BigInteger::BigInteger(const vector<double>& Coef, size_t len, long long other_base,short other_width)//用来接收一个数组,专门处理FFT的结果
:BASE(DefaultBase), WIDTH(DefaultWIDTH), s( new datatype)
{
	s->reserve(len >> 1);
	long long carry = 0;
	int i = 0;
	for (; i+1 < len; i += 2){
		long long x = carry + Coef[i] + Coef[i+1] * other_base;
		if (x > BASE)
		{
			carry = x / BASE;
			x %= BASE;
		}
		s->push_back(x%BASE);
	}
	if (i != len) carry += Coef[len - 1];
	while (carry) {
		s->push_back(carry%BASE);
		carry /= BASE;
	}

}

//移动构造函数
BigInteger::BigInteger(BigInteger&& rhs, const int Base, const int Width) NOEXCEPT  :
negative(rhs.negative), s(rhs.s), BASE(Base), WIDTH(Width)
{
	rhs.s = nullptr;//这样对其调用析构函数是安全的
}

BigInteger::BigInteger(const long long& num, const int Base, const int Width):\
BASE(Base), WIDTH(Width)
{
	s = new datatype();
	*this = num; 
}
BigInteger::BigInteger(const string& str, const int Base, const int Width) :\
BASE(Base), WIDTH(Width)
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
	do{
		s->push_back(num%BASE);
		num /= BASE;
	} while (num > 0);
	return *this;
}

BigInteger& BigInteger:: operator= (const string& str){
	    int str_len = str.length();
	    assert(str_len != 0);
		vec_from_str(str);
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

BigInteger BigInteger:: operator + (const BigInteger& rhs) const{
	BigInteger ret;
	eletype x = 0;
	char count = 0;
	int llen = length(), rlen = rhs.length();
	ret.reserve(max(llen,rlen));

	for (int i = 0; ; ++i){
	
		if (count == 0 && i >= llen && i >= rlen)
			break;
		x = count + (i < llen ? (*s)[i] : 0) + ( i < rlen ? rhs[i] : 0);
		count = x / BASE;
		ret.push_back(x%BASE);
	}
	ret.remove_redundant_zeros();
	return ret;
}
BigInteger  subtract(const BigInteger& lhs, const BigInteger& rhs, bool neg)
{
	if (!(lhs >= 0 && rhs >= 0))
	{
		cout << 1;
	}
		assert(lhs >= 0 && rhs >= 0);//暂时只处理正数
		assert(lhs.BASE == rhs.BASE);
		
	auto BASE = rhs.BASE;
	BigInteger ret;
	auto llen = lhs.length();
	auto rlen = rhs.length();

	if (lhs < rhs){
		return subtract(rhs, lhs, true);
	}

	ret.reserve(llen);
	int carry = 0;
	for (size_t i = 0; i < llen; ++i){
		int subtracted = carry + lhs[i];
		subtracted -= i < rlen ? rhs[i] : 0;
		carry = 0;
		if (subtracted < 0){
			subtracted += BASE;
			carry = -1;
		}
		ret.push_back(subtracted);
		assert(subtracted  >= 0);
	}
	if (neg)
		ret.set_resversable();
	return ret;
}
BigInteger BigInteger::operator - (const BigInteger& rhs) const
{
    //暂时只处理正数
	auto ret = subtract(*this, rhs);
	ret.remove_redundant_zeros();
	return ret;
}

#if 1
BigInteger BigInteger::operator *(const BigInteger& rhs) const
{

	//乘法的KR算法
	auto llen = length(), rlen = rhs.length();
	assert(llen > 0 && rlen > 0);

	//确保a大于b
	if (llen < rlen){
		return rhs* (*this);
      }
		
	if (llen== 1 && rlen == 1){
		return  long long((*s)[0]) * long long(rhs[0]);
		}
	if (*this == BigInteger("0") || rhs == BigInteger("0"))
		return BigInteger("0");

		BigInteger ret;
		
		if (rlen * 2 > llen){
			size_t sub_len = llen >> 1;   //低位的长度,确保低位长度一致
			
			size_t mid1 = llen - sub_len, mid2 = rlen - sub_len;//注意这里的trick
			auto a1 =  sub_integer_rough(0, mid1),      a0 = sub_integer_rough(mid1);
			auto b1 = rhs.sub_integer_rough(0, mid2), b0 = rhs.sub_integer_rough(mid2);
			size_t weight =  sub_len * WIDTH;//权重 = 8*sublen

			auto z2 = a1 * b1;// z2 = a1*b1
			auto z0 = a0 * b0;// z0 = a0*b0

			auto middle_res = (a1+ a0)*(b1+b0);//(a1+a0)(b1+b0)

			auto z1 = middle_res - z2 - z0;

			z2.append_zeros(weight<<1);//z2 * 10^(2*mid)
			z1.append_zeros(weight);
			ret = z2 + z1 + z0;
		}
		else{//b的位数太小,不足以分割
			size_t mid = llen >> 1;
			auto a1 = sub_integer_rough(0, mid), a0 = sub_integer_rough(mid);
			size_t weight = a0.length() * WIDTH;//权重 = 8*低位长度

			auto z1 = a1 * rhs;
			auto z0 = a0 * rhs;
			z1.append_zeros(weight);
			ret = z1 + z0;
		}
		ret.remove_redundant_zeros();
		return ret;
}
#endif

BigInteger BigInteger:: operator <<(size_t shift_num) const{

	auto len = length();
	int shift_amount = 36;//trick 一次最多移动36位  2^36 < (2^64-1)/(1+BASE)
	BigInteger ret(*this);
	ret.reserve(len + 10);

	for (int count = shift_num; count > 0;	count -= shift_amount){

		shift_amount = min(count, shift_amount);
		unsigned long long carry = 0;
		int i ;
		for ( i = 0; i < ret.length(); ++i)
		{
			unsigned long  long val = ( unsigned long long(ret[i]) << shift_amount) + carry;//long long转换,为了多移位
			ret[i] =  val%BASE;
			carry = val / BASE;
		}
		while (carry > 0){
			ret.push_back(carry % BASE);
			carry /= BASE;
		}
	}
	return ret;
}



BigInteger BigInteger:: operator >>(size_t shift_num) const{

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
	assert(start>=0 && start< _length && _length != 0);

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

void BigInteger::append_digits(const string& _str){
 	//assert(length() != 0);
	string str = *this;
	str += _str;
	this->vec_from_str(str);//避免重复拷贝_str
}

	void BigInteger::append_zeros(size_t num){
	append_digits(string(num, '0'));
	}

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

inline void swap(BigInteger& lhs, BigInteger& rhs){
	using std::swap;
	swap(lhs.s, rhs.s);
	swap(lhs.negative, rhs.negative);
}

ostream& operator << (ostream& out, const BigInteger& big){
	assert(big.length() > 0);
	out << big.str();

	return out;
}

//经过优化的atoi,比原来的快很多
inline int fast_atoi(const char* str, int len)
{
		int value_ = 0;
		switch (len) { // handle up to 10 digits, assume we're 32-bit
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


void BigInteger::vec_from_str(const string& str){

	s->clear();
	s->reserve(str.size() / WIDTH);
	int len = WIDTH;
	int start = str.length();
	auto ch = str.c_str();
	if (str[0] == '-'){
		negative = true;
		++ch;
		start -= 1;
	}
	while (start  > len){
		start -= len;
		s->push_back(fast_atoi(ch+start, len));
	}
	s->push_back(fast_atoi(ch, start));
}

  string  BigInteger::str_from_vec()  const
{
	string str;
	str.reserve(length()*WIDTH);

	if (negative)
		str += "-";

	 str += to_string(s->back());//最高位不需要补充0

	for (int i = length() - 2; i >= 0; --i){
		auto tmp = to_string((*s)[i]);
		if (tmp.length() < WIDTH)
			str += string(WIDTH - tmp.length(), '0');
		str += tmp;
	}
	return str;
}

void BigInteger::	set_resversable(){
	negative = !negative;
}

void inline BigInteger::remove_redundant_zeros()//去掉前置的0
 {

	//这里有多个细节;至少要保留一个0
	if (length() == 1) return;//不需要处理

	int non_zero = 0;
	int len = length();

	auto iter = s->rbegin();
	while (non_zero < len && *iter++ == 0){
		non_zero += 1;
	}

	if (non_zero)
	{
		if (non_zero == len )
			non_zero -= 1;
		auto beg = s->begin() + len - non_zero ;
		s->erase(beg,s->end());
	}
 }



