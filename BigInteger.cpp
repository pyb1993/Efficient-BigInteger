#include "BigInteger.h"
#include "StringBigdata.h"
#include "sstream"
#include "cassert"
using namespace std;

BigInteger::BigInteger(){
	s = new datatype();
}

BigInteger::BigInteger(const BigInteger& rhs) : s(new datatype(*rhs.s)), negative(rhs.negative)
{}

BigInteger::BigInteger(const long long& num){ 
	s = new datatype();
	*this = num; 
}
BigInteger::BigInteger(const string& str){
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
	return ret;
}
BigInteger  subtract(const BigInteger& lhs, const BigInteger& rhs, bool neg)
{
	assert(lhs >= 0 && rhs >= 0);//暂时只处理正数
	//首先消除前置的0;
	auto BASE = BigInteger::BASE;
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
	return subtract(*this, rhs);
}


#if 0
BigInteger BigInteger::operator *(const BigInteger& rhs) const{

	//乘法的分治算法
	auto llen = length(), rlen = rhs.length();
	assert(llen > 0 && rlen > 0);

	//确保a大于b
	if (llen < rlen){
		return rhs* (*this);
      }
		if (llen== 1 && rlen == 1){
			return long long( s[0] * rhs[0]);
		}

		BigInteger s;
		if (llen * 2 > rlen){
			size_t weight = llen >> 1;
			size_t mid1 = llen - weight, mid2 = rlen - weight;//注意这里的trick

			auto a1 =  sub_integer_rough(0, mid1),      a0 = sub_integer_rough(mid1);
			auto b1 = rhs.sub_integer_rough(0, mid2), b0 = rhs.sub_integer_rough(mid2);

			auto z2 = a1 * b1;// z2 = a1*b1
			auto z0 = a0 * b0;// z0 = a0*b0

			auto middle_res = (a1+ a0)*(b1+b0);//(a1+a0)(b1+b0)


			//auto z1 = Subtract(Subtract(middle_res, z2), z0);//z1 = middle_res - z2 - z0，

			//assert(z1[0] != '-' || !(cout << " (" <<a1+"+"+a0 << " ) * ("<<b1+"+"+b0 <<") = "<<middle_res<<endl<<\
						//middle_res<<" - "<< z2+"- "+z0 + " :(" +a0+"* "+b0+ ")"<<" = "<< z1<<endl));

			s = z2 + string(weight << 1, '0');// z2 * 10^(2*mid)
			s = add(s, (z1 + string(weight, '0')));// z1  * 10^mid
			s = add(s, z0);//z0
		}
		else{//b的位数太小,不足以计算
			size_t mid = a_len >> 1;
			size_t weight = a_len - mid;
			string a1 = a.substr(0, mid), a0 = a.substr(mid);
			auto z1 = Effective_multiply(a1, b);
			auto z0 = Effective_multiply(a0, b);
			s = z1 + string(weight, '0');
			s = add(s, z0);
		}
		return s;
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
	return BigInteger(Mod(s1,s2));
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

BigInteger BigInteger::sub_integer_rough(size_t start, size_t len)const{
	//注意这里的 start是逻辑上的,我们需要倒过来计算
	int new_start = length() - start - len;
	BigInteger ret;
	if (len+start > length() ){
		len = length();
	}
	auto iter    = s->begin() + new_start;
	auto e_iter = iter + len;
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