#ifndef BIGINTEGERPYB
#define BIGINTEGERPYB

#include "iostream"
#include "string"
#include "vector"
#include "algorithm"
#include "memory"
#include "map"
#include "climits"

#ifndef _MSC_VER
#define NOEXCEPT noexcept
#else
#define NOEXCEPT
#endif

using namespace std;
enum Radix { Bin, Hex, Dec };
class BigInteger{
public:
	using eletype = unsigned long;
	using maxtype = unsigned long long;
	using datatype = vector<eletype>;
	static const eletype Max = ULONG_MAX;//2^32 - 1
	static const unsigned long long BASE = (unsigned long long)ULONG_MAX + 1;//2^32
	static const int WIDTH = 32;
	static const int x;
	static map<size_t, BigInteger> TenPower;// = { convert_to_binary(string("111").begin(), 100000) };
	static map<size_t, map<size_t, BigInteger>> PreRepricoal;


	friend class FFTMultiPlier;
	friend inline void swap(BigInteger& lhs, BigInteger& rhs);
	friend ostream& operator << (ostream& out, const BigInteger& big);
	friend BigInteger operator -(const BigInteger& big);

	friend void subtract(const BigInteger& lhs, const BigInteger& r, BigInteger& ret, bool neg);
	friend void mul(const BigInteger&, const BigInteger&, BigInteger&);
	friend void mul(const BigInteger&, long, BigInteger&);
	friend void add(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret);
	friend void add(const BigInteger& lhs, long, BigInteger& ret);
	friend void div(const BigInteger& lhs, const BigInteger& rhs, BigInteger& q, BigInteger& ret, bool byconst);
	friend void knuthDiv(const BigInteger& dividend, const BigInteger& divisor, BigInteger& quotient, BigInteger& remainder);
	friend string convert_to_dec(const BigInteger& rhs);
	friend BigInteger find_appropriate_rep(size_t& k, size_t len);
	friend void _add(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret);
	friend void _add(const BigInteger& lhs, long, BigInteger& ret);
	friend void _subtract(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret);
	friend BigInteger fromDec(const string&, size_t sublen);
	friend BigInteger fromBinary(const string& str, Radix radix);


	BigInteger();
	BigInteger(const BigInteger&);
	BigInteger(unsigned long long& num);//可能溢出,必须用unsigned long long
	BigInteger(const long long& num);

	BigInteger(const string&, size_t sublen = 20000);
	BigInteger(const vector<double>&, size_t len, unsigned int other_width = 8);//用来接收一个数组,专门处理FFT的结果
	BigInteger(BigInteger&&) NOEXCEPT;
	~BigInteger();


	BigInteger& operator= (const BigInteger& rhs);
	BigInteger& operator= (long long num);
	BigInteger& operator= (unsigned long long& num);

	BigInteger& operator= (const string& str);
	BigInteger& operator= (BigInteger&&)NOEXCEPT;
	BigInteger operator + (const BigInteger&) const;
	BigInteger operator + (const long long &) const;

	BigInteger operator - (const BigInteger&) const;
	BigInteger operator * (const BigInteger&) const;
	BigInteger operator * (const long &) const;

	BigInteger operator / (const BigInteger&) const;
	BigInteger operator / (const long) const;
	BigInteger operator % (const BigInteger&) const;
	BigInteger operator % (const long) const;
	BigInteger operator <<(size_t num) const;
	BigInteger operator >>(size_t num) const;

	void operator *=  (const BigInteger&);//此处是否应该返回
	void operator *=  (long);//此处是否应该返回
	void operator +=  (const BigInteger&);
	void operator +=  (long);
	void operator <<= (size_t num);
	void operator >>= (size_t num);
	void operator -=  (const BigInteger&);

	bool operator <  (const BigInteger& other)const;
	bool operator >  (const BigInteger& other)const;
	bool operator >= (const BigInteger& other)const;
	bool operator == (const BigInteger&) const;
	bool operator != (const BigInteger&) const;
	bool operator != (long) const;
	bool operator == (long) const;

	eletype& operator[] (size_t index);
	eletype operator [] (size_t index) const;
	size_t length() const;
	size_t numBits() const;
	void resize(size_t size);
	void reserve(size_t size);
	BigInteger pow(const long long&) const;
	void push_back(const eletype& x);
	string str() const;
	string toString(Radix r) const;

	BigInteger sub_integer_rough(size_t start, unsigned int len = npos) const;
	BigInteger Inverse() const;
	void set_resversable();
	operator string () const;

private:
	string str_from_vec() const;
	
	
	inline void remove_redundant_zeros();
	datatype* s = nullptr;
	bool negative = false;
	static const size_t npos = -1;
	size_t substr_len = 20000;

	


};

inline void swap(BigInteger& lhs, BigInteger& rhs);
ostream& operator << (ostream& out, const BigInteger& big);
void add(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret);
void add(const BigInteger& lhs, long, BigInteger& ret);

void subtract(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret, bool neg = false);
void mul(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret);
void mul(const BigInteger& lhs, const BigInteger& rhs, long);
void div(const BigInteger& lhs, const BigInteger& rhs, BigInteger& q, BigInteger& ret, bool byconst = false);
void knuthDiv(const BigInteger& dividend, const BigInteger& divisor, BigInteger& quotient, BigInteger& remainder);
void divByConst(const BigInteger& lhs, const BigInteger& rhs, BigInteger& q, BigInteger& ret);
void _add(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret);
void _subtract(const BigInteger& lhs, const BigInteger& rhs, BigInteger& ret);
BigInteger fromDec(const string& str,size_t);
BigInteger fromBinary(const string& str, Radix radix);


inline unsigned long long fast_atoi(const char* str, unsigned int len);
inline unsigned long fast_btoi(const char* str, unsigned int len);
inline unsigned long fast_htoi(const char* str, unsigned int len);
inline string convert_to_binary(string::const_iterator _s, size_t);
inline string convert_to_binary(string s);
inline string convert_to_binary(unsigned long long, size_t);
inline void get_type(const string& str, Radix& radix, bool& negative, int& beg);

#endif