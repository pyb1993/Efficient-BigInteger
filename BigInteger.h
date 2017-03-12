#ifndef BIGINTEGERPYB
#define BIGINTEGERPYB

#include "iostream"
#include "string"
#include "vector"
#include "algorithm"
#include "memory"
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

	friend class FFTMultiPlier;

	friend inline void swap(BigInteger& lhs, BigInteger& rhs);
	friend ostream& operator << (ostream& out, const BigInteger& big);
	friend BigInteger operator -(const BigInteger& big);
	friend BigInteger subtract  (const BigInteger& lhs, const BigInteger& r,  bool neg = false);

	BigInteger();
	BigInteger(const BigInteger&);
	BigInteger( unsigned long long& num);//可能溢出,必须用unsigned long long
	BigInteger(const long long& num);

	BigInteger(const string&);
	BigInteger(const vector<double>&, size_t len, short other_width = 8);//用来接收一个数组,专门处理FFT的结果
	BigInteger(BigInteger&&) NOEXCEPT;
	~BigInteger();


	BigInteger& operator= (const BigInteger& rhs);
	BigInteger& operator= (long long num);
	BigInteger& operator= (unsigned long long& num);

	BigInteger& operator= (const string& str);
	BigInteger& operator= (BigInteger&&) NOEXCEPT;



	BigInteger operator + (const BigInteger&) const;
	BigInteger operator + (const long long &) const;

	BigInteger operator - (const BigInteger&) const;
	BigInteger operator * (const BigInteger&) const;
	BigInteger operator * (const long &) const;

	BigInteger operator / (const BigInteger&) const;
	BigInteger operator %(const BigInteger&) const;
	BigInteger operator <<(size_t num) const;
	BigInteger operator >>(size_t num) const;

	void operator *= (const BigInteger&);//此处是否应该返回
	void operator += (const BigInteger&);




	bool operator < (const BigInteger& other)const;
	bool operator > (const BigInteger& other)const;
	bool operator >= (const BigInteger& other)const;

	bool operator== (const BigInteger&) const;
	eletype& operator[] (size_t index);
	eletype operator[] (size_t index) const;
	size_t length() const;
	size_t digits_len()const;
	void resize(size_t size);
	void reserve(size_t size);
	void push_back(const eletype& x);
	string str() const;

	BigInteger sub_integer_rough(size_t start , unsigned int len = npos ) const;
	BigInteger sub_integer(size_t start, unsigned int   len = npos) const;
	BigInteger Inverse() const;
	void set_resversable();
	unsigned short highest_bitsnum() const;
	operator string () const;
    
private:
	void parse_binary(const string& str,Radix radix);
	string str_from_vec() const;
	inline void remove_redundant_zeros();
	datatype* s = nullptr;
	bool negative = false;
	static const size_t npos = -1;
	
};
inline void swap(BigInteger& lhs, BigInteger& rhs);
ostream& operator << (ostream& out, const BigInteger& big);
BigInteger subtract(const BigInteger& lhs, const BigInteger& r, \
const vector<BigInteger::eletype>&, bool neg = false);
inline unsigned long long fast_atoi(const char* str, int len);
inline string convert_to_binary(string::const_iterator _s, size_t len);
inline int get_min_binary_power(int x);
inline void bin_shift_multiply( BigInteger&, const string&);
inline void get_type(const string& str, Radix& radix, bool& negative, int& beg);
#endif