#ifndef BIGINTEGERPYB
#define BIGINTEGERPYB

#include "iostream"
#include "string"
#include "vector"
#include "algorithm"
#include "memory"

#ifndef _MSC_VER
#define NOEXCEPT noexcept
#else
#define NOEXCEPT
#endif

using namespace std;
class BigInteger{	
public:
	using eletype = int;
	using datatype = vector<eletype>;
	static const int DefaultBase = 100000000;
	static const int DefaultWIDTH = 8;
	const int BASE;
	const int WIDTH;

	friend class FFTMultiPlier;

	friend inline void swap(BigInteger& lhs, BigInteger& rhs);
	friend ostream& operator << (ostream& out, const BigInteger& big);
	friend BigInteger operator - (const BigInteger& big);
	friend BigInteger subtract     (const BigInteger& lhs, const BigInteger& r,  bool neg = false);

	BigInteger(const int Base = DefaultBase, const int WIDTH = DefaultWIDTH);
	BigInteger(const BigInteger&, const int Base = DefaultBase, const int WIDTH = DefaultWIDTH);
	BigInteger(const long long& num, const int Base = DefaultBase, const int WIDTH = DefaultWIDTH);
	BigInteger(const string&, const int Base = DefaultBase, const int WIDTH = DefaultWIDTH);
	BigInteger(const vector<double>&, size_t len, long long other_base = 10000,short other_width = 4);//用来接收一个数组,专门处理FFT的结果
	BigInteger(BigInteger&&, const int Base = DefaultBase, const int WIDTH = DefaultWIDTH) NOEXCEPT;
	~BigInteger();


	BigInteger& operator= (const BigInteger& rhs);
	BigInteger& operator= ( long long num);
	BigInteger& operator= (const string& str);
	BigInteger& operator= (BigInteger&&) NOEXCEPT;



	BigInteger operator + (const BigInteger&) const;
	BigInteger operator - (const BigInteger&) const;
	BigInteger operator * (const BigInteger&) const;
	BigInteger operator / (const BigInteger&) const;
	BigInteger operator %(const BigInteger&) const;
	BigInteger operator <<(size_t num) const;
	BigInteger operator >>(size_t num) const;


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
	void append_zeros(size_t);
	void append_digits(const string&);
	string         str() const;

	BigInteger sub_integer_rough(size_t start , unsigned int len = npos ) const;
	BigInteger sub_integer(size_t start, unsigned int   len = npos) const;
	BigInteger Inverse() const;
	void set_resversable();
	operator string () const;
    
private:

	void vec_from_str(const string&);
	string str_from_vec() const;
	void remove_redundant_zeros();
	datatype* s = nullptr;
	bool negative = false;
	static const size_t npos = -1;
	
};
inline void swap(BigInteger& lhs, BigInteger& rhs);
ostream& operator << (ostream& out, const BigInteger& big);
BigInteger operator - (const BigInteger&);
BigInteger subtract(const BigInteger& lhs, const BigInteger& r, \
const vector<BigInteger::eletype>&, bool neg = false);
inline  int fast_atoi(const char* str, int len);
#endif