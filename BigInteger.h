#ifndef BIGINTEGERPYB
#define BIGINTEGERPYB

#include "iostream"
#include "string"
#include "vector"
#include "algorithm"
#include "memory"
using namespace std;
class BigInteger{	
public:
	static const int BASE = 100000000;
	static const int WIDTH = 8;
	using eletype = int;
	using datatype = vector<eletype>;


	friend inline void swap(BigInteger& lhs, BigInteger& rhs);
	friend ostream& operator << (ostream& out, const BigInteger& big);
	friend BigInteger operator - (const BigInteger& big);
	friend BigInteger subtract     (const BigInteger& lhs, const BigInteger& r,  bool neg = false);

	BigInteger();
	BigInteger(const BigInteger&);
	BigInteger (const long long& num);
	BigInteger(const string&);
	~BigInteger();


	BigInteger& operator= (const BigInteger& rhs);
	BigInteger& operator= ( long long num);
	BigInteger& operator= (const string& str);
	BigInteger operator + (const BigInteger&) const;
	BigInteger operator - (const BigInteger&) const;
	BigInteger operator * (const BigInteger&) const;
	BigInteger operator %(const BigInteger&) const;
	


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

	BigInteger sub_integer_rough(size_t start , size_t len = -1 ) const;
	BigInteger sub_integer(size_t start, size_t   len = -1) const;
	
	void set_resversable();
	operator string () const;
    
private:

	void vec_from_str(const string&);
	 string str_from_vec() const;

	datatype* s = nullptr;
	bool negative = false;
};
inline void swap(BigInteger& lhs, BigInteger& rhs);
ostream& operator << (ostream& out, const BigInteger& big);
BigInteger operator - (const BigInteger&);
BigInteger subtract(const BigInteger& lhs, const BigInteger& r, \
const vector<BigInteger::eletype>&, bool neg = false);
inline  int fast_atoi(const char* str, int len);
#endif