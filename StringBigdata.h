#ifndef BIGNUMBER_ADD
#define BIGNUMBER_ADD
#include "iostream"
#include  "string"
#include "algorithm"

using namespace std;

const int M =7;
string add(const string& ans, const string& temp);
string multiply(const string& x,const string& y);
inline int Mod(string&, int&);//首先实现大数求余1 : A%B(A很大B可以用整数表示)
string Mod(const string&, const string&);//实现大数求余2 : A%B(A，B不可以用整数表示)
string Subtract(const  string& a, const string& b, bool reverse = false);//实现大数求减法 A-B
inline string remove_pre_zero(const string& a);
string convert_to_bin(const string&);//实现大数转换到二进制
string Effective_multiply(const string&, const string&);
int      Power_mod(int& x, int& y, int& z);
string Power_mod(const string& x, const string& y, const  string& z);//x,y,z都很大,比如x=5423123456789987654331452345;


#endif
