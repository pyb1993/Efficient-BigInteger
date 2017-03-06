#include "iostream"
#include "sstream"
#include "time.h"
#include "bitset"
#include "cassert"
#include "StringBigdata.h"
const unsigned long long LMAX = (1 << 64) - 1;
const unsigned long long Second_MAX =( LMAX / 10) - 10;
const string Second_MAX_STR = to_string(Second_MAX);
using namespace std;

#define string_to_int(x) (stoi(x))
#define int_to_string(s) (to_string(s))
#define CHOOSE 1


inline string remove_pre_zero(const string& a)
{
	auto t = a.find_first_not_of('0', 0);//去掉前置的0
	if (t == a.npos)
		return "0";
	else
		return a.substr(t);
}

int Mod(string& a, int & b)
{
	//遍历大数a的每一位不断求余更新
	int t = 0;
	for (auto& x:a)
	t = ((t * 10)%b + x - '0' ) % b;
	return t;
}
string Subtract(const string& a, const string& b,bool reverse)//暂时只处理两个正数
{
	//首先消除前置的0;
	string res;
	auto len_a = a.size();
	auto len_b = b.size();

	if (len_a < len_b||(len_a == len_b && a < b)){
		return Subtract(b,a,true);
	}

	res.resize(len_a+1);
	char count = 0;
	int tmp = 0;
	int pa , pb;
	size_t non_zero_pos = -1;
	
	for (pa = len_a - 1, pb = len_b - 1 ; pa >= 0; --pa, --pb)//计算和b相减
	{	
		char x = (pb >= 0) ? b[pb] : '0';
		tmp = a[pa] - x + count;
		count = 0;
		if (tmp < 0)
		{
			tmp += 10;
			count = -1;
		}
		assert(tmp >= 0);
		res[pa+1] = '0' + tmp;//注意这里是pa + 1,因为考虑可能的负号,多留了一个位置
		if (tmp != 0)
			non_zero_pos = pa + 1;
	}
	if (reverse){
		assert(non_zero_pos > 0);
		res[ --non_zero_pos ] = '-';
	}

	if (non_zero_pos == 0)//避免多余的substr
		return res;
	else if (non_zero_pos == -1)//说明a==b
		return "0";
	return res.substr(non_zero_pos);
}

	

	
string Mod(const string& a, const string& b)//实现大数求余的第二种情况,a,b都很大
{
		//默认a,b不存在前置的0
	    auto a_len = a.length(), b_len = b.length();
		if (a_len< b_len || (a_len == b_len && a < b))//a一定小于b
			return a;
		
		bool k = 0;
		string temp = a.substr(0, b_len);
		string old_ans, mod;
		string record[2] = {temp,temp};
		for (int i = b_len; i < a_len + 1; ++i)
		{
			while(1)//不断的试探
			{
				record[k] = Subtract(record[!k], b);// tmp - b
				if (record[k][0] == '-')
					break;
				k = !k;
			}
			if (i!=a_len)
			 record[!k] += a.substr(i, 1);
		}
		return record[!k];
}

string convert_to_bin(const string& _s)
{
	string binary;
	string s(_s);
	while (s != "0")//模拟不断的除2取余
	{
		int t = 0, old_t = 0;
		for (int i = 0; i != s.length(); ++i)//内循环计算一次除法和一次取余运算
		{
			t = (old_t * 10 + s[i] - '0') % 2;//计算余数
			s[i] = (s[i] - '0' + old_t * 10) / 2 + '0';//计算除以2的这一位的商
			old_t = t;//更新余数
		}
		binary += int_to_string(t);
		s = remove_pre_zero(s);//获得一次计算后的商
	}
	return string(binary.rbegin(), binary.rend());
}
int Power_mod(int& x, int& y, int& z)
{
	//(x^y)%z
	int mod = 1;
	int m = x%z;
	mod = ((y&1 )== 1) ? m : 1;
	while(y/=2)//只要y还不为0
	{
		auto temp = y & 1;//获得y最低位
		m = (m*m) % z;
		if (temp)//如果y二进制形式的这一位是1
			mod = (mod*m) % z;
	}
	return mod;
}
string Power_mod( const string& x,const string& y, const string& z)
{
	 //x = remove_pre_zero(x);
	 //y = remove_pre_zero(y);
	 
	//本算法采取了对大数据相乘的优化（n^1.42复杂度,分治求）,但是在快速求幂的
	//过程中耗时还是太多,简单的对大数乘除的优化不可能达到0.1ms的速度求出x^y%z的结果
	//现在大概200ms
	cout << "大型数据快速求幂: (x^y)%z \n";
	cout << "x的位数: " << x.length() << "\n y的位数" << y.length() << "\nz的位数" << z.length()<<endl;
//	cout<<	"\n(" << x << ")^(" << y << ")" << "%("<<z<< ")\n";

	string m = Mod(x, z);
	string y_bin = convert_to_bin(y);//获取y的二进制
	string mod = (y_bin.back() == '1') ? m : "1";
	for (size_t i = y_bin.length() - 2; i != -1; --i)
	{
		m = Mod(Effective_multiply(m, m), z);
		string ans;
		if (y_bin[i] == '1')
		mod = Mod(Effective_multiply(mod, m), z);
	}
	remove_pre_zero(mod);
	return mod;
}
#if CHOOSE
string Effective_multiply(const string& a, const string& b)
{
	//乘法的分治算法
	if (a.length() < b.length())//确保 a >=  b
		return Effective_multiply(string(b), string(a));
	if (a.length() <= 9){
		//cout << int_to_string(stoll(a)*stoll(b)) << " = " << a << "  X " << b <<"integer mul"<< endl;
		return  int_to_string(stoll(a)*stoll(b));
	}

	string s;
	if (b.length() > 4){
		string a1 = a.substr(0, a.length() / 2), a2 = a.substr(a.length() / 2);
		string b1 = b.substr(0, b.length() / 2), b2 = b.substr(b.length() / 2);
		int weight1 = a.length() - (a.length() / 2), weight2 = b.length() - (b.length() / 2);
	    s = Effective_multiply(a1, b1) + string(weight1 + weight2, '0');
		s = add(s, Effective_multiply(a1, b2) + string(weight1, '0'));
		s = add(s, Effective_multiply(a2, b1) + string(weight2, '0'));
		s = add(s, Effective_multiply(a2, b2));
	}
	else {
		string a1 = a.substr(0, a.length() / 2), a2 = a.substr(a.length() / 2);
		int weight1 = a.length() - (a.length() / 2);
	    s =  Effective_multiply(a1, b) + string(weight1, '0');
		auto s2 =  Effective_multiply(a2, b);
		s = add(s, s2);
	}
	return s;
}
#else
string Effective_multiply(const string& a, const string& b)//只需要 3次子递归
{


	//乘法的分治算法
	auto a_len = a.length(), b_len = b.length();
	assert(a_len > 0 && b_len > 0);
	if (a == "0" || b == "0")
		return "0";

	if (a_len < b_len)//确保a大于b
		return Effective_multiply(b, a);

	if (a_len< 9 && b_len < 9){
		   return  int_to_string(stoll(a)*stoll(b));
	}

	string s;
	if (b_len * 2 > a_len){
		size_t weight = a_len >> 1;
		size_t mid1 = a_len - weight, mid2 = b_len - weight;//注意这里的trick

		string a1 = a.substr(0, mid1), a0 = a.substr(mid1);
		string b1 = b.substr(0, mid2), b0 = b.substr(mid2);
	
		string z2 = Effective_multiply(a1, b1);// z2 = a1*b1
		string z0 = Effective_multiply(a0, b0);// z0 = a0*b0

		string middle_res =  Effective_multiply(add(a1, a0), add(b1, b0));//(a1+a0)(b1+b0)
		
		string z1;
		string tmp;
		//assert((middle_res.find('\0') == string::npos) && (z2.find('\0') == string::npos) && (z0.find('\0') == string::npos));

			 z1 = Subtract(Subtract(middle_res, z2), z0);//z1 = middle_res - z2 - z0，
			
		   //assert(z1[0] != '-' || !(cout << " (" <<a1+"+"+a0 << " ) * ("<<b1+"+"+b0 <<") = "<<middle_res<<endl<<\
			//middle_res<<" - "<< z2+"- "+z0 + " :(" +a0+"* "+b0+ ")"<<" = "<< z1<<endl));

		 s = z2 + string(weight<<1, '0');// z2 * 10^(2*mid)
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

string add(const string&  a, const string& b)
{
	auto a_len = a.size(), b_len = b.size();
	assert(!(a_len == 0 && b_len == 0));

	if (a_len < b_len)
		return add(b, a);
	
	string ans;
	size_t c = 0;//进位
	ans.resize(a_len + 1);
	int pos = ans.size() - 1;
	int i, j;

	for ( i = b_len - 1, j = a_len - 1; i >= 0; --i, --j)//循环计算每一位的加法
  {
		char bit_ans = c + (b[i] - '0') + (a[j] - '0');
		c = bit_ans / 10;//得到新进位
		ans[pos--] = bit_ans - c * 10 + '0';//ans.push_back(bit_ans - c * 10 + '0') ;//得到该位的值,char
  }
	for (; j >= 0; --j)
	{
		char bit_ans = a[j] - '0' + c;//得到该位的值,char
		c = bit_ans / 10;
		ans[pos--] = bit_ans - c * 10 + '0';//ans.push_back(bit_ans - c * 10 + ('0' - 0));//得到该位的值,char
	}
	if (c != 0)//处理最后一位进位
		ans[pos--] = '1';
	return ans.substr(pos+1,string::npos);
}





string multiply(const string& x, const string& y)
{
	string ans(x.size() + y.size() + 1, '0');
	string sub_x, sub_y;
	size_t pos_x = 0;//标记子串的位置
	while (pos_x < x.size())
	{
		sub_x = x.substr(pos_x, M);
		size_t pos_y = 0;
		while (pos_y < y.size())
		{
			sub_y = y.substr(pos_y, M);
			size_t weight = x.size() - pos_x - sub_x.size() + y.size() - pos_y - sub_y.size();//定义这一次子串的权重

			long long temp_result = stod(sub_x)*stod(sub_y);
			if (temp_result)
				add(ans, to_string(temp_result) + string(weight, '0'));//与结果累加
			pos_y += M;
		}
		pos_x += M;
	}
	return string(ans.find_first_not_of('0') + ans.begin(), ans.end());
}