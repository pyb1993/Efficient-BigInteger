#include "StringBigdata.h"
#include "BigInteger.h"
#include "iostream"
#include "sstream"
#include "time.h"
#include "cassert"

#define COMPARE 1
#define COMPUTE(x) do{clock_t start_time = clock(); \
	x; \
	clock_t end_time = clock(); \
	cout << "Running time is: " << static_cast<double>\
	(end_time - start_time) / CLOCKS_PER_SEC * 1000 << "ms" << endl; }while (0)  //输出运行时间

void Add_Test();
void Sub_Test();
void MultiPly_Test();
void Mod_Test();
void Power_Mod_Test();
void Test(){
	//Add_Test();
	Sub_Test();
	//MultiPly_Test();
	//Mod_Test();
	Power_Mod_Test();
}

void Add_Test(){
	string a6 = "366666666", b6 = "666666666";
	string a7 = "244444444", b7 = "444444444";
	//COMPUTE(cout << add(a6, b6) << endl);//1033333332
	//COMPUTE(cout << add(a7, b7) << endl);//

	if (COMPARE){
		string a8 = a6 + string(100758900, '5');
		string b8 = b6 + string(60054340, '7') + "437298423097458942307592308";
		COMPUTE(add(a8, b8));
	}

}

void Sub_Test(){
	string a1 = "5422296287557037040", b1 = "71185184863703704", c1 = "2962962958814814816";
//	COMPUTE(cout << Subtract(Subtract(a1, c1), b1) << endl);

//	COMPUTE(cout << Subtract(b1, b1) << endl);

	string a8 = string(100758900, '5');
	string b8 = string(60054340, '7') + "437298423097458942307592308";
	if (COMPARE)//比较性能
    	COMPUTE(Subtract(a8, b8));

}

void MultiPly_Test(){
	string a = string(10000, '3') + "172222222234566666666663456564365835634144443";// 1000 3
	string b = string(8000, '4') + "123455555555555255555555555552454232348897821";// 200 4
	string a2 = string(9000, '3') + "12387766544329887776544345678899";
	string b2 = string(9000, '2') + "54978723598347523985439875923085";
	string a3 = string(90000, '9') + "52345235252";
	string b3 = "532896572389755479872549875423985723987653298";
	string correct;
	//测试正确性
#if 1
	string a4 = "44444444444", b4 = "2222222222";
	correct = "98765432087901234568";
	COMPUTE(cout << (Effective_multiply(a4, b4) == correct) << endl);

	string a5 = "123456789123456789", b5 = "98765432123456789";//
	correct = "12193263126352689864349947750190521";
	COMPUTE(cout << (Effective_multiply(a5, b5) == correct) << endl);


	string a7 = "172222222234566666666663456564365835634144443", b7 = "123455555555555255555555555552454232348897821";
	correct = "21261790124980728703703303694873314334566791251130343554428157608965752534229925361958703";
	COMPUTE(cout << (Effective_multiply(a7, b7) == correct) << endl);

	string a8 = "158610744283305641116339861354533472329775122";
	correct = "25157368202104173113981398220939162641587521514899015555207154711683506072518799090114884";
	COMPUTE(cout << (Effective_multiply(a8, a8) == correct) << endl);

#endif

	//测试性能
# if 0
	COMPUTE(Effective_multiply(a, b));//10000*8000  12275

	COMPUTE(Effective_multiply(a2, b2));// 9000 * 9000  6623

	COMPUTE(cout << Effective_multiply(a3, b3));// 90000 * const  4533
#endif
}

void Mod_Test(){
	string a = "172222222234566666666663456564365835634144443";// 1000 3
	string b = "123455555555555255555555555552454232348897821";// 200 4
	COMPUTE(Mod(a, b));

	string a1 = string(1000, '3') + "158610744283305641116339861354533472329775122";
	string b1 = string(200, '4') + b;
	if (COMPARE)
	    COMPUTE(Mod(a1, b1));


}
void Power_Mod_Test(){

	string correct;
	string a = "172222222234566666666663456564365835634144443";// 1000 3
	string b = "123455555555555255555555555552454232348897821";// 200 4
	string c = "222222222222227777777777777777777777777722571";
	correct = "88130890624492142766705661885585287750833038";
	COMPUTE(cout << (Power_mod(a, b, c) == correct) << endl);// 


	string a1 = string(1000, '3') + a, b1 = string(200, '4') + b;
	//COMPUTE(cout << Power_mod(a1, b1, c) << endl);// 

	string a2 = string(10000, '3') + a, b2 = string(8000, '4') + b;
	COMPUTE(cout << Power_mod(a2, b2, c) << endl);// 
}

void Load_Test(){
#define EXPECT_STR(x) do {BigInteger y = x;assert(y.str() == string(x));}while(0)

	EXPECT_STR("123456789");
	EXPECT_STR("12345678");
	EXPECT_STR("12345678912345678910");
	EXPECT_STR("100000000000000000000");
	EXPECT_STR("-100000000000000000000");
	EXPECT_STR(string(100,'3')+"123456789");
}

void Big_Add_Test(){

#define EXPECT_INT(a,b,c) do {BigInteger x = a;BigInteger y = b;BigInteger z = c;\
	assert(z == (x + y)); }while (0)

	EXPECT_INT("123456", "123456", "246912");
	EXPECT_INT("123456123456", "123456", "123456246912");
	EXPECT_INT("123456123456123456", "123456", "123456123456246912");
	EXPECT_INT("54723895472398547230", "43528945723895723984", "98252841196294271214");
	EXPECT_INT("999999999999999999", "11111111111111", "1000011111111111110");
	if (COMPARE){
		BigInteger a8 = "45239872398472309" + string(100758900, '5');
		BigInteger b8 = "432972384022023949958" + string(60054340, '7') + "437298423097458942307592308";

	COMPUTE(a8 + b8);
}
}
#if 1
void Push_Front_Test(){
#define test_push(_t,x) do{	BigInteger t = _t;	t.append_zeros(x);\
	assert(t.str() == (string(_t)+string(x,'0')));}while (0)
    

	test_push("123456", 5);
	test_push("345232", 100);
	test_push("234523865748937209", 999);
	test_push(string(1000, '4') + "542",1);
	test_push("1111111", 1);
	BigInteger s("0000000000000000");
	s.append_digits("1");
}
#endif

void sub_integer_rough_Test(){
#define test_sub(x,s,e,res) do{BigInteger t = x; \
	auto sub = t.sub_integer_rough(s, e); \
	assert(sub.str() == res); }while (0)

	BigInteger u = "123456789100000000";
	test_sub("123456789100000000", 0, 1,"12");
    test_sub("123456789100000000999999999", 0, 1, "123");
	test_sub("123456789100000000", 2, 1, "0");
	test_sub("123456789100000000", 0, 1, "12");

}

void Operator_Subtract_Test(){

#define EXPECT_INT_SUB(a,b,c) do {BigInteger x = a;BigInteger y = b;BigInteger z = c;\
	assert(z == (x - y)); }while (0)

	
	EXPECT_INT_SUB("123456", "123456", "0");
	EXPECT_INT_SUB("123456123456", "123456", "123456000000");
	EXPECT_INT_SUB("123456123456123456", "123456", "123456123456000000");
	EXPECT_INT_SUB("54723895472398547230", "43528945723895723984", "11194949748502823246");
	EXPECT_INT_SUB("999999999999999999", "11111111111111", "999988888888888888");
	EXPECT_INT_SUB("1", "123456789987654321", "-123456789987654320");
	EXPECT_INT_SUB("5422296287557037040", "71185184863703704", "5351111102693333336");
	BigInteger a8 = "45239872398472309" + string(100758900, '5');
	BigInteger b8 = "432972384022023949958" + string(60054340, '7') + "437298423097458942307592308";
	if (COMPARE){
		COMPUTE(a8 + b8);
	}

}

void Digits_Test(){
#define test_digits(x,s,e,res) do{BigInteger t = x; \
	auto sub = t.sub_integer(s, e); \
	assert(sub.str() == res); }while (0)

	test_digits("123456", 2, 3, "345");
	test_digits("1234567891234567", 0, 9, "123456789");
	test_digits("123456999999999999999",7 , -1, "99999999999999");
	
}
void Operator_Mod_Test(){

#define check_mod(x,y,z) do{assert((x%y)==z);}while(0)
	BigInteger a = "172222222234566666666663456564365835634144443";// 1000 3
	BigInteger b = "123455555555555255555555555552454232348897821";// 200 4
	BigInteger correct = "48766666679011411111107901011911603285246622";
	check_mod(a,b,correct);

	
#if 1
	BigInteger a1 = string(1000, '3') + "158610744283305641116339861354533472329775122";
	BigInteger b1 = (string(200, '4') + "123455555555555255555555555552454232348897821");
	correct = "90691965790704901677330261963876131847549\
58818036214446122552271032357689245532028384930706999117863789625062896214097\
065882559894377367640182030209506368739503353301237508903141281380208333333333\
333158610744283305641116339861354533472329775122";
	check_mod(a1, b1, correct);
	if (COMPARE)
	COMPUTE(a1 % b1);
#endif

	check_mod(BigInteger("20000001122334455"), BigInteger("1000000000"), BigInteger("122334455"));
}

int  main()
{
	/**String Based Test Case*/
	//Test();
	BigInteger x;
	Load_Test();
	Add_Test();
	Big_Add_Test();
	Push_Front_Test();
	Digits_Test();
	sub_integer_rough_Test();
	Sub_Test();
	Operator_Subtract_Test();
	Mod_Test();
	Operator_Mod_Test();
	return 0;
}