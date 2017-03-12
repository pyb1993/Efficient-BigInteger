#include "StringBigdata.h"
#include "BigInteger.h"
#include "iostream"
#include "sstream"
#include "time.h"
#include "cassert"
#include "FFTFunctor.h"
#include "FFTMultiplier.h"
#define COMPARE 1
#define COMPUTE(x) do{clock_t start_time = clock(); \
	x; \
	clock_t end_time = clock(); \
	cout << "Running time is: " << static_cast<double>\
	(end_time - start_time) / CLOCKS_PER_SEC * 1000 << "ms" << endl; }while (0)  //输出运行时间
#define EXPECT_INT(result,correct) do{assert(result == correct);}while(0)

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

#if 0
	if (COMPARE){
		string a8 = a6 + string(100758900, '5');
		string b8 = b6 + string(60054340, '7') + "437298423097458942307592308";
		COMPUTE(add(a8, b8));
	}
#endif

}

void Sub_Test(){
	string a1 = "5422296287557037040", b1 = "71185184863703704", c1 = "2962962958814814816";
//	COMPUTE(cout << Subtract(Subtract(a1, c1), b1) << endl);

//	COMPUTE(cout << Subtract(b1, b1) << endl);

	
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
	//COMPUTE(cout << (Effective_multiply(a4, b4) == correct) << endl);

	string a5 = "123456789123456789", b5 = "98765432123456789";//
	correct = "12193263126352689864349947750190521";
	//COMPUTE(cout << (Effective_multiply(a5, b5) == correct) << endl);


	string a7 = "172222222234566666666663456564365835634144443", b7 = "123455555555555255555555555552454232348897821";
	correct = "21261790124980728703703303694873314334566791251130343554428157608965752534229925361958703";
	//COMPUTE(cout << (Effective_multiply(a7, b7) == correct) << endl);

	string a8 = "158610744283305641116339861354533472329775122";
	correct = "25157368202104173113981398220939162641587521514899015555207154711683506072518799090114884";
	//COMPUTE(cout << (Effective_multiply(a8, a8) == correct) << endl);

#endif

	//测试性能
	if (COMPARE){
		COMPUTE(Effective_multiply(a, b));//10000*8000  12275
		COMPUTE(Effective_multiply(a2, b2));// 9000 * 9000  6623
		COMPUTE(Effective_multiply(a3, b3));// 90000 * const  4533
	}
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
	//COMPUTE(cout << (Power_mod(a, b, c) == correct) << endl);// 


	string a1 = string(1000, '3') + a, b1 = string(200, '4') + b;
	//COMPUTE(cout << Power_mod(a1, b1, c) << endl);// 

	string a2 = string(10000, '3') + a, b2 = string(8000, '4') + b;
	//COMPUTE(cout << Power_mod(a2, b2, c) << endl);// 
}

void Load_Test(){
#define EXPECT_STR(x) do {BigInteger y = x;assert(y.str() == string(x));}while(0)
#define EXPECT_BIN_STR(x,c) do {BigInteger y = x;assert(y.str()==c);}while(0)
#if 1

	/*****测试10进制***********************************/
	EXPECT_STR("123456789");
	EXPECT_STR("12345678");
	EXPECT_STR("12345678912345678910");
	EXPECT_STR("100000000000000000000");
	EXPECT_STR("-100000000000000000000");
	EXPECT_STR(string(100,'3')+"123456789");
	EXPECT_STR("100000000000000000000000000000000000000000000000000000000");
	EXPECT_STR(string(6, '1'));
	EXPECT_STR(string(11, '1'));
	EXPECT_STR(string(1023, '3'));
	EXPECT_STR("-123456789987654320");
	EXPECT_STR((string(6888, '1') + "523452345"));
	EXPECT_STR((string(11888, '5') + "523537289452345"));

	/*****测试二进制*******************************************/
	//cout << BigInteger("0b1111100001111111000011111111111111111111");
	EXPECT_BIN_STR("0b1111100001111111000011111111111111111111","1067283644415");
	EXPECT_BIN_STR("0b11110000110101010101001010101001010101010","2068742230698");
	EXPECT_BIN_STR("0b1111000011010101010100","3945812");
	/*****测试16进制**********************/
	EXPECT_BIN_STR("0xABCDE2203163623452454", "12981160722937192563156052");
	EXPECT_BIN_STR("0xABCDEF758329057582203163623452454", "3653876239062442701971792070265375368276");
	EXPECT_BIN_STR("0x00000", "0");
	EXPECT_BIN_STR("0xF", "15");
#endif

	if (!COMPARE)
		return;

	do {
		cout << "test parsing decimal string\n";

		/***********测试10进制效率**********************/
		string a = "12434213431" + string(10000, '3');
		string a2 = "12434213431" + string(50000, '3');
		string a3 = "12434213431" + string(100000, '3');
		string a4 = "12434213431" + string(500000, '3');

		BigInteger b;

		COMPUTE(b = BigInteger(a));
		COMPUTE(b = BigInteger(a2));
		COMPUTE(b = BigInteger(a3));
		COMPUTE(b = BigInteger(a4););
		cout << b.length() << endl;
	} while (0);
}

void Big_Add_Test(){

#define EXPECT_BIG_ADD(a,b,c) do {BigInteger x = a;BigInteger y = b;BigInteger z = c;\
	EXPECT_INT(x + y, z); }while (0)
#define  EXPECT_INT_ADD(a,b,c) do {EXPECT_INT(BigInteger(a)+b,BigInteger(c));}while(0)

	/*测试与long运算*/
	EXPECT_INT_ADD("123456789011111111", 1,"123456789011111112");
	EXPECT_INT_ADD("123456789054234511111111", 1342424431, "123456789054235853535542");
	EXPECT_INT_ADD("999999999999999999999999999", ((1 << 31) - 1), "1000000000000000002147483646");

	/*测试Big+*/

	EXPECT_BIG_ADD("123456", "123456", "246912");
	EXPECT_BIG_ADD("123456123456", "123456", "123456246912");
	EXPECT_BIG_ADD("123456123456123456", "123456", "123456123456246912");
	EXPECT_BIG_ADD("54723895472398547230", "43528945723895723984", "98252841196294271214");
	EXPECT_BIG_ADD("999999999999999999", "11111111111111", "1000011111111111110");
	EXPECT_BIG_ADD("15681","7442833","7458514");
	EXPECT_BIG_ADD("95437298723452938475389443298903485909254234524575382543245689999999999999999999122", \
				   "75432729875438254372894572389657252345234523453987582395743987583429754298574442833",\
		"170870028598891192848284015688560738254488757978562964938989677583429754298574441955");

	
		if (COMPARE){
			do{
				cout << "test BigInteger operator + :\n";
				BigInteger a = "0x45239872398472309" + string(100758900, '5');
				BigInteger b = "0x432972384022023949958" + string(60054340, '7') + "437298423097458942307592308";
				BigInteger a1 = "0x45239872398472309" + string(100000, '5');
				BigInteger b1 = "0x432972384022023949958" + string(100000, '7') + "437298423097458942307592308";

				COMPUTE(a + b);
				COMPUTE(for (int i = 0; i < 100000; ++i) auto c = a1 + b1;);
			} while (0);


			do{
				cout << "test string add  :" << endl;
				string a = string(100758900, '5');
				string b = string(60054340, '7') + "437298423097458942307592308";
				BigInteger a1 = "45239872398472309" + string(100000, '5');
				BigInteger b1 = "432972384022023949958" + string(100000, '7') + "437298423097458942307592308";

				COMPUTE(add(a, b));
				//COMPUTE(for (int i = 0; i < 100000; ++i) add(a1,b1);); too slow

			}while(0);
		}
	
}



void sub_integer_rough_Test(){
#define test_sub(x,s,e,res) do{BigInteger t = x; \
	auto sub = t.sub_integer_rough(s, e); \
	assert(sub.str() == res); }while (0)

	test_sub("123456789100000000", 0, 1,"28744523");
    test_sub("123456789100000000999999999", 0, 1, "6692605");
	//test_sub("1234567891000000000000000000000000000000", 4, 1, "0");
	//test_sub("123456789100000000", 0, 1, "12");

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
	EXPECT_INT_SUB("71185184863703704", "70185184863745604", "999999999958100");

	if (COMPARE)//比较性能
	{
		do{
			cout << "test operator -\n";
			BigInteger a = "0x4523" + string(100758900, '5');
			BigInteger b = "0x5234532" + string(60054340, '7') + "437298423097458942307592308";
			BigInteger a1 = "0x4523532452" + string(100030, '5');
			BigInteger b1 = "0x5234532" + string(100000, '7') + "437298423097458942307592308";
			COMPUTE(auto c = a - b;);
			COMPUTE(for (int i = 0; i < 100000;++i) auto c = a1 - b1;);

		} while (0);


		do{
			cout << "test string subtract\n";
			string a = string(100758900, '5');
			string b = string(60054340, '7') + "437298423097458942307592308";
			COMPUTE(Subtract(a, b));
		} while (0);

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

void Operator_Mul_Test(){
#define EXPECT_BIG_MUL(a,b,c) do {EXPECT_INT((BigInteger(a)*BigInteger(b)),BigInteger(c));} while(0)
#define EXPECT_INT_MUL(a,b,c) do {EXPECT_INT((BigInteger(a)* b ),BigInteger(c));}while(0)


	/*BIG + INT */
	EXPECT_INT_MUL("12345578329532984572398", 1, "12345578329532984572398");
	EXPECT_INT_MUL("12345", 12345, "152399025");
	EXPECT_INT_MUL("7458514", 7458514, "55629431088196");
	EXPECT_INT_MUL("123456789123456789", 2147483647, \
		"265121435753750918488629483");

	/*BIG + BIG*/
	EXPECT_BIG_MUL("1", "12345578329532984572398", "12345578329532984572398");
	EXPECT_BIG_MUL("34252354231", "0", "0");
	EXPECT_BIG_MUL("12345", "12345", "152399025");
	EXPECT_BIG_MUL("7458514", "7458514", "55629431088196");
	EXPECT_BIG_MUL("12345987654321", "12345543212345", "152417924085497789759792745");
	EXPECT_BIG_MUL("44444444444", "2222222222", "98765432087901234568");
	EXPECT_BIG_MUL("123456789123456789", "98765432123456789", \
		"12193263126352689864349947750190521");

	EXPECT_BIG_MUL("172222222234566666666663456564365835634144443", \
		"123455555555555255555555555552454232348897821", \
		"21261790124980728703703303694873314334566791251130343554428157608965752534229925361958703");

	string a8 = "158610744283305641116339861354533472329775122";
	EXPECT_BIG_MUL(a8, a8, \
		"25157368202104173113981398220939162641587521514899015555207154711683506072518799090114884");

	if (COMPARE){
		do{
			cout << "test Big operator mul:\n";
			BigInteger a = "0x532"+string(50000, '3') + "172222222234566666666663456564365835634144443";// 1000 3
			BigInteger b = "0x6263"+string(50000, '4') + "123455555555555255555555555552454232348897821";// 200 4
			BigInteger a2 = "0x532"+string(9000, '3') + "12387766544329887776544345678899";
			BigInteger b2 = "0x535423A2"+string(9000, '2') + "54978723598347523985439875923085";
			BigInteger a3 = "0x532D"+string(90000, '9') + "52345235252";
			BigInteger b3 = "0x532896572389755479872549875423985723987653298";

			COMPUTE(a = a*b);// 50000 * 45000
			COMPUTE(a2 = a2*b2);// 9000 * 9000
			COMPUTE(a3 = a3*b3);// 90000 * const

		

		} while (0);

	 }
}

void Operator_Shif_Test(){
#define EXPECT_INT_SHIFT(B,num,C) do{\
	auto _y = BigInteger(B)<<num;EXPECT_INT(_y,BigInteger(C));}while(0)

	EXPECT_INT_SHIFT("1", 2, "4");
	EXPECT_INT_SHIFT("1111222212", 1, "2222444424");
	EXPECT_INT_SHIFT("0", 2, "0");
	EXPECT_INT_SHIFT("15342424532453254", 1, "30684849064906508");
	EXPECT_INT_SHIFT("15342424532453254", 10, "15710642721232132096");
	EXPECT_INT_SHIFT("15342424532453254", 12, "62842570884928528384");
	EXPECT_INT_SHIFT("15342424532453254", 50, "17274034351829107757228014698496");
	EXPECT_INT_SHIFT("-1", 1, "-2");
}

void FFT_Test(){
#define EXPECT_FFT_MUL(x,y,correct) do{\
	FFTMultiPlier f(BigInteger(x), BigInteger(y)); \
	assert(f.MulWithFFT() == BigInteger(correct));}while (0)

	EXPECT_FFT_MUL("1", "4294967296","4294967296");
	EXPECT_FFT_MUL("123456", "123456", "15241383936");
	EXPECT_FFT_MUL("123456789", "123456789", "15241578750190521");
	EXPECT_FFT_MUL("1", "12345578329532984572398", "12345578329532984572398");

	auto f = FFTMultiPlier(BigInteger("0"), BigInteger("34252354231"));
	auto ans =  (f.MulWithFFT());

	EXPECT_FFT_MUL("34252354231", "0", "0");
	EXPECT_FFT_MUL("12345", "12345", "152399025");
	EXPECT_FFT_MUL("7458514", "7458514", "55629431088196");
	EXPECT_FFT_MUL("12345987654321", "12345543212345", "152417924085497789759792745");
	EXPECT_FFT_MUL("44444444444", "2222222222", "98765432087901234568");
	EXPECT_FFT_MUL("123456789123456789", "98765432123456789", \
		"12193263126352689864349947750190521");
	EXPECT_FFT_MUL("172222222234566666666663456564365835634144443", \
		"123455555555555255555555555552454232348897821", \
		"21261790124980728703703303694873314334566791251130343554428157608965752534229925361958703");

	string a3 = "158610744283305641116339861354533472329775122";
	EXPECT_FFT_MUL("158610744283305641116339861354533472329775122", a3, \
		"25157368202104173113981398220939162641587521514899015555207154711683506072518799090114884");
	BigInteger a4 = "1" + string(2000, '0');
	EXPECT_FFT_MUL(a4, ("3" + string(2000, '0')), ("3" + string(4000, '0')));

	if (COMPARE){
		do{
			cout << "test FFT operator for big mul:\n";
			BigInteger a = "0x532" + string(500000, '3') + "172222222234566666666663456564365835634144443";// 1000 3
			BigInteger b = "0x6263" + string(500000, '4') + "123455555555555255555555555552454232348897821";// 200 4
		
			BigInteger c(0);
			COMPUTE(auto f = FFTMultiPlier(a, b);   c = f.MulWithFFT(););// 50000 * 45000
			cout << c.length();
		} while (0);

		do{
			cout << "test Big operator And FFT for small mul:\n";
			BigInteger a0 = "0x532" + string(1, '3') + "12387766544329887776544345678899";
			BigInteger b0 = "0x535423A2" + string(1, '2') + "54978723598347523985439875923085";

			BigInteger a1 = "0x532" + string(1000, '3') + "12387766544329887776544345678899";
			BigInteger b1 = "0x535423A2" + string(1000, '2') + "54978723598347523985439875923085";

			BigInteger a2 = "0x532" + string(5000, '3') + "172222222234566666666663456564365835634144443";// 1000 3
			BigInteger b2 = "0x6263" + string(5000, '4') + "123455555555555255555555555552454232348897821";// 200 4
			
			BigInteger a3 = "0x532D" + string(50000, '9') + "52345235252";
			BigInteger b3 = "0x532896572389755479872549875423985723987653298";

			BigInteger a4 = "0x532D" + string(150000, '9') + "52345235252";
			BigInteger b4 = "0x532896572389755" + string(500000,'7');

			cout << "Big Operator *\n";
			COMPUTE(a0*b0);// 100 100
			COMPUTE(a1*b1);// 1000 * 1000
			COMPUTE(a2*b2);// 5000 * 5000
			COMPUTE(a3*b3);// 50000 * 50000
			COMPUTE(a4*b4);//100000 * const

			cout << "FFT  *\n";
			BigInteger c(0);
			COMPUTE(auto f = FFTMultiPlier(a0, b0); c = f.MulWithFFT());// 100 100
			COMPUTE(auto f = FFTMultiPlier(a1, b1); c = f.MulWithFFT());// 100 100
			COMPUTE(auto f = FFTMultiPlier(a2, b2); c = f.MulWithFFT());// 100 100
			COMPUTE(auto f = FFTMultiPlier(a3, b3); c = f.MulWithFFT());// 100 100
			COMPUTE(auto f = FFTMultiPlier(a4, b4); c = f.MulWithFFT());// 100 100

		
		} while (0);
	}

}

void Operator_addeq_Test(){

#define EXPECT_BIG_ADD(a,b,c) do {BigInteger x = a;BigInteger y = b;BigInteger z = c;\
	EXPECT_INT(x + y, z); }while (0)

#define EXPECT_BIG_ADDEQ(a,b,c) do{BigInteger sum = a;sum += BigInteger(b);\
	EXPECT_INT(sum, BigInteger(c));} while (0)

	EXPECT_BIG_ADDEQ("123456", "123456", "246912");
	EXPECT_BIG_ADDEQ("123456123456", "123456", "123456246912");
	EXPECT_BIG_ADDEQ("123456123456123456", "123456", "123456123456246912");
	EXPECT_BIG_ADDEQ("54723895472398547230", "43528945723895723984", "98252841196294271214");
	EXPECT_BIG_ADDEQ("999999999999999999", "11111111111111", "1000011111111111110");
	EXPECT_BIG_ADDEQ("15681", "7442833", "7458514");
	EXPECT_BIG_ADDEQ("95437298723452938475389443298903485909254234524575382543245689999999999999999999122", \
		"75432729875438254372894572389657252345234523453987582395743987583429754298574442833", \
		"170870028598891192848284015688560738254488757978562964938989677583429754298574441955");

	if (COMPARE){
		do{
			cout << "test BigInteger operator + :\n";
			BigInteger _a = "0x45239872398472309" + string(100758900, '5');
			BigInteger _b = "0x432972384022023949958" + string(60054340, '7') + "437298423097458942307592308";
			BigInteger _a1 = "0x45239872398472309" + string(100000, '5');
			BigInteger _b1 = "0x432972384022023949958" + string(100000, '7') + "437298423097458942307592308";

			auto a = _a, b = _b, a1 = _a1, b1 = _b1;
			COMPUTE(a = a + b);
			COMPUTE(for (int i = 0; i < 100000; ++i) {a1 = a1 + b1;});

			a = _a;b = _b;a1 = _a1;	b1 = _b1;
			cout << "compare with operator +=\n";
			COMPUTE(a += b);
			COMPUTE(for (int i = 0; i < 100000; ++i) { a1 += b1;});

		} while (0);
	}
	}


int  main()
{
	/**String Based Test Case*/
	//Test();

	Load_Test();

#if 0

	Digits_Test();
	Push_Front_Test();
	Mod_Test();
	Sub_Test();

	Operator_Mod_Test();
	sub_integer_rough_Test();
	Operator_Shif_Test();
	Big_Add_Test();
	Operator_Subtract_Test();
#endif
	//MultiPly_Test();
	//Operator_Mul_Test();
	//FFT_Test();
	//Operator_addeq_Test();
	return 0;
}