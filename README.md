# Efficient-BigInteger
A efficient libarary for BigInteger.
这个项目用来比较高精度大数运算的一些算法的优劣,在语言层面上进行了力所能及的优化,并编写了一定的测试用例确保正确性。
整个项目脉络清晰,代码简洁。
我也会给出自己设计时遇到的障碍和解决的方法。

目前进展:已经实现10进制非递归/递归的FFT高精度乘法,KR乘法(用于计算位数中等的数据),长除法,加法,取模,减法,取反,移位,到字符串的相互转换
实现了各种构造函数和内存管理。
对比了用string作为载体，正序表示大整数的方法。

接下来要将整个项目进行一次改动,变成BASE 2的表示办法.这样可以方便的实现进制转换以及移位操作的实现(10进制下效率低),另外由于牛顿迭代的变种实现整除(基于2^k)需要用到大量的移位,如果继续基于10进制则性能很难得到提升。同时加减也会得到很大提高,因为%可以变成&(2^k - 1),进位可以变成移位！

记录一个坑:long long在两个unsigned long乘的时候溢出,开始没有发现,因为在10进制下没有支持二进制,所以就用int作为元素。现在base为2^32的时候就必须考虑无符号的问题了...调了一下午才找到问题:用 unsigned long long作强制转换就好了。

记录一个瓶颈:采用二进制以后,由于一般用户使用十进制,那么就需要进制转换,现在的10进制转二进制算法在3000个10进制的情况下大概20ms,但是如果增加到10000位就要1秒左右了,100000位10进制就几乎算不出来了。。。当然用户输入16进制的话速度很快。

进展:昨天折腾了一天,尝试了大量的方法,现在只能做到对50位万的10进制在30s左右转换完,1万位的大概在40-50ms左右,比昨天快了几十倍,但是还是比不上java的BigInteger(50万位)只要15s,而BigInt这个库函数只要4s就可以搞定,不过我还有很多地方可以优化,现在没有实现+=,* = 等运算符,所以会浪费一次拷贝。另外可以考虑把小数乘法转换成大数乘法。Go On!!!!!!!!!!!!!!.

进展:现在测试10万位加法进行10万次,一共耗6400ms左右。对比 java 的BigInteger 6.394s,我的代码没有太多安全性上的考虑,即便在这样的情况下速度还是没有优势,还得继续加油!
进展:测试10万位减法10万次,一共大约7000ms左右,比java BigInteger稍微i慢一些。如果**类型不一致导致的类型提升**会造成一定影响。所以要尽量避免。
但是**BigInt**这个库真的太快了...比我的快了几乎一倍。。
