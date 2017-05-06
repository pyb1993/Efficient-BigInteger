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

进展:测试了一下非递归的FFT,50万位的乘法只需要0.5s,几乎和BigInt持平,比BigInteger快很多倍.但是由于误差,所以还是要配合KR算法突破位数的限制。

进展:提高了10进制到2进制转换的效率,方法是将原来的*10变成大整数乘法,同时利用FFT优化大整数乘法,所以速度变快。现在转换50万位的10进制数只需要25s,当然比起BigInteger还是差不少,更不要说BigInt了。继续加油!

进展:对FFT做了一些小优化,现在乘法进步到4.7s左右了,进制转换进步到21s左右，有了不错的提高.同时针对+=作了特殊优化,原因是很多操作都有a = a+b这样的操作,现在我们直接避免+=操作返回任何东西,而是在a的基础上直接计算,这样避免了一次拷贝,经过测试a+=b比a = a+b提高了40%的速度。但是代价就是代码写的很难看,因为所有的代码都要避免拷贝。

进展：现在parsing部分已经可以达到7.6s的时间,已经达到了和BigInt一样的水平(当然我没有太多工程上的考虑)，主要优化的办法就是递归计算10进制的大整数，最开始取２００００位的大整数，然后我们现在需要计算两个部分，第一10^(20000)的大整数表示以及那个两万位的大整数，然后接下来就开始递归求解20000万位的大整数，注意10^20000可以通过预处理先计算好，递归的时候我们降低一次取的标准，每次取10000位，然后执行同样的操作，最后就可以降低到一个能够直接求解的地步！

接下来主要是修改toString这个模块，现在采取double dabble算法进行转换，对于大数据来说太慢了。我准备使用原来的１０进展大数类来进行转换。
暂时没有去理toString,而是测试了除法，发现20万位/10万位需10s,需要的时间太久了，４０万/２０万位的则很久都计算不出来结果，几种常见的大数类都只要３－５ｓ就能计算出来。

##增加knuth algorithm D; 

最终结果:parsing:50 万位10进制 转 BigInteger  7.6s
        +       :10万位加法 执行 10万次 4s
                 1亿位加法,46ms
                 
        -       :1亿位减法 41ms
                 10万减法执行10万次,6s
        
        *       :50万位乘法,0.46s
                 直接计算50000的阶乘,2.43 s
        
        /        :40万位/20万位 3.9s
        
        toString: 转换为3万位10进制,0.336s,转换为21万位10进制,2.78s
        
 最后发现BigInt之所以很快有部分原因是它很多地方都没有返回最终结果,而是依靠将结果写入参数,这和BigInteger很相似,不能直接用+,/,*,/,%,<<,>>等符号，不是很方便,但是节约了一次拷贝。
