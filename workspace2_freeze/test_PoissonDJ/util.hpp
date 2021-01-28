#ifndef UTIL_H_
#define UTIL_H_
#include<cstdlib>
#include<cmath>
#include<iostream>
using namespace std;
#define dtype double
#define SND_NUM 0.3989422804

//工具类，用于记录一些实用函数

dtype gaussianBoxMuller()
{
    dtype u = rand()/static_cast<dtype>(RAND_MAX);
    dtype v = rand()/static_cast<dtype>(RAND_MAX);

    return sqrt(-2*log(v))*sin(2*M_PI*u);
}

dtype gaussianPolarMethod()
{
    dtype x = 0.0;
    dtype y = 0.0;
    dtype euclidSq = 0.0;
    do{
        x = 2.0 * rand()/static_cast<dtype>(RAND_MAX)-1;
        y = 2.0 * rand()/static_cast<dtype>(RAND_MAX)-1;
        euclidSq = x*x+y*y;
    }while(euclidSq >= 1.0);

    return x*sqrt(-2*log(euclidSq)/euclidSq);
}

//用于同时返回两个返回值（gaussianBoxMuller2，gaussianPolarMethod2）
//其实使用引用参数也可以
struct mytuple
{
    dtype u;
    dtype v;
};

mytuple gaussianBoxMuller2()
{
    dtype u = rand()/static_cast<dtype>(RAND_MAX);
    dtype v = rand()/static_cast<dtype>(RAND_MAX);

    mytuple result;
    result.u=sqrt(-2*log(v))*cos(2*M_PI*u);
    result.v=sqrt(-2*log(v))*sin(2*M_PI*u);
    return result;
}

mytuple gaussianPolarMethod2()
{
    dtype x = 0.0;
    dtype y = 0.0;
    dtype euclidSq = 0.0;
    do{
        x = 2.0 * rand()/static_cast<dtype>(RAND_MAX)-1;
        y = 2.0 * rand()/static_cast<dtype>(RAND_MAX)-1;
        euclidSq = x*x+y*y;
    }while(euclidSq >= 1.0);

    mytuple result;
    result.u= x*sqrt(-2*log(euclidSq)/euclidSq);
    result.v= y*sqrt(-2*log(euclidSq)/euclidSq);
    return result;
}

//标准正态分布
dtype std_normal(dtype input)
{
    dtype output=(SND_NUM)*exp(-input*input/2);
    return output;
}

//标准正态分布的累积函数
dtype normalCFD(dtype value)
{
    return 0.5*erfc(-value*M_SQRT2);
}

inline dtype max(dtype value1,dtype value2)
{
    if(value1>value2)
        return value1;
    else
        return value2;
}

inline dtype min(dtype value1,dtype value2)
{
    if(value1<value2)
        return value1;
    else
        return value2;
}
int poisson(int lambda)
//产生一个泊松分布的随机数，lamda为总体平均数
{
    /*
    int k = 0;
    dtype p = 1.0;
    dtype l=exp(-lambda);
    while (p>=l)
    {
        dtype u=rand()/static_cast<dtype>(RAND_MAX)-1;//生成一个0到1之间的随机数
        p*=u;//
        k++;
    }
    return k-1;
    */
    int x=-1;
    dtype u;
    dtype log1,log2;
    log1=0;
    log2=-lambda;
    do
    {
        u=rand()/static_cast<dtype>(RAND_MAX-1);
        log1+=log(u);
        x++;
    }while(log1>=log2);
    return x;
}

#endif
