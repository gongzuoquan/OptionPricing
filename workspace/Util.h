#ifndef UTIL_H_
#define UTIL_H_
#include<cstdlib>
#include<cmath>

#define dtype float
#define SND_NUM 0.3989422804
#define MT_N 624
//工具类，用于记录一些实用函数

//MT19937 待验证
class MT19937{ 
    bool isInit;
    unsigned int MT[MT_N];
    int index;
    void srand(int seed)
    {
        index=0;
        isInit=1;
        MT[0]=seed;
        for (int i=1;i<MT_N;i++)
        {
            unsigned int t=1812433253*(MT[i-1])^(MT[i-1]>>30)+i;
            MT[i]=t&0xffffffff;
        }
    }
    void generate()
    {
        for(int i=0;i<MT_N;i++)
        {
            unsigned int y=(MT[i]&0x80000000)+(MT[(i+1)%MT_N]&0x7fffffff);
            MT[i]=MT[(i+397)%MT_N]^(y>>1);
            if(y&1)
                MT[i]^=2567483615;
        }
    }
    unsigned int Rand()
    {
        if (!isInit)
            srand((int)time(NULL));
        if (index==0)
            generate();
        unsigned int y=MT[index];
        y=y^(y>>11);
        y=y^((y<<7)&2636928640);
        y=y^((y<<15)&4022730752);
        y=y^(y>>18);
        index=(index+1)%MT_N;
        return y;
    }
};

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

dtype std_normal(dtype input)
{
    dtype output=(SND_NUM)*exp(-input*input/2);
    return output;
}

dtype MoroInvCND(dtype P){
    const dtype a1 = 2.50662823884;
    const dtype a2 = -18.61500062529;
    const dtype a3 = 41.39119773534;
    const dtype a4 = -25.44106049637;
    const dtype b1 = -8.4735109309;
    const dtype b2 = 23.08336743743;
    const dtype b3 = -21.06224101826;
    const dtype b4 = 3.13082909833;
    const dtype c1 = 0.337475482272615;
    const dtype c2 = 0.976169019091719;
    const dtype c3 = 0.160797971491821;
    const dtype c4 = 2.76438810333863E-02;
    const dtype c5 = 3.8405729373609E-03;
    const dtype c6 = 3.951896511919E-04;
    const dtype c7 = 3.21767881768E-05;
    const dtype c8 = 2.888167364E-07;
    const dtype c9 = 3.960315187E-07;
    dtype y, z;

    y = P - 0.5;  //偏移到中点为原点
    //Moro算法以0.08和0.92为分界点，中间(0.08<x<0.92)用Beasley算法，两边用截断的切比雪夫序列
    if(fabs(y) < 0.42){          
        z = y * y;
        z = y * (((a4 * z + a3) * z + a2) * z + a1) / ((((b4 * z + b3) * z + b2) * z + b1) * z + 1);  //嵌套乘法
    }
    else
    {
        if(y > 0)
            z = log(-log(1.0 - P));
        else
            z = log(-log(P));

        z = c1 + z * (c2 + z * (c3 + z * (c4 + z * (c5 + z * (c6 + z * (c7 + z * (c8 + z * c9)))))));
        if(y < 0) z = -z;
    }

    return z;
}
dtype normalCFD(dtype value)
{
    return 0.5*erfc(-value*M_SQRT2);
}

void least_square_second_order(dtype *data)
{
    
    return;
}
#endif
