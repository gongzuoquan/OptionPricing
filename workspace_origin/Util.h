#ifndef UTIL_H_
#define UTIL_H_
#include<cstdlib>

#define dtype float
//工具类，用于记录一些实用函数
dtype gaussianBoxMuller()
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
#endif
