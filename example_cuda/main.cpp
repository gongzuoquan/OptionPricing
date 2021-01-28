//2020.12.13
//版本：cuda单个gpu+helper头文件
//规模：2^20 * 200
//期权类型：亚式期权，看涨
#include<iostream>
#include<cmath>
#include"common.h"
#include<time.h>
using namespace std;
#define N (int)pow(2,20)
//#define N (int)pow(2,16)
#define M 300

extern "C" dtype MC_pricing_gpu(TOption&,const int,const int,unsigned int);
int main()
{
    clock_t start=clock();
    unsigned int seed=123456UL;
    TOption option(100.0,100.0,1.0,0.05,0.2);
    dtype price=MC_pricing_gpu(option,N,M,seed);
    cout<<"price: "<<price<<endl;
    clock_t end=clock();
 
    dtype time=(end-start)/(dtype)CLOCKS_PER_SEC;
    cout<<"time: "<<time<<" s"<<endl;
    return 0;
}
