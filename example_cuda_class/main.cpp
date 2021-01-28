//2020.12.31
//版本：cuda单个gpu+Option类
//规模：2^20 * 300
//期权类型：亚式期权，看涨
#include<iostream>
#include<cmath>
#include"option.h"
#include<time.h>
using namespace std;
#define N (int)pow(2,20)
#define M 300

int main()
{
    clock_t start=clock();
    unsigned int seed=123456UL;
    Option option(100.0,100.0,1.0,0.05,0.2);
    option.MC_pricing_gpu(N,M,seed);
    cout<<"price: "<<option.price<<endl;
    clock_t end=clock();
 
    dtype time=(end-start)/(dtype)CLOCKS_PER_SEC;
    cout<<"time: "<<time<<" s"<<endl;
    return 0;
}
