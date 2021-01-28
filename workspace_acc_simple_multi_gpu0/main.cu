#include<stdio.h>
#include"cuda_runtime.h"
#include"nccl.h"
#include"option.hpp"
#include<iostream>
#define DNUM 2
void test_option(int simN,int time_step)
{
    auto start=system_clock::now();
    //dtype S0=100,K=100,T=0.5,V=0.5,R=0.1;
    //dtype S0=164,K=166,T=0.0959,V=0.29,R=0.0521;

    //dtype S0=50,K=60,T=1,V=0.4,R=0.1;
    dtype S0=100,K=100,T=1.0,V=0.2,R=0.05;

    std::cout<<"S0: "<<S0<<std::endl;
    std::cout<<"K: "<<K<<std::endl;
    std::cout<<"V: "<<V<<std::endl;
    std::cout<<"R: "<<R<<std::endl;
    std::cout<<"T: "<<T<<std::endl;

    std::cout<<"模拟次数："<<simN<<std::endl;
    std::cout<<"时间步数："<<time_step<<std::endl;

    //欧式期权
    Option option(S0,K,T,R,V);
    option.MC_pricing(simN,time_step);
    std::cout<<"看涨期权价格为："<<option.price<<std::endl;
    std::cout<< "欧式期权误差为："<<option.error<<std::endl;
    std::cout << "欧式期权95%置信区间为 ["<< option.price-1.96*option.error<<" ; "<< option.price+1.96*option.error<<" ]"<< std::endl;
    std::cout<<"使用 Euler Maruyama 方法"<<endl;
    std::cout<<std::endl;

    auto end=system_clock::now();
    auto durat=duration_cast<microseconds>(end-start);
    std::cout<<"total cost "<<double(durat.count())*microseconds::period::num/microseconds::period::den<<endl;

    return;
}
int main(int argc, char **argv)
{
    srand(random_seed);
    
    if(argc>=2)
    {
        simN=atoi(argv[1]);
    }
    if(argc>=3)
    {
        time_step=atoi(argv[2]);
    }
    for(int i=0;i<TYPE_NUM;i++)
    {
        test_flag[i]=false;
    }
    ncclComm_t comms[DNUM];
    int devs[DNUM]={0,1};

    test_option(simN,time_step);

    return 0;
}
