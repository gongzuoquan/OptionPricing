#include<iostream>
#include<cmath>
#include<cstdlib>
#include<chrono>
//#include"AsianOption.h"
#include<openacc_curand.h>
#include"newOption.h"
using namespace std;
using namespace chrono;

void test_option(int simN)
{
    auto start=system_clock::now();
    //dtype S0=100,K=100,T=0.5,V=0.5,R=0.1;
    dtype S0=164,K=166,T=0.0959,V=0.29,R=0.0521;
    int time_step=1000;
    AsianOption option(S0,K,T,R,V);
    //Option option(S0,K,T,R,V);
    cout<<"S0: "<<option.S0<<endl;
    cout<<"K: "<<option.K<<endl;
    cout<<"V: "<<option.V<<endl;
    cout<<"R: "<<option.R<<endl;
    cout<<"T: "<<option.T<<endl;

    cout<<"模拟次数："<<simN<<endl;

    option.MC_pricing(simN,time_step);
    cout<<"看涨期权价格为："<<option.price<<endl;
    option.MC_pricing(simN,time_step,true);
    cout<<"看跌期权价格为："<<option.price<<endl;
    
    dtype cprice=option.BS_return();
    cout<<"BS看涨期权价格："<<cprice<<endl;
    dtype pprice=option.BS_return(true);
    cout<<"BS看跌期权价格："<<pprice<<endl;

    //option.MC_pricing(simN,time_step,"asia");
    //cout<<"看涨期权价格为："<<option.price<<endl;
    //option.MC_pricing(simN,time_step,true,"asia");
    //cout<<"看跌期权价格为："<<option.price<<endl;

    auto end=system_clock::now();
    auto durat=duration_cast<microseconds>(end-start);
    cout<<"total cost "<<double(durat.count())*microseconds::period::num/microseconds::period::den<<endl;

    return;
}
int main(int argc, char **argv)
{
    srand(123456);
    int simN=100000;
    if(argc>=2)
    {
        simN=atoi(argv[1]);
    }
    test_option(simN);

    return 0;
}
