#include<iostream>
#include<cstdlib>
#include"tree.hpp"
#include<cmath>
#include<cstdlib>
#include<chrono>
#include<openacc_curand.h>
#include<unistd.h>
using namespace std;
using namespace chrono;
#define TYPE_NUM 6
enum option_type{euro,asian,lookup,parisian,american,bermuda};
bool test_flag[TYPE_NUM];
int simN=1048576;
int time_step=300;
int random_seed=1234;

void test_option(int simN,int time_step)
{
    //test_flag[euro]=true;
    test_flag[asian]=true;
    //test_flag[lookup]=true;
    //test_flag[parisian]=true;

    auto start=system_clock::now();
    //dtype S0=100,K=100,T=0.5,V=0.5,R=0.1;
    //dtype S0=164,K=166,T=0.0959,V=0.29,R=0.0521;

    //dtype S0=50,K=60,T=1,V=0.4,R=0.1;
    dtype S0=100,K=100,T=1.0,V=0.2,R=0.05;

    cout<<"S0: "<<S0<<endl;
    cout<<"K: "<<K<<endl;
    cout<<"V: "<<V<<endl;
    cout<<"R: "<<R<<endl;
    cout<<"T: "<<T<<endl;

    cout<<"模拟次数："<<simN<<endl;
    cout<<"时间步数："<<time_step<<endl;

    auto option_start=system_clock::now();
    auto option_end=system_clock::now();
    auto option_durat=duration_cast<microseconds>(option_end-option_start);

    //欧式期权
    if(test_flag[euro])
    {
        Option option(S0,K,T,R,V);
        //warm up
        option.MC_pricing(simN,time_step);
        option_start=system_clock::now();
        option.MC_pricing(simN,time_step);
        cout<<">欧式看涨期权价格为："<<option.price<<endl;
        cout<< ">欧式看涨期权误差为："<<option.error<<endl;
        cout << ">欧式看涨期权95%置信区间为 [ "<< option.price-1.96*option.error<<" ; "<< option.price+1.96*option.error<<" ]"<< endl;
        //cout<<"使用 Euler Maruyama 方法"<<endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">欧式看涨期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;

        option_start=system_clock::now();
        option.MC_pricing(simN,time_step,true);
        cout<<">欧式看跌期权价格为："<<option.price<<endl;
        cout<< ">欧式看跌期权误差为："<<option.error<<endl;
        cout << ">欧式看跌期权95%置信区间为 [ "<< option.price-1.96*option.error<<" ; "<< option.price+1.96*option.error<<" ]"<<endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">欧式看跌期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den;
        cout<<endl;    

        /*
        // Black-Scholes定价
        dtype cprice=option.BS_return();
        cout<<"BS看涨期权价格："<<cprice<<endl;
        dtype pprice=option.BS_return(true);
        cout<<"BS看跌期权价格："<<pprice<<endl;
        cout<<endl;
        */
    }

    //亚式期权
    if(test_flag[asian])
    {
        AsianOption aoption(S0,K,T,R,V);
        option_start=system_clock::now();
        aoption.MC_pricing(simN,time_step);
        cout<<">浮动亚式看涨期权价格为："<<aoption.price<<endl;
        cout<< ">浮动亚式看涨期权误差为："<<aoption.error<<endl;
        cout << ">浮动亚式看涨期权95%置信区间为 [ "<< aoption.price-1.96*aoption.error<<" ; "<< aoption.price+1.96*aoption.error<<" ]"<< endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">浮动亚式看涨期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;

        option_start=system_clock::now();
        aoption.MC_pricing(simN,time_step,true);
        cout<<">浮动亚式看跌期权价格为："<<aoption.price<<endl;
        cout<< ">浮动亚式看跌期权误差为："<<aoption.error<<endl;
        cout << ">浮动亚式看跌期权95%置信区间为 [ "<< aoption.price-1.96*aoption.error<<" ; "<< aoption.price+1.96*aoption.error<<" ]"<< endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">浮动亚式看跌期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;

        option_start=system_clock::now();
        aoption.MC_pricing(simN,time_step,false,true);
        cout<<">固定亚式看涨期权价格为："<<aoption.price<<endl;
        cout<< ">固定亚式看涨期权误差为："<<aoption.error<<endl;
        cout << ">固定亚式看涨期权95%置信区间为 [ "<< aoption.price-1.96*aoption.error<<" ; "<< aoption.price+1.96*aoption.error<<" ]"<< endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">固定亚式看涨期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;

        option_start=system_clock::now();
        aoption.MC_pricing(simN,time_step,true,true);
        cout<<">固定亚式看跌期权价格为："<<aoption.price<<endl;
        cout<< ">固定亚式看跌期权误差为："<<aoption.error<<endl;
        cout << ">固定亚式看跌期权95%置信区间为 [ "<< aoption.price-1.96*aoption.error<<" ; "<< aoption.price+1.96*aoption.error<<" ]"<< endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">固定亚式看跌期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;
    }

    //回望期权
    if(test_flag[lookup])
    {
        LookbackOption loption(S0,K,T,R,V);

        option_start=system_clock::now();
        loption.MC_pricing(simN,time_step);
        cout<<">浮动价格回望看涨期权价格为："<<loption.price<<endl;
        cout<< ">浮动价格回望看涨期权误差为："<<loption.error<<endl;
        cout << ">浮动价格回望看涨期权95%置信区间为 [ "<< loption.price-1.96*loption.error<<" ; "<< loption.price+1.96*loption.error<<" ]"<< endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">浮动价格看涨期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;

        option_start=system_clock::now();
        loption.MC_pricing(simN,time_step,true);
        cout<<">浮动价格回望看跌期权价格为："<<loption.price<<endl;
        cout<< ">浮动价格回望看跌期权误差为："<<loption.error<<endl;
        cout << ">浮动价格回望看跌期权95%置信区间为 [ "<< loption.price-1.96*loption.error<<" ; "<< loption.price+1.96*loption.error<<" ]"<< endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">浮动价格回望看跌期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;
        
        option_start=system_clock::now();
        loption.MC_pricing(simN,time_step,false,true);
        cout<<">固定价格回望看涨期权价格为："<<loption.price<<endl;
        cout<< ">固定价格回望看涨期权误差为："<<loption.error<<endl;
        cout << ">固定价格回望看涨期权95%置信区间为 [ "<< loption.price-1.96*loption.error<<" ; "<< loption.price+1.96*loption.error<<" ]"<< endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">固定价格回望看涨期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;

        option_start=system_clock::now();
        loption.MC_pricing(simN,time_step,true,true);
        cout<<">固定价格回望看跌期权价格为："<<loption.price<<endl;
        cout<< ">固定价格回望看跌期权误差为："<<loption.error<<endl;
        cout << ">固定价格回望看跌期权95%置信区间为 [ "<< loption.price-1.96*loption.error<<" ; "<< loption.price+1.96*loption.error<<" ]"<< endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">固定价格回望看跌期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;

    }
    
    //欧式巴黎期权
    if(test_flag[parisian])
    {
        //巴黎期权分两个MC定价函数：
        //结尾为in的是敲入期权；结尾为out的是敲出期权
        //格式为：
        //MC_pricing_in/out(模拟次数-int, 时间步数-int, 障碍价格-dtype, 障碍窗口时间-dtype, 是否为看跌-bool, 是否为向上类型-bool)
        dtype pS0=12.0,pK=10.0,pT=1.0,pV=0.2,pR=0.05;
        dtype barrier_price=12.0;
        dtype barrier_condition_time=0.1;

        //dtype pS0=16.0,pK=8.0,pT=1.0,pV=0.2,pR=0.05;
        //dtype barrier_price=17.0;
        //dtype barrier_condition_time=0.1;

        ParisianOption poption(pS0,pK,pT,pR,pV,barrier_price,barrier_condition_time,true);
        //ParisianOption poption(S0,K,T,R,V);
        poption.display();
        cout<<endl;

        option_start=system_clock::now();
        poption.MC_pricing_in(simN,time_step);
        cout<<">向上敲入巴黎看涨期权价格为："<<poption.price<<endl;
        cout<< ">向上敲入巴黎看涨期权误差为："<<poption.error<<endl;
        cout<<  ">向上敲入巴黎看涨期权95%置信区间为 [ "<< poption.price-1.96*poption.error<<" ; "<< poption.price+1.96*poption.error<<" ]"<<endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">向上敲入巴黎看涨期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;

        //poption.MC_pricing_in(simN,time_step,true);
        //cout<<"向上敲入巴黎看跌期权价格为："<<poption.price<<endl;
        //cout<< "向上敲入巴黎看跌期权误差为："<<poption.error<<endl;
        //cout<<  "向上敲入巴黎看跌期权95%置信区间为 [ "<< poption.price-1.96*poption.error<<" ; "<< poption.price+1.96*poption.error<<" ]"<<endl;
        //cout<<endl;

        option_start=system_clock::now();
        poption.MC_pricing_out(simN,time_step);
        cout<<">向上敲出巴黎看涨期权价格为："<<poption.price<<endl;
        cout<< ">向上敲出巴黎看涨期权误差为："<<poption.error<<endl;
        cout << ">向上敲出巴黎看涨期权95%置信区间为 [ "<< poption.price-1.96*poption.error<<" ; "<< poption.price+1.96*poption.error<<" ]"<< endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">向上敲出巴黎看涨期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;

        //poption.MC_pricing_out(simN,time_step,true);
        //cout<<"向上敲出巴黎看跌期权价格为："<<poption.price<<endl;
        //cout<< "向上敲出巴黎看跌期权误差为："<<poption.error<<endl;
        //cout << "向上敲出巴黎看跌期权95%置信区间为 [ "<< poption.price-1.96*poption.error<<" ; "<< poption.price+1.96*poption.error<<" ]"<< endl;
        //cout<<endl;
        
    }


    if(test_flag[parisian])
    {
        dtype pS0=12.0,pK=10.0,pT=1.0,pV=0.2,pR=0.05;
        dtype barrier_price=12.0;
        dtype barrier_condition_time=0.1;

        //dtype pS0=16.0,pK=8.0,pT=1.0,pV=0.2,pR=0.05;
        //dtype barrier_price=17.0;
        //dtype barrier_condition_time=0.1;

        ParisianOption poption(pS0,pK,pT,pR,pV,barrier_price,barrier_condition_time,false);
        poption.display();
        cout<<endl;

        option_start=system_clock::now();
        poption.MC_pricing_in(simN,time_step);
        cout<<">向下敲入巴黎看涨期权价格为："<<poption.price<<endl;
        cout<< ">向下敲入巴黎看涨期权误差为："<<poption.error<<endl;
        cout<<  ">向下敲入巴黎看涨期权95%置信区间为 [ "<< poption.price-1.96*poption.error<<" ; "<< poption.price+1.96*poption.error<<" ]"<< endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">向下敲入巴黎看涨期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;

        //poption.MC_pricing_in(simN,time_step,true);
        //cout<<"向下敲入巴黎看跌期权价格为："<<poption.price<<endl;
        //cout<< "向下敲入巴黎看跌期权误差为："<<poption.error<<endl;
        //cout<<  "向下敲入巴黎看跌期权95%置信区间为 ["<< poption.price-1.96*poption.error<<" ; "<< poption.price+1.96*poption.error<<" ]"<< endl;
        //cout<<endl;
        
        option_start=system_clock::now();
        poption.MC_pricing_out(simN,time_step);
        cout<<">向下敲出巴黎看涨期权价格为："<<poption.price<<endl;
        cout<< ">向下敲出巴黎看涨期权误差为："<<poption.error<<endl;
        cout<<   ">向下敲出巴黎看涨期权95%置信区间为 [ "<< poption.price-1.96*poption.error<<" ; "<< poption.price+1.96*poption.error<<" ]"<< endl;
        option_end=system_clock::now();
        option_durat=duration_cast<microseconds>(option_end-option_start);
        cout<<">向下敲出巴黎看涨期权运行时间： "<<double(option_durat.count())*microseconds::period::num/microseconds::period::den<<endl;
        cout<<endl;

        //poption.MC_pricing_out(simN,time_step,true);
        //cout<<"向下敲出巴黎看跌期权价格为："<<poption.price<<endl;
        //cout<< "向下敲出巴黎看跌期权误差为："<<poption.error<<endl;
        //cout<<   "向下敲出巴黎看跌期权95%置信区间为 [ "<< poption.price-1.96*poption.error<<" ; "<< poption.price+1.96*poption.error<<" ]"<< endl;
        //cout<<endl;
    }

    //auto end=system_clock::now();
    //auto durat=duration_cast<microseconds>(end-start);
    //cout<<"total cost "<<double(durat.count())*microseconds::period::num/microseconds::period::den<<endl;

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
    test_option(simN,time_step);

    return 0;
}
/*
int main(int argc, char **argv)
{
    int sim_num=100000;
    if (argc>=2)
        sim_num=atoi(argv[1]);
    cout<<"simulation number: "<<sim_num<<endl;
    cout<<endl;
    
    //int time_step=100;
    //if (argc>=2)
    //    time_step=atoi(argv[1]);
    //cout<<"time step: "<<time_step<<endl;
    //cout<<endl;

    ////BitreeOption bt(50,52,0.5,0.05,0,1.2,0.8);
    //BitreeOption bt(50,50,0.5,0.05,0.3,1.2,0.8);
    ////BitreeOption bt(20,21,1.5,0.12,0,1.1,0.9);
    //bt.pricing(time_step,true);
    //cout<<"the price is "<<bt.price<<endl;
    ////bt.replication_sim(4.19,0.62817);

    Option op(50,50,2,.05,.3);
    op.MC_pricing(sim_num,200,true);
    cout<<"the EuropeanOption price is "<<op.price<<endl;
    cout<<endl;

    AsianOption asop(50,50,2,.05,.3);
    asop.MC_pricing(sim_num,200,true);
    cout<<"the AsianOption price is "<<asop.price<<endl;
    cout<<endl;

    LookbackOption lbop(50,50,2,.05,.3);
    lbop.MC_pricing(sim_num,200,true);
    cout<<"the LookbackOption price is "<<lbop.price<<endl;
    cout<<endl;

    return 0;
}*/
