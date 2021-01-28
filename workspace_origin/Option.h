#ifndef OPTION_H_
#define OPTION_H_

/*
 * 2020.3.20
 * 在Prof.Deng的指导下 ，创建Computional Finance Code Library
 * 此版为创始，致力于亚式期权的MC定价算法
 * Option 期权类为此版的核心类
 * 本文件将创建 Option 期权类
 * 定义期权的基本参数（如：起始价格 S0，协议周期 T ，等等）
 */

#include"Asset.h"

class Payoff
{
    public:
        dtype K;
    public:
        Payoff(dtype in_K)
        {
            K=in_K;
            return;
        }
        Payoff()
        {
            K=0.0;
            return;
        }

        dtype getPayoff(dtype in_ST)
        {
            dtype payoff = in_ST-K;
            if(payoff<=0.0)
                return 0.0;
            else
                return payoff;
        }
};

class Option
{
    public:
        Asset asset;
        unsigned int T;
        dtype S0;
        dtype price;
    public:
        Option(const dtype &in_S0,
                const unsigned int &in_T,
                const unsigned int &in_time_step
                )
        {
            asset = Asset(in_S0,in_T,in_time_step);
            S0 = in_S0;
            T = in_T; // 期权周期通常以月计
            price = 0.0;
        }

        Option(const Asset &in_asset)
        {
            asset=in_asset;
            price = 0.0;
        }
        ~Option()
        {
            //delete(asset);
            return;
        }

        
        template<class P>
        void pricing_MC_simulation(const unsigned int sim_num,
                P payoff,
                bool use_gpu=true
                )
        {
            unsigned int iter;
            //dtype *payoff_array = new dtype(sim_num);
            dtype *payoff_array = (dtype*)malloc(sizeof(dtype)*sim_num);
            dtype summary = 0.0;
           
            for(iter=0;iter<sim_num;iter++)
            {
                dtype ST;
                if(use_gpu)
                    ST=asset.MC_simul_unit(); //这里我们希望设计一个Asset的函数，用于生成asset的一次随机模拟
                else
                    ST=asset.MC_simul_unit_gpu(); //这里我们希望设计一个Asset的函数，用于生成asset的一次随机模拟

                //payoff_array[iter] = getEuroPayoff(ST);
                payoff_array[iter] = payoff.getPayoff(ST);
                summary += payoff_array[iter];
            }
            
            price = summary/static_cast<dtype>(sim_num);
            //delete(payoff_array);
            free(payoff_array);
            return;
        }
};
#endif
