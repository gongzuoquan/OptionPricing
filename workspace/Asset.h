#ifndef ASSET_H_
#define ASSET_H_
#include<cmath>
#include"Util.h"
//#include<openacc.h>
//#include</public/apps/pgi/linux86-64/19.10/include/openacc_curand.h>
#include<openacc_curand.h>

/*
 * 2020.3.20
 * 本文件将创建 Asset 资产类
 * 该类将在主类Option 期权类中使用
 * 含义是期权的标的资产（实际情况下为股票、外汇、利率等）
 * 资产的特性是其价格会随时间变动（例如：股票）
 */

/*
class Asset
{
    public:
        static dtype rate;
        static dtype volarity;
        dtype S0; //资产起始价格
        dtype *S; //记录随时间变化的资产价格
        unsigned int T; //记录时间
        unsigned int time_step; //时间轴上的采样数目
    public:
        explicit Asset()
        {
            S0 = 0.0;
            T = 0;
            time_step=0;
#pragma acc enter data copyin(this)
        }
        Asset(dtype in_S0,unsigned int in_T,unsigned int in_time_step)
        {
            if(in_time_step<=0)
            {
                std::cout<<"输入参数错误，资产类（Asset）对象初始化失败"<<std::endl;
            }
            S0 = in_S0;
            T = in_T;
            time_step=in_time_step;
            //S = new dtype[in_T];
            //S[0] = S0;
#pragma acc enter data copyin(this)
        }

        ~Asset()
        {
            //delete(S);
            return ;
        }


        static void setRate(dtype in_rate)
        {
            rate = in_rate;
            return ;
        }

        static void setVolarity(dtype in_volarity)
        {
            volarity = in_volarity;
            return ;
        }

        // MC算法单元的GPU版本
        //#pragma acc routine seq
        dtype MC_simul_unit_gpu()
        {
            dtype asset_price = S0;
            dtype dt=static_cast<dtype>(T)/static_cast<dtype>(time_step);
            dtype rand_num;
            dtype sum = 0.0;

            unsigned long long seed=1235ULL;
            //unsigned long long seed=rand();
            unsigned long long seq=0ULL;
            unsigned long long offset=0ULL;
            curandState_t state;
            curand_init(seed,seq,offset,&state);

            for(unsigned long long t=0;t<time_step;t++)
            {
                rand_num=curand_normal_double(&state);

                // 欧拉-丸山 Method
                //asset_price = asset_price*(1+rate*dt+volarity*rand_num*sqrt(dt));
                asset_price = asset_price*exp( (rate-0.5*volarity*volarity)*dt+volarity*rand_num*sqrt(dt));
                sum += asset_price;
            }

            return sum/static_cast<dtype>(time_step) ;
        }
};
*/

class Asset
{
    public:
        dtype S0; //资产起始价格
        dtype *S; //记录随时间变化的资产价格
        unsigned int T; //记录时间
        unsigned int time_step; //时间轴上的采样数目
        dtype rate; //利率
        dtype volarity; //波动率
    public:
        explicit Asset()
        {
            S0 = 0.0;
            T = 0;
            time_step=0;
#pragma acc enter data copyin(this)
        }
        Asset(dtype in_S0,unsigned int in_T,unsigned int in_time_step)
        {
            if(in_time_step<=0)
            {
                std::cout<<"输入参数错误，资产类（Asset）对象初始化失败"<<std::endl;
            }
#pragma acc enter data copyin(this)
            S0 = in_S0;
#pragma acc update device(S0)
            T = in_T;
#pragma acc update device(T)
            time_step=in_time_step;
#pragma acc update device(time_step)
            //S = new dtype[in_T];
            //S[0] = S0;
        }

        ~Asset()
        {
            //delete(S);
#pragma acc exit data delete(this)
            return ;
        }

        dtype MC_simul_unit_cpu()
        {
            // 生成随机数链，正态分布
            // 一共有sim_num条模拟路径，
            // 每条路径理论上是连续时间上的随机波动曲线，但我们将其处理为离散的采样点
            // 通常假设中这条链是连续时间马尔科夫过程，即当前价格依赖于之前的价格的相关公式
            // 数学上处理为随机微分方程，数值分析中使用欧拉-丸山方法可以进行处理
            // 根据Black-Scholes公式，模拟出整条链上的数据
            // 注：B-S公式的使用需要提供无风险利率rate（r），波动率volatility（v）
            dtype asset_price = S0;
            dtype dt=static_cast<dtype>(T)/static_cast<dtype>(time_step);
            dtype rand_num;
            dtype sum = 0.0;

            for(int t=0;t<time_step;t++)
            {
                rand_num = gaussianBoxMuller();
                // 欧拉-丸山 Method
                //asset_price = asset_price*(1+rate*dt+volarity*rand_num*sqrt(dt));
                asset_price = asset_price*exp( (rate-0.5*volarity*volarity)*dt+volarity*rand_num*sqrt(dt));
                sum += asset_price;
            }

            //std::cout<<sum/static_cast<dtype>(time_step)<<std::endl;
            return sum/static_cast<dtype>(time_step) ;
        }

        /*******************************************************************/


        // MC算法单元的GPU版本
        //#pragma acc routine seq
        dtype MC_simul_gpu(const int n)
        {
            
//            dtype rate=0.05;
//            dtype volarity=0.2;
            dtype asset_price = S0;
            dtype dt=static_cast<dtype>(T)/static_cast<dtype>(time_step);
            dtype rand_num;
            dtype sum = 0.0;

            unsigned long long seed=12345ULL;
            //unsigned long long seed=rand();
            unsigned long long seq=0ULL;
            unsigned long long offset=0ULL;
            #pragma acc parallel loop 
            for(unsigned int iter=0;iter<n;iter++)
            {
                curandState_t state;
                curand_init(seed,seq,offset,&state);

                #pragma acc loop seq
                for(unsigned long long t=0;t<time_step;t++)
                {
                    rand_num=curand_normal_double(&state);

                    // 欧拉-丸山 Method
                    //asset_price = asset_price*(1+rate*dt+volarity*rand_num*sqrt(dt));
                    asset_price = asset_price*exp( (rate-0.5*volarity*volarity)*dt+volarity*rand_num*sqrt(dt));
                    sum += asset_price;
                }
            }

            return sum/static_cast<dtype>(time_step) ;
        }
};

#endif 
