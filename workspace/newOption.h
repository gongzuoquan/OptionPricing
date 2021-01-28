#ifndef OPTION_H_
#define OPTION_H_
#include<iostream>
#include"Util.h"
#include<cmath>
#include<string>
using namespace std;

/*2020.7.8
 *我们将建立一个Option类，
 *作为所有期权的父类，所有香草期权、奇异期权都继承于该类
 *(EuroOption,AmerOption,AsianOption)
*/
//一、欧式期权
class Option
{
    public:
        dtype S0; //开盘价
        dtype K; //敲定价
        dtype T; //成熟期
        dtype R; //无风险利率
        dtype V; //波动率
        dtype price; //*期权价格
        dtype discount; //折现因子

    public:
        Option(){
            S0=0.0,K=0.0,T=0.0,R=0.0,V=0.0;
            price=0.0;
            discount=0.0;
        }
        Option(const dtype &in_S0,
                dtype &in_K,
                const dtype &in_T,
                const dtype in_R,
                const dtype in_V):S0(in_S0),K(in_K),T(in_T),R(in_R),V(in_V)
        {
            price=0.0;
            discount=exp(-R*T);
            cout<<"option initialization"<<endl;
        }
        Option(Option &in_option)
        {
            S0=in_option.S0; //开盘价
            K=in_option.K; //敲定价
            T=in_option.T; //成熟期
            R=in_option.R; //无风险利率
            V=in_option.V; //波动率
            price=in_option.price; //*期权价格
            discount=in_option.discount; //折现因子
        }

        ~Option(){}
        
        dtype payoff_call(const dtype &St)
        {
            if(St > K)
                return St-K;
            else
                return 0;

        }

        dtype payoff_put(const dtype &St)
        {
            if(K > St)
                return K-St;
            else
                return 0;

        }

        dtype BS_return(const bool &is_put=false)
        {
            dtype d1=(log(S0/K)+(R+V*V/2)*T)/(V*sqrt(T));
            //dtype d2=log(S0/K)+(R-V*V/2)*T/(V*sqrt(T));
            dtype d2=d1-V*sqrt(T);
            dtype C=0;

            dtype Nd1=0.5*erfc(-d1*M_SQRT1_2);
            dtype Nd2=0.5*erfc(-d2*M_SQRT1_2);

            C=(!is_put)?(S0*Nd1-K*discount*Nd2):(K*discount*(1-Nd2)-S0*(1-Nd1));
            return C;
        }

        void MC_pricing(
                const int &sim_num,
                const dtype &time_step,
                const bool is_input=false
                )
        {
            if(time_step<=0 || sim_num<=0)
            {
                cerr<<"MC_sim Function got error parameters."<<endl;
                cerr<<"TIME_STEP and SIM_NUM must be positive."<<endl;
                return;
            }
            dtype asset_price;
            dtype dt = T/static_cast<dtype>(time_step);
            dtype rand_num;
            //dtype sum ;
            //dtype expect_sum = 0.0;
            dtype payoff_sum = 0.0;
            for(int iter=0;iter<sim_num;iter++)
            {
                asset_price = S0;
                //sum = 0.0;
                for(int t=1;t<time_step;t++)
                {
                    rand_num = gaussianBoxMuller();
                    asset_price = asset_price*exp((R-0.5*V*V)*dt+V*rand_num*sqrt(dt));

                }
                dtype payoff;
                if (!is_input)
                    payoff=payoff_call(asset_price);
                else
                    payoff=payoff_put(asset_price);
                payoff_sum+=payoff;
            }
            price = payoff_sum/static_cast<dtype>(sim_num);
            price *= discount;

            return;
        }
};

//亚式期权
class AsianOption:public Option
{
    public:
        AsianOption(const dtype &in_S0,
                dtype &in_K,
                const dtype &in_T,
                const dtype in_R,
                const dtype in_V):Option(in_S0,in_K,in_T,in_R,in_V)
        {
            cout<<"AsianOption initialization"<<endl;
        }
        AsianOption(AsianOption &in_option)
        {
            S0=in_option.S0; //开盘价
            K=in_option.K; //敲定价
            T=in_option.T; //成熟期
            R=in_option.R; //无风险利率
            V=in_option.V; //波动率
            price=in_option.price; //*期权价格
            discount=in_option.discount; //折现因子
        }
    public:
        dtype MC_pricing(
                const int &sim_num,
                const dtype &time_step,
                const bool is_input=false
                )
        {
            dtype asset_price;
            dtype dt = T/static_cast<dtype>(time_step);
            dtype rand_num;
            dtype sum ;
            //dtype expect_sum = 0.0;
            dtype payoff_sum = 0.0;
            for(int iter=0;iter<sim_num;iter++)
            {
                asset_price = S0;
                sum = 0.0;
                for(int t=1;t<time_step;t++)
                {
                    rand_num = gaussianBoxMuller();
                    asset_price = asset_price*exp((R-0.5*V*V)*dt+V*rand_num*sqrt(dt));
                    sum += asset_price;
                }
                //expect_sum+=sum/static_cast<dtype>(time_step);
                dtype average=sum/static_cast<dtype>(time_step);
                dtype payoff;
                if (!is_input)
                    payoff=payoff_call(average);
                else
                    payoff=payoff_put(average);
                //cout<<payoff<<endl;
                payoff_sum+=payoff;
            }
            //price = expect_sum/static_cast<dtype>(sim_num);
            //cout<<"dt: "<<dt<<endl;
            price = payoff_sum/static_cast<dtype>(sim_num);
            price *= discount;

            return price;
        }
};

//美式期权
class AmericanOption:public Option
{
    public:
        AmericanOption(const dtype &in_S0,
                dtype &in_K,
                const dtype &in_T,
                const dtype in_R,
                const dtype in_V):Option(in_S0,in_K,in_T,in_R,in_V)
        {
            cout<<"AmericanOption initialization"<<endl;
        }

    public:
        dtype MC_pricing(
                const int &sim_num,
                const dtype &time_step,
                const bool is_input=false
                )
        {
            dtype asset_price;
            dtype dt = T/static_cast<dtype>(time_step);
            dtype rand_num;
            dtype sum ;
            //dtype expect_sum = 0.0;
            dtype payoff_sum = 0.0;
            for(int iter=0;iter<sim_num;iter++)
            {
                asset_price = S0;
                sum = 0.0;
                for(int t=1;t<time_step;t++)
                {
                    rand_num = gaussianBoxMuller();
                    asset_price = asset_price*exp((R-0.5*V*V)*dt+V*rand_num*sqrt(dt));
                    sum += asset_price;
                }
                //expect_sum+=sum/static_cast<dtype>(time_step);
                dtype average=sum/static_cast<dtype>(time_step);
                dtype payoff;
                if (!is_input)
                    payoff=payoff_call(average);
                else
                    payoff=payoff_put(average);
                //cout<<payoff<<endl;
                payoff_sum+=payoff;
            }
            //price = expect_sum/static_cast<dtype>(sim_num);
            //cout<<"dt: "<<dt<<endl;
            price = payoff_sum/static_cast<dtype>(sim_num);
            price *= discount;

            return price;
        }
};
#endif
