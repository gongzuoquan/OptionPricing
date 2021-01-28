#ifndef OPTION_CPU_
#define OPTION_CPU_
#include<iostream>
#include<omp.h>
#include<openacc.h>
#include"util.hpp"
#include<cmath>
#include<string>
#include<random>
#include<openacc_curand.h>
#include<trng/yarn2.hpp>
#include<trng/normal_dist.hpp>
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

        //随MC方法更改的变量
        dtype error; //误差
        dtype dt; //时间微元
        dtype var[10]; //参数池
        //1.对于 Poisson 跳跃扩散模型：alpha，std，lambda
        //2.

        dtype (Option::*asset_func)(const dtype&,const dtype&); //资产路径生成函数
        unsigned long long seed;

    public:
        Option(){
            S0=0.0,K=0.0,T=0.0,R=0.0,V=0.0;
            price=0.0;
            discount=0.0;
        }
        Option( const dtype &in_S0,
                const dtype &in_K,
                const dtype &in_T,
                const dtype in_R,
                const dtype in_V):S0(in_S0),K(in_K),T(in_T),R(in_R),V(in_V)
        {
            price=0.0;
            discount=exp(-R*T);
            asset_func=&Option::asset_EulerMaruyama_method;
            seed=123ULL;

            cout<<"\noption initialization"<<endl;
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

        // * 资产价格生成方法
        //1.Euler—Maruyama方法
        //#pragma acc routine 
        inline dtype asset_EulerMaruyama_method(const dtype &asset_price,const dtype &random_num)
        {
            //dtype rand_num = gaussianBoxMuller();
            return asset_price*(1+R*dt+V*random_num*sqrt(dt));
        }

        //2.Milstein方法
        //#pragma acc routine 
        inline dtype asset_Milstein_method(const dtype &asset_price,const dtype &random_num)
        {
            return asset_price*(1+R*dt+V*random_num*sqrt(dt)+0.5*V*V*(random_num*random_num-1)*dt);
            //return asset_price*exp((R-0.5*V*V)*dt+V*rand_num*sqrt(dt));
        }

        /*
        //3.Poisson跳跃扩散模型——Test
        inline dtype asset_PoissonJumpDiffusion_method(const dtype &asset_price)
        {
            dtype lambda=var[0]; //Poisson分布参数lambda
            dtype jump_mu=var[1]; //跳跃幅度所服从的正态分布的均值
            dtype jump_sigma=var[2]; //跳跃幅度所服从的正态分布的均值
            dtype kappa=exp(jump_mu)-1.0;

            dtype normal_rand = gaussianBoxMuller();
            dtype poisson_rand = poisson(lambda);

            dtype sum_jump_range=0.0;
            for(int i=0;i<poisson_rand;i++)
            {
                sum_jump_range+=(jump_mu-0.5*jump_sigma*jump_sigma)+jump_sigma*gaussianBoxMuller();
            }

            return asset_price*exp((R-0.5*V*V-lambda*kappa)*dt+V*normal_rand*sqrt(dt)+sum_jump_range);
        }
        // 在使用Poisson跳跃扩散模型前需要对模型参数进行初始化
        int set_for_PoissonJD(         
                const dtype lambda, 
                const dtype jump_mu, 
                const dtype jump_sigma
                )
        {
            var[0]=lambda;
            var[1]=jump_mu;
            var[2]=jump_sigma;
            return 0;
        }
        */

        // 该函数用于切换资产价格生成方法
        void change_asset_func(dtype (Option::*in_func)(const dtype&,const dtype&))
        {
            asset_func=in_func;
            return;
        }

        void display()
        {
            cout<<"S0: "<<S0<<endl;
            cout<<"K: "<<K<<endl;
            cout<<"V: "<<V<<endl;
            cout<<"R: "<<R<<endl;
            cout<<"T: "<<T<<endl;
        }
        
        //1.看涨期权
        #pragma acc routine 
        dtype payoff_call(const dtype &St)
        {
            return max(St-K,0);
        }

        //2.看跌期权
        #pragma acc routine 
        dtype payoff_put(const dtype &St)
        {
            return max(K-St,0);
        }

        //Black-Scholes解析解(call/put)
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

        //MC误差计算
        void compute_error(const dtype *payoff_array,int sim_num)
        {
            dtype variance=0.0;
            for(int iter=0;iter<sim_num;iter++)
            {
                variance+=(payoff_array[iter]-price)*(payoff_array[iter]-price);
            }
            variance/=sim_num;
            error=sqrt(variance*exp(-2*R*T)/sim_num);

            return;
        }

        /*
        //多GPU定价部分：
        // * Monte Carlo定价
        void MC_pricing(
                const int &sim_num,
                const int &time_step,
                const bool &is_input=false
                )
        {
            if(time_step<=0 || sim_num<=0)
            {
                cerr<<"MC_sim Function got error parameters."<<endl;
                cerr<<"TIME_STEP and SIM_NUM must be positive."<<endl;
                return;
            }
            //unsigned long long seq=0ULL;
            //unsigned long long offset=0ULL;
            long long seq=0;
            long long offset=0;
            dt = T/static_cast<dtype>(time_step);
            dtype payoff_sum = 0.0;
            dtype *payoff_array = new dtype[sim_num];
            
            int gpu_n=acc_get_num_devices(acc_device_nvidia);
            printf("Number of GPUs = %d\n",gpu_n);
            if (gpu_n<=0) 
            {
                printf("error: no nvidia.\n");
                exit(1);
            }

            #pragma omp parallel num_threads(gpu_n)
            {
                acc_set_device_num(omp_get_thread_num()+1, acc_device_nvidia);
                printf("GPU Device #%d\n",acc_get_current_cuda_device());

                int gpu_id=acc_get_current_cuda_device()-1;
                int gpu_num=acc_get_num_devices(acc_device_nvidia);

                int len=sim_num/gpu_num;
                int rem=sim_num%gpu_num;

                int start,end;
                if(gpu_id<rem)
                {
                    start=gpu_id*(len+1);
                    end=start+len;
                }
                else{
                    start=gpu_id*len+rem;
                    end=start+len-1;
                }
                cout<<start<<" "<<end<<endl;

                #pragma acc enter data create(payoff_array[start:end-start+1])
                {
                    #pragma acc parallel
                    {
                        curandState_t state;
                        #pragma acc loop reduction(+:payoff_sum) private(state)
                        for(unsigned long long iter=start;iter<end;iter++)
                        {
                            dtype asset_price = S0;
                            seed+=iter;
                            curand_init(seed,seq,offset,&state);

                            #pragma acc loop seq
                            for(int t=1;t<time_step;t++)
                            {
                                //asset_price = asset_EulerMaruyama_method(asset_price);
                                double randn=curand_normal_double(&state);
                                asset_price=asset_price*exp((R-0.5*V*V)*dt+V*randn*sqrt(dt));
                                //asset_price = (this->*asset_func)(asset_price,randn);
                            }
                            dtype payoff;
                            if (!is_input)
                                payoff=payoff_call(asset_price);
                            else
                                payoff=payoff_put(asset_price);
                            payoff_sum+=payoff;
                            payoff_array[iter]=payoff;
                        }
                    }
                }
                acc_wait_all();
                #pragma acc exit data copyout(payoff_array[start:end-start+1])
            }
            
            price = payoff_sum/static_cast<dtype>(sim_num);
            price *= discount;

            compute_error(payoff_array,sim_num);
            delete payoff_array;
            return;
        }
        */
/*
        // * Monte Carlo定价
        void MC_pricing(
                const int &sim_num,
                const int &time_step,
                const bool &is_input=false
                )
        {
            if(time_step<=0 || sim_num<=0)
            {
                cerr<<"MC_sim Function got error parameters."<<endl;
                cerr<<"TIME_STEP and SIM_NUM must be positive."<<endl;
                return;
            }
            //unsigned long long seq=0ULL;
            //unsigned long long offset=0ULL;
            long long seq=0;
            long long offset=0;
            dt = T/static_cast<dtype>(time_step);
            dtype payoff_sum = 0.0;
            dtype restrict *payoff_array = new dtype[sim_num];

            #pragma acc parallel copyout(payoff_array[0:sim_num])
            {
                curandState_t state;
                #pragma acc loop reduction(+:payoff_sum) private(state)
                for(unsigned long long iter=0;iter<sim_num;iter++)
                {
                    dtype asset_price = S0;
                    seed+=iter;
                    curand_init(seed,seq,offset,&state);

                    #pragma acc loop seq
                    for(int t=1;t<time_step;t++)
                    {
                        //asset_price = asset_EulerMaruyama_method(asset_price);
                        double randn=curand_normal_double(&state);
                        asset_price=asset_price*exp((R-0.5*V*V)*dt+V*randn*sqrt(dt));
                        //asset_price = (this->*asset_func)(asset_price,randn);
                    }
                    dtype payoff;
                    if (!is_input)
                        payoff=payoff_call(asset_price);
                    else
                        payoff=payoff_put(asset_price);
                    payoff_sum+=payoff;
                    payoff_array[iter]=payoff;
                }
            }
            
            price = payoff_sum/static_cast<dtype>(sim_num);
            price *= discount;

            compute_error(payoff_array,sim_num);
            delete payoff_array;
            return;
        }
*/
        // * Monte Carlo定价
        void MC_pricing(
                const int &sim_num,
                const int &time_step,
                const bool &is_input=false
                )
        {
            if(time_step<=0 || sim_num<=0)
            {
                cerr<<"MC_sim Function got error parameters."<<endl;
                cerr<<"TIME_STEP and SIM_NUM must be positive."<<endl;
                return;
            }
            dt = T/static_cast<dtype>(time_step);
            dtype payoff_sum = 0.0;
            dtype *payoff_array = new dtype[sim_num];
            unsigned long seed = 12345UL; //<--add

            #pragma omp parallel reduction(+:payoff_sum) default(none) shared(sim_num,time_step,S0,K,T,V,R,dt,seed,is_input,payoff_array) //<--add
            {
                trng::yarn2 rng; //<--add
                rng.seed(seed); //<--add
                int size=omp_get_num_threads(); //<--add
                int rank=omp_get_thread_num(); //<--add

                trng::normal_dist<> normal(0,1); //<--add
                rng.jump(2*(rank*sim_num/size)); //<--add

                int len=sim_num/size;
                int rem=sim_num%size;
                int start,end;
                if(rank<rem)
                {
                    start=rank*(len+1);
                    end=start+len;
                }
                else{
                    start=rank*len+rem;
                    end=start+len;
                }
                //cout<<start<<" : "<<end<<endl;
                //printf("范围：%d ~ %d\n",start,end);

                for(int  iter=start;iter<end;iter++) //<--add
                {
                    dtype asset_price = S0;
                    for(int t=1;t<time_step;t++)
                    {
                        //asset_price = asset_EulerMaruyama_method(asset_price);
                        //double randn=curand_normal_double(&state);
                        double randn=normal(rng);
                        asset_price=asset_price*exp((R-0.5*V*V)*dt+V*randn*sqrt(dt));
                        //asset_price = (this->*asset_func)(asset_price,randn);
                    }
                    dtype payoff;
                    if (!is_input)
                        payoff=payoff_call(asset_price);
                    else
                        payoff=payoff_put(asset_price);
                    payoff_sum+=payoff;
                    payoff_array[iter]=payoff;
                }
            }
            
            price = payoff_sum/static_cast<dtype>(sim_num);
            price *= discount;

            compute_error(payoff_array,sim_num);
            delete payoff_array;
            return;
        }
};

//二、亚式期权
class AsianOption:public Option
{
    public:
        AsianOption(
                const dtype &in_S0,
                const dtype &in_K,
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
            price=in_option.price; // *期权价格
            discount=in_option.discount; //折现因子
        }
    public:
        //1.浮动看涨亚式期权
        dtype payoff_call_float(const dtype &St,const dtype &average)
        {
            return max(St-average,0);
        }

        //2.浮动看跌亚式期权
        dtype payoff_put_float(const dtype &St,const dtype &average)
        {
            return max(average-St,0);
        }

        //3.固定看涨亚式期权
        dtype payoff_call_fixed(const dtype &average)
        {
            return max(average-K,0);
        }

        //4.固定看跌亚式期权
        dtype payoff_put_fixed(const dtype &average)
        {
            return max(K-average,0);
        }

        // * Monte Carlo定价
        void MC_pricing(
                const int &sim_num,
                const int &time_step,
                const bool &is_input=false,
                const bool &is_fixed=false
                )
        {
            if(time_step<=0 || sim_num<=0)
            {
                cerr<<"MC_sim Function got error parameters."<<endl;
                cerr<<"TIME_STEP and SIM_NUM must be positive."<<endl;
                return;
            }
            dt = T/static_cast<dtype>(time_step);
            dtype payoff_sum = 0.0;
            dtype *payoff_array = new dtype[sim_num];
            unsigned long seed =12345UL;

            #pragma omp parallel reduction(+:payoff_sum) default(none) shared(sim_num,time_step,S0,K,T,V,R,dt,seed,is_input,payoff_array,is_fixed)
            {
                trng::yarn2 rng;
                rng.seed(seed);
                int size=omp_get_num_threads();
                int rank=omp_get_thread_num();

                trng::normal_dist<> normal(0,1);
                rng.jump(2*(rank*sim_num/size));

                int len=sim_num/size;
                int rem=sim_num%size;
                int start,end;
                if(rank<rem)
                {
                    start=rank*(len+1);
                    end=start+len;
                }
                else{
                    start=rank*len+rem;
                    end=start+len;
                }

                for(int iter=start;iter<end;iter++)
                {
                    dtype asset_price = S0;
                    dtype sum = 0.0;
                    for(int t=1;t<time_step;t++)
                    {
                        //asset_price = asset_EulerMaruyama_method(asset_price);
                        double randn=normal(rng);
                        asset_price=asset_price*exp((R-0.5*V*V)*dt+V*randn*sqrt(dt));
                        //asset_price = (this->*asset_func)(asset_price);
                        sum += asset_price;
                    }
                    //expect_sum+=sum/static_cast<dtype>(time_step);
                    dtype average=sum/static_cast<dtype>(time_step);
                    dtype payoff;
                    if (!is_input)
                    {
                        if(!is_fixed)
                            payoff=payoff_call_float(asset_price,average);
                        else
                            payoff=payoff_call_fixed(average);
                    }
                    else
                     {
                        if(!is_fixed)
                            payoff=payoff_put_float(asset_price,average);
                        else
                            payoff=payoff_put_fixed(average);
                    }
                    //cout<<payoff<<endl;
                    payoff_sum+=payoff;
                    payoff_array[iter]=payoff;
                }

            }
                price = payoff_sum/static_cast<dtype>(sim_num);
                price *= discount;

                compute_error(payoff_array,sim_num);
                delete payoff_array;
                return;
            }

};

//三、回望期权
class LookbackOption:public Option
{
    public:
        LookbackOption(
                const dtype &in_S0,
                const dtype &in_K,
                const dtype &in_T,
                const dtype in_R,
                const dtype in_V):Option(in_S0,in_K,in_T,in_R,in_V)
        {
            cout<<"LookbackOption initialization"<<endl;
        }
        LookbackOption(LookbackOption &in_option)
        {
            S0=in_option.S0; //开盘价
            K=in_option.K; //敲定价
            T=in_option.T; //成熟期
            R=in_option.R; //无风险利率
            V=in_option.V; //波动率
            price=in_option.price; // *期权价格
            discount=in_option.discount; //折现因子
        }
    public:
        //1.浮动看涨回望期权
        dtype payoff_call_float(const dtype &St,const dtype &Smin)
        {
            return max(St-Smin,0);
        }

        //2.浮动看跌回望期权
        dtype payoff_put_float(const dtype &St,const dtype &Smax)
        {
            return max(Smax-St,0);
        }

        //3.固定看涨回望期权
        dtype payoff_call_fixed(const dtype &Smax)
        {
            return max(Smax-K,0);
        }

        //4.固定看跌回望期权
        dtype payoff_put_fixed(const dtype &Smin)
        {
            return max(K-Smin,0);
        }

        // * Monte Carlo定价
        void MC_pricing(
                const int &sim_num,
                const int &time_step,
                const bool &is_input=false,
                const bool &is_fixed=false
                )
        {
            if(time_step<=0 || sim_num<=0)
            {
                cerr<<"MC_sim Function got error parameters."<<endl;
                cerr<<"TIME_STEP and SIM_NUM must be positive."<<endl;
                return;
            }
            dt = T/static_cast<dtype>(time_step);
            dtype payoff_sum = 0.0;
            dtype *payoff_array = new dtype[sim_num];
            unsigned long seed = 12345UL;

            #pragma omp parallel reduction(+:payoff_sum) default(none) shared(sim_num,time_step,S0,K,T,V,R,dt,seed,is_input,payoff_array,is_fixed)
            {
                trng::yarn2 rng;
                rng.seed(seed);
                int size=omp_get_num_threads();
                int rank=omp_get_thread_num();

                trng::normal_dist<> normal(0,1);
                rng.jump(2*(rank*sim_num/size));

                int len=sim_num/size;
                int rem=sim_num%size;
                int start,end;
                if(rank<rem)
                {
                    start=rank*(len+1);
                    end=start+len;
                }
                else{
                    start=rank*len+rem;
                    end=start+len;
                }

                for(int iter=start;iter<end;iter++)
                {
                    dtype asset_price = S0;
                    dtype max_price=S0;
                    dtype min_price=S0;
                    for(int t=1;t<time_step;t++)
                    {
                        //asset_price = (this->*asset_func)(asset_price);
                        double randn=normal(rng);
                        asset_price=asset_price*exp((R-0.5*V*V)*dt+V*randn*sqrt(dt));
                        max_price=max(asset_price,max_price);
                        min_price=min(asset_price,min_price);
                    }

                    dtype payoff;
                    if (!is_input)
                    {
                        if(!is_fixed)//float price
                            payoff=payoff_call_float(asset_price,min_price);
                        else //fixed price
                            payoff=payoff_call_fixed(max_price);
                    }
                    else
                     {
                        if(!is_fixed)//float price
                            payoff=payoff_put_float(asset_price,max_price);
                        else //fixed price
                            payoff=payoff_put_fixed(min_price);
                    }
                    //cout<<payoff<<endl;
                    payoff_sum+=payoff;
                    payoff_array[iter]=payoff;
                }
             
            }
                price = payoff_sum/static_cast<dtype>(sim_num);
                price *= discount;

                compute_error(payoff_array,sim_num);
                delete payoff_array;
                return;
            }
};

//四、巴黎期权（欧式）
class ParisianOption:public Option
{
    public:
        dtype barrier_price; //触发障碍状态的价格
        dtype barrier_condition_time; //障碍形成条件
        bool is_up_type; //是否为向上类型
    public:
        ParisianOption(
                const dtype &in_S0,
                const dtype &in_K,
                const dtype &in_T,
                const dtype &in_R,
                const dtype &in_V,
                const dtype &in_barrier_price, //触发障碍状态的价格
                const dtype &in_barrier_condition_time, //障碍形成条件
                const bool &in_is_up_type=false //是否为向上类型
                ):Option(in_S0,in_K,in_T,in_R,in_V)
        {
            barrier_price=in_barrier_price; //触发障碍状态的价格
            barrier_condition_time=in_barrier_condition_time; //障碍形成条件
            is_up_type=in_is_up_type; //是否为向上类型
            cout<<"ParisianOption initialization"<<endl;
        }
        ParisianOption(ParisianOption &in_option)
        {
            S0=in_option.S0; //开盘价
            K=in_option.K; //敲定价
            T=in_option.T; //成熟期
            R=in_option.R; //无风险利率
            V=in_option.V; //波动率
            price=in_option.price; // *期权价格
            discount=in_option.discount; //折现因子

            barrier_price=in_option.barrier_price; //触发障碍状态的价格
            barrier_condition_time=in_option.barrier_condition_time; //障碍形成条件
            is_up_type=in_option.is_up_type; //是否为向上类型
        }
    public:
        // * Monte Carlo定价
        // 巴黎期权分为 向上/向下-敲入/敲出-看涨/看跌

        //1.巴黎期权敲入类型MC模拟
        void MC_pricing_in( //敲入类型
                const int &sim_num,
                const int &time_step,
                const bool &is_input=false
                )
        {
            if(time_step<=0 || sim_num<=0)
            {
                cerr<<"MC_sim Function got error parameters."<<endl;
                cerr<<"TIME_STEP and SIM_NUM must be positive."<<endl;
                return;
            }
            dt = T/static_cast<dtype>(time_step);
            dtype payoff_sum = 0.0;
            dtype *payoff_array = new dtype[sim_num];
            unsigned long seed = 12345UL;

            #pragma omp parallel reduction(+:payoff_sum) default(none) shared(sim_num,time_step,S0,K,T,V,R,dt,seed,is_input,payoff_array,barrier_price,is_up_type,barrier_condition_time)
            {
                trng::yarn2 rng;
                rng.seed(seed);
                int size=omp_get_num_threads();
                int rank=omp_get_thread_num();

                trng::normal_dist<> normal(0,1);
                rng.jump(2*(rank*sim_num/size));

                int len=sim_num/size;
                int rem=sim_num%size;
                int start,end;
                if(rank<rem)
                {
                    start=rank*(len+1);
                    end=start+len;
                }
                else{
                    start=rank*len+rem;
                    end=start+len;
                }

                for(int iter=start;iter<end;iter++)
                {
                    dtype asset_price = S0;
                    bool is_in_barrier_state=false; //是否处于障碍状态，即处于触发资产价格逾越障碍价格的状态
                    dtype in_barrier_state_time=0.0; //处于障碍状态的累积时间，一旦该时间大于障碍条件时间，则触发敲入/敲出事件
                    bool is_in_entry_state=false; //处于无效状态
                    for(int t=1;t<time_step;t++)
                    {
                        double randn=normal(rng);
                        asset_price=asset_price*exp((R-0.5*V*V)*dt+V*randn*sqrt(dt));
                        
                        if(!is_in_entry_state) //如果当前尚未期权敲入
                        {
                            if((asset_price>barrier_price && is_up_type)||(asset_price<barrier_price && !is_up_type))//触发障碍条件
                            {
                                if(!is_in_barrier_state)
                                {
                                    is_in_barrier_state=true;
                                }
                                in_barrier_state_time+=dt;
                                if(in_barrier_state_time>barrier_condition_time)
                                {
                                    is_in_entry_state=true;
                                }
                            }
                            else
                            {
                                if(is_in_barrier_state) //如果未进入敲入状态，而退出了障碍状态，则障碍累计时间归零
                                {
                                    is_in_barrier_state=false;
                                    in_barrier_state_time=0.0;
                                }
                            }
                        }
                    }
                    
                    dtype payoff=0.0;
                    if(is_in_entry_state)
                    {
                        if (!is_input)
                        {
                            payoff=payoff_call(asset_price);
                        }
                        else
                        {
                            payoff=payoff_put(asset_price);
                        }
                    }
                    
                    //cout<<payoff<<endl;
                    payoff_sum+=payoff;
                    payoff_array[iter]=payoff;
                }

            }
                price = payoff_sum/static_cast<dtype>(sim_num);
                price *= discount;

                compute_error(payoff_array,sim_num);
                delete payoff_array;
                return;
            }

        //2.巴黎期权敲出类型MC模拟
        void MC_pricing_out( //敲出类型
                const int &sim_num,
                const int &time_step,
                const bool &is_input=false
                )
        {
            if(time_step<=0 || sim_num<=0)
            {
                cerr<<"MC_sim Function got error parameters."<<endl;
                cerr<<"TIME_STEP and SIM_NUM must be positive."<<endl;
                return;
            }
            dt = T/static_cast<dtype>(time_step);
            dtype payoff_sum = 0.0;
            dtype *payoff_array = new dtype[sim_num];
            unsigned long seed = 12345UL;

            #pragma omp parallel reduction(+:payoff_sum) default(none) shared(sim_num,time_step,S0,K,T,V,R,dt,seed,is_input,payoff_array,barrier_price,is_up_type,barrier_condition_time)
            {
                trng::yarn2 rng;
                rng.seed(seed);
                int size=omp_get_num_threads();
                int rank=omp_get_thread_num();

                trng::normal_dist<> normal(0,1);
                rng.jump(2*(rank*sim_num/size));

                int len=sim_num/size;
                int rem=sim_num%size;
                int start,end;
                if(rank<rem)
                {
                    start=rank*(len+1);
                    end=start+len;
                }
                else{
                    start=rank*len+rem;
                    end=start+len;
                }

                for(int iter=start;iter<end;iter++)
                {
                    dtype asset_price = S0;
                    bool is_in_barrier_state=false; //是否处于障碍状态，即处于触发资产价格逾越障碍价格的状态
                    dtype in_barrier_state_time=0.0; //处于障碍状态的累积时间，一旦该时间大于障碍条件时间，则触发敲入/敲出事件
                    bool is_in_entry_state=true; //处于有效状态
                    for(int t=1;t<time_step;t++)
                    {
                        double randn=normal(rng);
                        asset_price=asset_price*exp((R-0.5*V*V)*dt+V*randn*sqrt(dt));
                        
                        if(is_in_entry_state) //如果当前尚未期权敲出
                        {
                            if((asset_price>barrier_price && is_up_type)||(asset_price<barrier_price && !is_up_type))//触发障碍条件
                            {
                                if(!is_in_barrier_state)
                                {
                                    is_in_barrier_state=true;
                                }
                                in_barrier_state_time+=dt;
                                if(in_barrier_state_time>barrier_condition_time)
                                {
                                    is_in_entry_state=false;
                                }
                            }
                            else
                            {
                                if(is_in_barrier_state) //如果未进入敲出状态，而退出了障碍状态，则障碍累计时间归零
                                {
                                    is_in_barrier_state=false;
                                    in_barrier_state_time=0.0;
                                }
                            }
                        }
                    }
                    
                    dtype payoff=0.0;
                    if(is_in_entry_state)
                    {
                        if (!is_input)
                        {
                            payoff=payoff_call(asset_price);
                        }
                        else
                        {
                            payoff=payoff_put(asset_price);
                        }
                    }
                    
                    //cout<<payoff<<endl;
                    payoff_sum+=payoff;
                    payoff_array[iter]=payoff;
                }

            }
                price = payoff_sum/static_cast<dtype>(sim_num);
                price *= discount;

                compute_error(payoff_array,sim_num);
                delete payoff_array;
                return;
            }

};

#endif
