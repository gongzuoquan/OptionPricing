#ifndef OPTION_CPU_
#define OPTION_CPU_
#include<iostream>
#include"util.h"
#include<random>
#include<fstream>
//#include<cmath>
//#include<string>
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

        //dtype var[10]; //参数池
        //对于 Poisson 跳跃扩散模型：
        dtype lambda,jump_mu,jump_sigma;

        //随FDM方法更改的变量
        dtype dS;
        dtype Smin;
        dtype Smax;

        string save_path;
        string load_path;
        int num_save;

        dtype (Option::*asset_func)(const dtype&); //资产路径生成函数
        unsigned long long seed; //curand 随机种子，用于并行

    public:
        Option();
        Option( const dtype &in_S0,
                const dtype &in_K,
                const dtype &in_T,
                const dtype in_R,
                const dtype in_V);
        
        Option(Option &in_option);

        ~Option(){}

        void set_S0(dtype in_S0);
        
        // * 资产价格生成方法
        //1.Euler—Maruyama方法
        inline dtype asset_EulerMaruyama_method(const dtype &asset_price)
        {
            dtype rand_num = gaussianBoxMuller();
            return asset_price*(1+R*dt+V*rand_num*sqrt(dt));
        }

        //2.Milstein方法
        inline dtype asset_Milstein_method(const dtype &asset_price)
        {
            dtype rand_num = gaussianBoxMuller();
            //return asset_price*(1+R*dt+V*rand_num*sqrt(dt)+0.5*V*V*(rand_num*rand_num-1)*dt);
            return asset_price*exp((R-0.5*V*V)*dt+V*rand_num*sqrt(dt));
        }

        //3.Poisson跳跃扩散模型
        inline dtype asset_PoissonJumpDiffusion_method(const dtype &asset_price)
        {
            //dtype lambda=var[0]*dt;
            //dtype jump_mu=var[1];
            //dtype jump_sigma=var[2];
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
                );
        
        // 该函数用于切换资产价格生成方法
        void change_asset_func(dtype (Option::*in_func)(const dtype&));

        void display();
                
        //1.看涨期权
        dtype payoff_call(const dtype &St);
        
        //2.看跌期权
        dtype payoff_put(const dtype &St);
        
        //Black-Scholes解析解(call/put)
        dtype BS_return(const bool &is_put=false);
                
        //保存MC路径
        void save_stock_prices(dtype *stock_price, int num_save, int time_step);
                
        //从csv文件读取MC路径
        void load_stock_prices(dtype *stock_price,int &sim_num,int &time_step);
        
        //MC误差计算
        void compute_error(const dtype *payoff_array,int sim_num);
        
        // * Monte Carlo定价
        void MC_pricing(
                const int &sim_num,
                const int &time_step,
                const bool &is_input=false
                );
        
        void show_grid(dtype *grid,int sim_num,int time_step);

        void show_fdm_all(dtype *grid,int sim_num,int time_step);
        
        // * FDM 定价
        void FDM_pricing(
                const int &sim_num,
                const int &time_step,
                const dtype in_Smin, //设置网格边界，价格上涨最小值
                const dtype in_Smax, //设置网格边界，价格上涨最大值
                const bool &is_put=false
                );
};

/*
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
        dtype MC_pricing(
                const int &sim_num,
                const int &time_step,
                const bool &is_input=false,
                const bool &is_fixed=false
                )
        {
            dt = T/static_cast<dtype>(time_step);
            dtype payoff_sum = 0.0;
            dtype *payoff_array = new dtype[sim_num];
            for(int iter=0;iter<sim_num;iter++)
            {
                dtype asset_price = S0;
                dtype sum = 0.0;
                for(int t=1;t<time_step;t++)
                {
                    //asset_price = asset_EulerMaruyama_method(asset_price);
                    asset_price = (this->*asset_func)(asset_price);
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
            price = payoff_sum/static_cast<dtype>(sim_num);
            price *= discount;

            compute_error(payoff_array,sim_num);
            delete payoff_array;
            return price;
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
        dtype MC_pricing(
                const int &sim_num,
                const int &time_step,
                const bool &is_input=false,
                const bool &is_fixed=false
                )
        {
            dt = T/static_cast<dtype>(time_step);
            dtype payoff_sum = 0.0;
            dtype *payoff_array = new dtype[sim_num];
            for(int iter=0;iter<sim_num;iter++)
            {
                dtype asset_price = S0;
                dtype max_price=S0;
                dtype min_price=S0;
                for(int t=1;t<time_step;t++)
                {
                    asset_price = (this->*asset_func)(asset_price);
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
            price = payoff_sum/static_cast<dtype>(sim_num);
            price *= discount;

            compute_error(payoff_array,sim_num);
            delete payoff_array;
            return price;
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
        dtype MC_pricing_in( //敲入类型
                const int &sim_num,
                const int &time_step,
                const bool &is_input=false
                )
        {
            dt = T/static_cast<dtype>(time_step);
            dtype payoff_sum = 0.0;
            dtype *payoff_array = new dtype[sim_num];
            for(int iter=0;iter<sim_num;iter++)
            {
                dtype asset_price = S0;
                bool is_in_barrier_state=false; //是否处于障碍状态，即处于触发资产价格逾越障碍价格的状态
                dtype in_barrier_state_time=0.0; //处于障碍状态的累积时间，一旦该时间大于障碍条件时间，则触发敲入/敲出事件
                bool is_in_entry_state=false; //处于无效状态
                for(int t=1;t<time_step;t++)
                {
                    asset_price = (this->*asset_func)(asset_price);
                    
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
            price = payoff_sum/static_cast<dtype>(sim_num);
            price *= discount;

            compute_error(payoff_array,sim_num);
            delete payoff_array;
            return price;
        }

        //2.巴黎期权敲出类型MC模拟
        dtype MC_pricing_out( //敲出类型
                const int &sim_num,
                const int &time_step,
                const bool &is_input=false
                )
        {
            dt = T/static_cast<dtype>(time_step);
            dtype payoff_sum = 0.0;
            dtype *payoff_array = new dtype[sim_num];
            for(int iter=0;iter<sim_num;iter++)
            {
                dtype asset_price = S0;
                bool is_in_barrier_state=false; //是否处于障碍状态，即处于触发资产价格逾越障碍价格的状态
                dtype in_barrier_state_time=0.0; //处于障碍状态的累积时间，一旦该时间大于障碍条件时间，则触发敲入/敲出事件
                bool is_in_entry_state=true; //处于有效状态
                for(int t=1;t<time_step;t++)
                {
                    asset_price = (this->*asset_func)(asset_price);
                    
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
            price = payoff_sum/static_cast<dtype>(sim_num);
            price *= discount;

            compute_error(payoff_array,sim_num);
            delete payoff_array;
            return price;
        }

};

//五、美式期权 
class AmericanOption:public Option
{
    public:
        AmericanOption(
                const dtype &in_S0,
                const dtype &in_K,
                const dtype &in_T,
                const dtype in_R,
                const dtype in_V):Option(in_S0,in_K,in_T,in_R,in_V)
        {
            cout<<"AmericanOption initialization"<<endl;
        }
        AmericanOption(AmericanOption &in_option)
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
        // * Monte Carlo定价
        // Least Square Monte Carlo定价
        dtype MC_pricing(
                int sim_num,
                int time_step,
                const bool &is_input=false
                )
        {
            time_step=time_step+1;
            dtype *asset_mat=new dtype[sim_num*time_step]; //存储整个模拟过程中的所有股票价格

            dt = T/static_cast<dtype>(time_step);
            dtype *payoff_mat=new dtype[sim_num*time_step]; //payoff矩阵记录了所有路径所有时间点上的payoff
            dtype *payoff_not_zero_x=new dtype[sim_num]; //记录payoff不为零的点，用于最小二乘拟合
            dtype *payoff_not_zero_y=new dtype[sim_num]; //记录payoff不为零的点，用于最小二乘拟合

            //1.价格生成

            // 从文件中加载资产价格路径，默认路径为 load_path="save_MC_stock_path.csv";
            //load_stock_prices(asset_mat,sim_num,time_step);
            //dt = T/static_cast<dtype>(time_step);
            
            for(int iter=0;iter<sim_num;iter++)
            {
                dtype asset_price = S0;
                asset_mat[iter*time_step]=asset_price;
                for(int t=1;t<time_step;t++)
                {
                    //价格生成，如果调用了load这里以下两行需要注释
                    asset_price = (this->*asset_func)(asset_price);
                    asset_mat[iter*time_step+t]=asset_price; //所有资产价格都进行储存

                    //调用load之后执行
                    //asset_price=asset_mat[iter*time_step+t];
                    
                    dtype payoff;
                    if (!is_input)
                        payoff=payoff_call(asset_price);
                    else
                        payoff=payoff_put(asset_price);
                    payoff_mat[iter*time_step+t]=payoff;
                }
                //每条路径的每个时间节点都进行计算
            }

            //保存资产价格路径，默认路径为 save_path="save_MC_stock_path.csv";
            //save_stock_prices(asset_mat,20,time_step);
            
            //这里我们还需要记录行权时刻
            int *exec_time=new int[sim_num];
            for(int e=0;e<sim_num;e++)
            {
                exec_time[e]=time_step-1;
            }

            //2.下面开始反向回溯价格路径，比较每个时刻执行方案与持有方案的收益，以确定决策
            //我们要比较的是当前立即执行所获得的收益（存储在payoff_mat），和持有价值的收益（需要使用拟合）
            //持有价值如何拟合：X——当前股票价格（收益不为零的点），Y——下一时段的payoff值
            //使用该拟合函数生成该时间段的持有收益

            for(int t=time_step-2;t>=0;t--)
            {
                //先要判断payoff是否为零，不为零的提取出来，作为拟合点
                int num_non_zero=0;
                //持有价值如何拟合：X——当前股票价格（收益不为零的点），Y——下一时段的payoff值
                for(int i=0;i<sim_num;i++)
                {
                    if(payoff_mat[i*time_step+t]>0)
                    {
                        payoff_not_zero_x[num_non_zero]=asset_mat[i*time_step+t]; //当前点的股票价格
                        payoff_not_zero_y[num_non_zero]=payoff_mat[i*time_step+exec_time[i]]*exp(-R*dt*(exec_time[i]-t));  // ** 要注意这里折现时的操作，时间要从执行时间点开始往回算
                        num_non_zero++;
                    }
                }

                LeastSquareEigen ls(payoff_not_zero_x,payoff_not_zero_y,num_non_zero); //最小二乘拟合

                //使用该拟合函数生成该时间段的持有收益
                for(int i=0;i<sim_num;i++)
                {
                    if(payoff_mat[i*time_step+t]>0)
                    {
                        dtype continuation_payoff=ls.get_Y(asset_mat[i*time_step+t]);
                        if(continuation_payoff<=payoff_mat[i*time_step+t])
                        {
                            exec_time[i]=t;
                        }
                    }
                }

            }

            //cout<<"执行时刻："<<endl;
            dtype payoff_sum=0.0;
            dtype *payoff_array=new dtype[sim_num];
            for(int i=0;i<sim_num;i++)
            {
                payoff_array[i]=payoff_mat[i*time_step+exec_time[i]]*exp(-R*dt*exec_time[i]);
                payoff_sum+=payoff_array[i];
            }

            price = payoff_sum/static_cast<dtype>(sim_num);

            compute_error(payoff_array,sim_num);
            cout<<"定价结果："<<price<<endl;

            delete[] payoff_array;
            delete[] payoff_not_zero_x;
            delete[] payoff_not_zero_y;

            delete[] asset_mat;
            delete[] payoff_mat;
            delete[] exec_time;

            return price;
        }
};
*/

#endif
