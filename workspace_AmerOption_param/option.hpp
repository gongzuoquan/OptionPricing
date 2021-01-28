#ifndef OPTION_CPU_
#define OPTION_CPU_
#include<iostream>
#include"util.hpp"
#include<random>
#include<fstream>
#include"param.h"
//#include<cmath>
//#include<string>
using namespace std;

/*2020.11.18
 * 此份代码为重要组分
 * 包含 
 * 1.NCCL多GPU并行 --NO
 * 2.美式期权LSMC方法 --OK
 * 3.多维资产 --NO
 * 4.对偶变量法 --NO
 * 5.Heston随机波动率模型 --NO
 * 6.cublas进行最小二乘和楚列斯基分解 --NO
*/
//一、期权模板类
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

        string save_path;
        string load_path;
        int num_save;


        //新增参数
        dtype rho; //该参数用于生成对偶变量，为对偶变量的相关系数
        dtype kappa,theta,epslon; //随机波动率方程：dv(t)=kappa(theta-v(t))dt+epslon*sqrt(v(t))dWt

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
            save_path="save_MC_stock_path.csv";
            load_path="save_MC_stock_path.csv";
            num_save=20;

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
        Option( Params params)
        {
            S0=params.S0;
            K=params.K;
            T=params.T;
            R=params.R;
            V=params.V;
            price=0.0;
            discount=exp(-R*T);
            save_path="save_MC_stock_path.csv";
            load_path="save_MC_stock_path.csv";
            num_save=20;

            cout<<"\noption initialization"<<endl;
        }
        ~Option(){}

        void set_S0(dtype in_S0)
        {
            S0=in_S0;
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
        dtype payoff_call(const dtype &St)
        {
            return max(St-K,0);
        }

        //2.看跌期权
        dtype payoff_put(const dtype &St)
        {
            return max(K-St,0);
        }
                
        //保存MC路径
        void save_stock_prices(dtype *stock_price, int num_save, int time_step)
        {
            cout<<num_save<<" "<<time_step<<endl;
            ofstream fsave(save_path, ios::out);
            for(int i=0;i<num_save;i++)
            {
                for(int t=0;t<time_step;t++)
                {
                    fsave<<stock_price[i*time_step+t];
                    if(t!=time_step-1)
                        fsave<<",";
                }
                fsave<<endl;
            }
            fsave.close();
            return ;
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

};

//二、美式期权 
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
        AmericanOption(
                Params params):Option(params)
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
        // Least Square Monte Carlo定价
        dtype MC_pricing(
                int sim_num,
                int time_step,
                dtype rho,
                const bool &is_input=false
                )
        {
            time_step=time_step+1;
            dtype *asset_mat=new dtype[sim_num*time_step]; //存储整个模拟过程中的所有股票价格

            dt = T/static_cast<dtype>(time_step);
            dtype *payoff_mat=new dtype[sim_num*time_step]; //payoff矩阵记录了所有路径所有时间点上的payoff
            dtype *payoff_not_zero_x=new dtype[sim_num]; //记录payoff不为零的点，用于最小二乘拟合
            dtype *payoff_not_zero_y=new dtype[sim_num]; //记录payoff不为零的点，用于最小二乘拟合
            
            //初始化相关系数矩阵
            dtype corr[4];
            corr[0]=1.0,corr[1]=rho,corr[2]=rho,corr[3]=1.0;
            // | 1.0   rho |
            // | rho   1.0 |
            dtype L[4];
            L[0]=0.0,L[1]=0.0,L[2]=0.0,L[3]=0.0;
            cholesky_decomposition(corr,L,2,2);

            //1.价格生成
            for(int iter=0;iter<sim_num;iter++)
            {
                dtype asset_price = S0;
                asset_mat[iter*time_step]=asset_price;
                for(int t=1;t<time_step;t++)
                {
                    dtype rand_num = gaussianBoxMuller();
                    asset_price=asset_price*exp((R-0.5*V*V)*dt+V*rand_num*sqrt(dt));
                    asset_mat[iter*time_step+t]=asset_price; //所有资产价格都进行储存

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
                        payoff_not_zero_y[num_non_zero]=payoff_mat[i*time_step+exec_time[i]]*exp(-R*dt*(exec_time[i]-t));  //** 要注意这里折现时的操作，时间要从执行时间点开始往回算
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
                }//

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

#endif
