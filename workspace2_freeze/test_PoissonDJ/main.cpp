#include<iostream>
#include"util.hpp"
#include<cmath>
#include<random>
#include<time.h>
using namespace std;

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
        //对于 Poisson 跳跃扩散模型：alpha，std，lambda

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

        // * 资产价格生成方法
        //Milstein方法
        inline dtype asset_Milstein_method(const dtype &asset_price)
        {
            dtype rand_num = gaussianBoxMuller();
            //return asset_price*(1+R*dt+V*rand_num*sqrt(dt)+0.5*V*V*(rand_num*rand_num-1)*dt);
            return asset_price*exp((R-0.5*V*V)*dt+V*rand_num*sqrt(dt));
        }

        //Poisson跳跃扩散模型——Test
        inline dtype asset_PoissonJumpDiffusion_method(const dtype &asset_price)
        {
            dtype lambda=var[0]*dt;
            dtype jump_mu=var[1];
            dtype jump_sigma=var[2];
            dtype kappa=exp(jump_mu)-1.0;

            dtype normal_rand = gaussianBoxMuller();
            dtype poisson_rand = poisson(lambda);

            dtype sum_jump_range=0.0;
            for(int i=0;i<poisson_rand;i++)
            {
                sum_jump_range+=jump_mu+jump_sigma*gaussianBoxMuller();
                //sum_jump_range+=(jump_mu-0.5*jump_sigma*jump_sigma)+jump_sigma*gaussianBoxMuller();
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
            for(int iter=0;iter<sim_num;iter++)
            {
                dtype asset_price = S0;
                for(int t=1;t<time_step;t++)
                {
                    asset_price = asset_Milstein_method(asset_price);
                    //asset_price = asset_PoissonJumpDiffusion_method(asset_price);
                }
                dtype payoff;
                if (!is_input)
                    payoff=payoff_call(asset_price);
                else
                    payoff=payoff_put(asset_price);
                payoff_sum+=payoff;
                payoff_array[iter]=payoff;
            }
            price = payoff_sum/static_cast<dtype>(sim_num);
            price *= discount;

            compute_error(payoff_array,sim_num);
            delete payoff_array;
            return;
        }
};
int main(int argc, char **argv)
{
    srand(time(NULL));
    int sim_num=100000;
    int time_step=200;
    if(argc>=2)
    {
        sim_num=atoi(argv[1]);
    }

    dtype S0=100.0,K=100.0,T=1.0,V=0.2,R=0.05;
    Option option(S0,K,T,R,V);

    //lambda,mu,sigma
    option.set_for_PoissonJD(1.0,1.083287,0.4);
    option.MC_pricing(sim_num,time_step,true);
    cout<<"使用 Poisson 跳扩散 方法"<<endl;
    cout<<"看涨期权价格为："<<option.price<<endl;
    cout<< "欧式期权误差为："<<option.error<<endl;
    cout << "欧式期权95%置信区间为 ["<< option.price-1.96*option.error<<" ; "<< option.price+1.96*option.error<<" ]"<< endl;
    cout<<endl;

    return 0;
}


