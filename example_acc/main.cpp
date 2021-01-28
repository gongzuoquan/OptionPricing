//2020.12.27
//版本：acc单个gpu+nvprof
//规模：2^20 * 300
//期权类型：欧式看涨
//标准结果：10.45
#include<iostream>
#include<cmath>
#include<ctime>
#include<openacc_curand.h>
using namespace std;
#define dtype float
int sim_num=pow(2,20);
int time_step=300;
int main(int argc,char **argv)
{
    clock_t start=clock();
    dtype S0=100.0;
    dtype K=100.0;
    dtype R=0.05;
    dtype V=0.2;
    dtype T=1.0;
    dtype dt=T/time_step;
    dtype sum=0.0;

    long long seed=12345ULL;
    long long seq=0ULL;
    long long offset=0ULL;

    #pragma acc parallel
    {
        curandState_t state;
        #pragma acc loop reduction(+:sum) private(state)
        for(int i=0;i<sim_num;i++)
        {
            dtype St=S0;
            seed=i;
            curand_init(seed,seq,offset,&state);
            #pragma acc loop seq
            for(int t=0;t<time_step;t++)
            {
                dtype randn=curand_normal(&state);
                St=St*exp((R-0.5*V*V)*dt+V*sqrt(dt)*randn);
            }
            dtype payoff;
            if(St-K>0)
                payoff=St-K;
            else
                payoff=0.0;
            sum+=payoff;
        }
    }

    dtype average=sum/sim_num;
    dtype price=average*exp(-R*T);
    cout<<"该期权定价结果为："<<price<<endl;

    clock_t end=clock();
    cout<<"程序耗时："<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
    return 0;
}
