//2020.12.28
//测试TRNG+openMP
//欧式看涨
#include<iostream>
#include<cstdio>
#include<cmath>
#include<omp.h>
#include<trng/yarn2.hpp>
//#include<trng/uniform01_dist.hpp> //0到1间的均匀分布
#include<trng/normal_dist.hpp> //正态分布
#include<chrono>
using namespace std;
using namespace chrono;
#define dtype double

int main()
{
    auto start=system_clock::now();
    dtype S0=100.0;
    dtype K=100.0;
    dtype R=0.05;
    dtype V=0.2;
    dtype T=1.0;

    int sim_num=pow(2,20);
    //int sim_num=100;
    //const int sim_num=1000001; //需要加上const，否则编译不通过
    const int time_step=300;
    dtype dt=T/time_step;
    dtype sum=0.0;
    const unsigned long seed = 12345UL;
    
    //omp_set_num_threads(40);
    #pragma omp parallel reduction(+:sum) default(none) shared(sim_num,time_step,S0,K,T,V,R,dt)
    {
        trng::yarn2 rng; //随机数生成器(引擎)
        rng.seed(seed);

        int size = omp_get_num_threads(); //线程总数
        int rank=omp_get_thread_num(); //当前线程id
        //printf("总线程数：%d\n",size);

        trng::normal_dist<> normal(0,1);      //正态分布
        rng.jump(2*(rank*sim_num/size)); //jump ahead

        // throw random points into square
        for(long i=rank*sim_num/size;i<(rank+1)*sim_num/size;++i)
        {
            dtype St=S0;
            for(int t=0;t<time_step;t++)
            {
                St=St*exp((R-0.5*V*V)*dt+sqrt(dt)*V*normal(rng));
            }
            dtype payoff=(St-K>0)?(St-K):0.0;
            //printf("%.3f\n",payoff);
            //cout<<payoff<<endl;
            sum+=payoff;
        }
    }
    dtype average=sum/sim_num;
    dtype price=average*exp(-R*T);
    cout<<"average: "<<average<<endl;
    cout<<price<<endl;

    auto end=system_clock::now();
    auto durat=duration_cast<microseconds>(end-start);
    cout<<"程序耗时："<<(double)(durat.count())*microseconds::period::num/microseconds::period::den<<endl;
    return 0;
}
