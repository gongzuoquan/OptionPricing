//2020.12.27
//版本：acc+omp多个gpu+nvprof
//规模：2^20 * 300
//期权类型：欧式看涨
//标准结果：10.45
#include<iostream>
#include<cmath>
#include<ctime>
#include<openacc_curand.h>
#include<openacc.h>
#include<omp.h>
#include<chrono>
using namespace std;
using namespace chrono;
#define dtype float

int sim_num=pow(2,20);
int time_step=300;
int main(int argc,char **argv)
{
    //clock_t start=clock();
    auto start=system_clock::now();
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
    dtype *payoff_array=new dtype[sim_num];

    int gpu_n=acc_get_num_devices(acc_device_nvidia);
    printf("Number of GPUs = %d\n",gpu_n);
    if (gpu_n<=0)
    {
        printf("error: no nvidia.\n");
        exit(1);
    }

    #pragma omp parallel num_threads(gpu_n)
    {
        acc_set_device_num(omp_get_thread_num()+1,acc_device_nvidia);
        printf("GPU Device #%d\n",acc_get_current_cuda_device());

        int gpu_num=acc_get_num_devices(acc_device_nvidia);
        int gpu_id=acc_get_current_cuda_device()-1;
        printf("gpu id: %d\n",gpu_id);

        int len=sim_num/gpu_num;
        int rem=sim_num%gpu_num;
        int start,end;
        //负载均衡，综合考虑数据局部性:
        if(gpu_id<rem)
        {
            start=gpu_id*(len+1);
            end=start+len;
        }
        else
        {
            start=gpu_id*len+rem;
            end=start+len-1;
        }
        cout<<start<<" "<<end<<endl;

        #pragma acc enter data create(payoff_array[start:end-start+1])
        {
            #pragma acc parallel
            {
                curandState_t state;
                #pragma acc loop reduction(+:sum) private(state) independent
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
            }//end acc parallel    
        }
        acc_wait_all();
        #pragma acc exit data copyout(payoff_array[start:end-start+1])

    }//end omp parallel

    dtype average=sum/sim_num;
    dtype price=average*exp(-R*T);
    cout<<"该期权定价结果为："<<price<<endl;

    dtype var=0.0;
    for(int i=0;i<sim_num;i++)
    {
        var+=(payoff_array[i]-average)*(payoff_array[i]-average);
    }

    var/=sim_num;
    dtype error=sqrt(var*exp(-2*R*T)/sim_num);
    cout<<"95%置信区间为："<<"["<<price-error<<","<<price+error<<"]"<<endl;

    //clock_t end=clock();
    //cout<<"程序耗时："<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
    auto end = system_clock::now();
    auto durat=duration_cast<microseconds>(end-start);
    cout<<"程序耗时："<<(double)(durat.count())*microseconds::period::num/microseconds::period::den<<endl;
    delete payoff_array;

    return 0;
}
