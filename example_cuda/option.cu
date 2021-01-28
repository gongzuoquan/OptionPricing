#include<curand_kernel.h>
#include<cuda_runtime.h>
#include"helper_cuda.h"
#include"common.h"
#include<stdio.h>
/*2020.12.23
  经过使用nvprof进行分析，发现大部分的时间消耗在init_rand函数上
  因此我们需要对init_rand函数进行调整
 */

__constant__ int time_step;
__constant__ dtype S0;
__constant__ dtype K;
__constant__ dtype T;
__constant__ dtype R;
__constant__ dtype V;
__constant__ dtype dt;

__global__ void init_rand(curandState *const states,const unsigned int seed)
{
    unsigned int tid=blockIdx.x*blockDim.x+threadIdx.x;
    curand_init(seed,tid,0,&states[tid]);
}

__global__ void MC_call(curandState *states,dtype *const array)
{
    unsigned int tid=blockIdx.x*blockDim.x+threadIdx.x;
    curandState local_state=states[tid];
    dtype asset_price=S0;
    for(int i=0;i<time_step;i++)
    {
        asset_price=asset_price*exp((R-0.5*V*V)*dt+V*sqrt(dt)*curand_normal(&local_state));
    }
    //dtype payoff=(asset_price-K)>0.0?(asset_price-K):0.0;
    dtype payoff=(K-asset_price)>-9.0?(K-asset_price):0.0;
    array[tid]=payoff;

    return ;
}

extern "C" dtype MC_pricing_gpu(TOption &option,const int n,const int m,unsigned int seed)
{
    option.display();
    dtype h_dt=option.T/m;
    cudaMemcpyToSymbol( S0,&option.S0,1*sizeof(int),0,cudaMemcpyHostToDevice )  ;
    cudaMemcpyToSymbol( K ,&option.K ,1*sizeof(int),0,cudaMemcpyHostToDevice )  ;
    cudaMemcpyToSymbol( T ,&option.T ,1*sizeof(int),0,cudaMemcpyHostToDevice )  ;
    cudaMemcpyToSymbol( R ,&option.R ,1*sizeof(int),0,cudaMemcpyHostToDevice )  ;
    cudaMemcpyToSymbol( V ,&option.V ,1*sizeof(int),0,cudaMemcpyHostToDevice )  ;
    cudaMemcpyToSymbol( dt,&h_dt,1*sizeof(int),0,cudaMemcpyHostToDevice )  ;

    printf("shape: %d x %d\n",n,m);
    dtype *gpu_array;
    curandState *states;
    checkCudaErrors(cudaMalloc((void**)&states,sizeof(curandState)*n));
    checkCudaErrors(cudaMalloc((void**)&gpu_array,sizeof(dtype)*n));
   
    size_t blocksize=1024;
    dim3 block,grid;
    block.x=blocksize;
    grid.x=(n+blocksize-1)/blocksize;

    //unsigned int seed=12345UL;
    cudaMemcpyToSymbol( time_step,&m,1*sizeof(int),0,cudaMemcpyHostToDevice )  ;
    init_rand<<<grid,block>>>(states,seed);
    MC_call<<<grid,block>>>(states,gpu_array);

    dtype *h_array=(dtype*)malloc(sizeof(dtype)*n);
    cudaMemcpy(h_array,gpu_array,n*sizeof(dtype),cudaMemcpyDeviceToHost);
    dtype sum=0.0;
    for(int i=0;i<n;i++)
        sum+=h_array[i];
    dtype price = sum/n*exp(-option.R*option.T);


    cudaFree(gpu_array);
    cudaFree(states);
    free(h_array);

    return price;
}
