#include<iostream>
#include<cmath>
#include<cstdlib>
//#include"AsianOption.h"
#include<openacc_curand.h>
#include"Asset.h"
using namespace std;
#define dtype float 
#define N 100000

dtype Asset::rate=0.05;
dtype Asset::volarity=0.2;
int main()
{
    Asset asset(100,1,100);
    dtype MC_mean=0;
    #pragma acc data copyin(Asset::rate,Asset::volarity)
    {
        //#pragma acc parallel loop reduction(+:MC_mean)
        //for(int i=0;i<N;i++)
        //    MC_mean+=asset.MC_simul_unit_gpu();
        MC_mean=asset.MC_simul_gpu(N);

    }
    cout<<MC_mean/N<<endl;

    return 0;
}
